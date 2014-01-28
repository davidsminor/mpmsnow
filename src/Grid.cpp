#include "Grid.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>

#include <GL/gl.h>

#include <Eigen/Geometry>
#include <Eigen/SVD>

using namespace Eigen;

#define GRAVITY -9.8f
#define COULOMBFRICTION 0.5f
#define DECLARE_WEIGHTARRAY( NAME ) float buf_##NAME[12]; float * NAME[] = { &buf_##NAME[1], &buf_##NAME[5], &buf_##NAME[9] };

Grid::Grid( const ParticleData& d, float gridH, float timeStep, const ConstituativeModel& model )
	: m_gridH( gridH ), m_timeStep( timeStep ), m_constituativeModel( model )
{
	// work out the physical size of the grid:
	m_xmin = m_ymin = m_zmin = 1.e10;
	m_xmax = m_ymax = m_zmax = -1.e10;

	for( size_t i=0; i < d.particleX.size(); ++i )
	{
		minMax( d.particleX[i][0], m_xmin, m_xmax );
		minMax( d.particleX[i][1], m_ymin, m_ymax );
		minMax( d.particleX[i][2], m_zmin, m_zmax );
	}
	
	// calculate grid dimensions and quantize bounding box:
	m_nx = fixDim( m_xmin, m_xmax );
	m_ny = fixDim( m_ymin, m_ymax );
	m_nz = fixDim( m_zmin, m_zmax );
			
	// little array with indexes going from -1 to store shape function weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	Vector3i particleCell;
	
	// calculate masses:
	m_gridMasses.resize( m_nx * m_ny * m_nz );
	m_gridMasses.setZero();

	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		Vector3f x = d.particleX[p];
		cellAndWeights( x, particleCell, w );
		// splat the masses:
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					float weight = w[0][i] * w[1][j] * w[2][k];
					m_gridMasses[ idx ] +=
						d.particleM[p] * weight;
				}
			}
		}
	}
	
	// grid masses can end up less than zero due to numerical issues in the shape functions, so clamp 'em:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		if( m_gridMasses[i] < 0 )
		{
			m_gridMasses[i] = 0;
		}
	}

	// calculate velocities:
	m_gridVelocities.resize( m_nx * m_ny * m_nz * 3 );
	m_gridVelocities.setZero();

	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w );
		
		// splat the velocities:
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					if( m_gridMasses[idx] > 0 )
					{
						float particleMass = d.particleM[p];
						float gridCellMass = m_gridMasses[idx];
						float overallWeight = w[0][i] * w[1][j] * w[2][k] *
							( particleMass / gridCellMass );

						m_gridVelocities.segment<3>( 3 * idx ) += overallWeight * d.particleV[p];
					}
				}
			}
		}
	}		
}


void Grid::draw() const
{
	glColor3f( 0,0.3f,0 );
	glBegin( GL_LINES );

	// xy
	for( int i=0; i <= m_nx; ++i )
	{
		for( int j=0; j <= m_ny; ++j )
		{
			glVertex3f( m_xmin + i * m_gridH, m_ymin + j * m_gridH, m_zmin );
			glVertex3f( m_xmin + i * m_gridH, m_ymin + j * m_gridH, m_zmax );
		}
	}
	// zy
	for( int i=0; i <= m_nz; ++i )
	{
		for( int j=0; j <= m_ny; ++j )
		{
			glVertex3f( m_xmin, m_ymin + j * m_gridH, m_zmin + i * m_gridH );
			glVertex3f( m_xmax, m_ymin + j * m_gridH, m_zmin + i * m_gridH );
		}
	}

	// xz
	for( int i=0; i <= m_nx; ++i )
	{
		for( int j=0; j <= m_nz; ++j )
		{
			glVertex3f( m_xmin + i * m_gridH, m_ymin, m_zmin + j * m_gridH );
			glVertex3f( m_xmin + i * m_gridH, m_ymax, m_zmin + j * m_gridH );
		}
	}
	glEnd();
}

void Grid::computeDensities( ParticleData& d ) const
{
	d.particleDensities.resize( d.particleX.size(), 0 );
	
	// little array with indexes going from -1 to store shape function weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	Vector3i particleCell;

	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w );
		
		// transfer densities back onto the particles:
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					
					// accumulate the particle's density:
					d.particleDensities[p] += w[0][i] * w[1][j] * w[2][k] * m_gridMasses[ idx ] / ( m_gridH * m_gridH * m_gridH );
					
				}
			}
		}
	}
}

void Grid::updateParticleVelocities( ParticleData& d )
{
	// little array with indexes going from -1 to store shape function weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	Vector3i particleCell;

	// this is, like, totally doing things FLIP style. The paper recommends a combination of FLIP and PIC...
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w );

		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					d.particleV[p] += w[0][i] * w[1][j] * w[2][k] *
						( m_gridVelocities.segment<3>( 3 * idx ) - m_prevGridVelocities.segment<3>( 3 * idx ) );
				}
			}
		}
	}
}

void Grid::applyImplicitUpdateMatrix(
	const ParticleData& d,
	const std::vector<CollisionObject*>& collisionObjects,
	const VectorXf& vNPlusOne,
	VectorXf& result ) const
{

	// This method computes the forward momenta in this frame in terms of the velocities
	// in the next frame:
	// m * v^(n+1) - m_timeStep * dF(v^(n+1) * m_timeStep)
	
	// work out force differentials when you perturb the grid positions by v * m_timeStep:
	VectorXf df( vNPlusOne.size() );

	// so: this method effectively applies a symmetric matrix which is quite diagonally
	// dominant, with the diagonals largely controlled by the masses. Unfortunately, if
	// the masses vary significantly form cell to cell (which they almost always do),
	// this makes the eigenvalues SUCK, in that they range from something like 1.e-9 to 20.
	// The conjugate gradient solver hates this, so instead we transform the problem so that
	// those masses move off the main diagonal, but the matrix remains symmetric. This means
	// we need to divide both the input and the output of this function by the square roots
	// of the masses:
	
	// divide input:
	VectorXf vTransformed = vNPlusOne;
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		if( m_gridMasses[i] != 0 )
		{
			vTransformed.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
	}

	calculateForceDifferentials( d, m_timeStep * vTransformed, df );
	
	result.resize( vTransformed.size() );
	for( int i=0; i < m_nx; ++i )
	{
		for( int j=0; j < m_ny; ++j )
		{
			for( int k=0; k < m_nz; ++k )
			{
				int idx = coordsToIndex( i, j, k );
				result.segment<3>( 3 * idx ) = m_gridMasses[ idx ] * vTransformed.segment<3>( 3 * idx ) - m_timeStep * df.segment<3>( 3 * idx );
			}
		}
	}
	
	// divide output:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		if( m_gridMasses[i] != 0 )
		{
			result.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
	}
}

void Grid::calculateForceDifferentials( const ParticleData& d, const VectorXf& dx, VectorXf& df ) const
{
	df.setZero();

	// little array with indexes going from -1 to store shape function derivative weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	DECLARE_WEIGHTARRAY( dw );
	Vector3i particleCell;
	
	for( size_t p = 0; p < d.particleF.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w, dw );
		
		// work out deformation gradient differential for this particle when grid nodes are
		// moved by their respective v * Dt
		Matrix3f dFp = Matrix3f::Zero();
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					Vector3f weightGrad( dw[0][i] * w[1][j] * w[2][k], w[0][i] * dw[1][j] * w[2][k], w[0][i] * w[1][j] * dw[2][k] );
					Vector3f deltaX = dx.segment<3>( 3 * idx );
					dFp += deltaX * weightGrad.transpose() * d.particleF[p];
				}
			}
		}

		Matrix3f forceMatrix;
		m_constituativeModel.forceDifferentialDensity( forceMatrix, dFp, d, p );
		forceMatrix = d.particleVolumes[p] * forceMatrix * d.particleF[p].transpose();

		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					Vector3f weightGrad( dw[0][i] * w[1][j] * w[2][k], w[0][i] * dw[1][j] * w[2][k], w[0][i] * w[1][j] * dw[2][k] );
					
					// add on difference in velocity due to this force:
					df.segment<3>( 3 * idx ) -= forceMatrix * weightGrad;
				}
			}
		}
	}
}

void Grid::calculateForces( const ParticleData& d, VectorXf& forces ) const
{
	// little array with indexes going from -1 to store shape function derivative weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	DECLARE_WEIGHTARRAY( dw );
	Vector3i particleCell;
	
	// start with gravity:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		forces.segment<3>(3 * i) = m_gridMasses[i] * Vector3f( 0, GRAVITY, 0 );
	}
	
	// add on internal forces:
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w, dw );
		
		Matrix3f dEdF;
		m_constituativeModel.dEnergyDensitydF( dEdF, d, p );
		
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					Vector3f weightGrad( dw[0][i] * w[1][j] * w[2][k], w[0][i] * dw[1][j] * w[2][k], w[0][i] * w[1][j] * dw[2][k] );
					forces.segment<3>( 3 * idx ) -= d.particleVolumes[p] * dEdF * d.particleF[p].transpose() * weightGrad;
				}
			}
		}
	}
}

float Grid::calculateEnergy( const ParticleData& d ) const
{
	float e = 0;
	for( size_t p=0; p < d.particleF.size(); ++p )
	{
		e += d.particleVolumes[p] * m_constituativeModel.energyDensity( d, p );
	}
	return e;
}

//#define DIAGONALPRECONDITIONED 1

unsigned Grid::matrixTexture( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects ) const
{
	VectorXf x( m_gridVelocities.size() );
	x.setZero();
	
	VectorXf b( m_gridVelocities.size() );
	MatrixXf M( m_gridVelocities.size(), m_gridVelocities.size() );
	
	for( int i=0; i < x.size(); ++i )
	{
		std::cerr << i << " of " << x.size() << " mass = " << m_gridMasses[i/3] << std::endl;

		x[i] = 1;
		applyImplicitUpdateMatrix(d, collisionObjects, x, b );
		x[i] = 0;

#ifdef DIAGONALPRECONDITIONED

		for( int j=0; j < x.size() / 3; ++j )
		{
			if( m_gridMasses[j] != 0 )
			{
				b.segment<3>( 3 * j ) /= sqrt( m_gridMasses[j] );
			}
		}
		if( m_gridMasses[i / 3] != 0 )
		{
			M.block( 0, i, m_gridVelocities.size(), 1 ) = b / sqrt( m_gridMasses[i / 3] );
		}
		else
		{
			M.block( 0, i, m_gridVelocities.size(), 1 ) = b;
		}
#else
		M.block( 0, i, m_gridVelocities.size(), 1 ) = b;
#endif
	}

	MatrixXf shouldBeZero = M.transpose() - M;
	std::cerr << shouldBeZero.maxCoeff() << " - " << shouldBeZero.minCoeff() << std::endl;
	
	float maxM = M.maxCoeff();
	float minM = M.minCoeff();
	float norm = fabs( minM );
	if( fabs( maxM ) > norm )
	{
		norm = fabs( maxM );
	}

	std::vector<unsigned char> sparsity;
	int texW = int( ceil( (float)x.size() / 2 ) * 2 );
	int texH = texW;

	for(int j=0; j < texH; ++j)
	{
		for(int i=0; i < texW; ++i)
		{
			if( i >= x.size() || j >= x.size() )
			{
				sparsity.push_back( 0 );
				continue;
			}
			
			sparsity.push_back( (unsigned char)( fabs( M(i,j) ) * 255 / norm ) );
		}
	}

	GLuint matrixTexture;
	
	glEnable( GL_TEXTURE_2D );
	glGenTextures( 1, &matrixTexture );
	glBindTexture(GL_TEXTURE_2D, matrixTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, texW, texH, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, &sparsity[0] );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable( GL_TEXTURE_2D );
	
	std::cerr << "compute svd!!" << std::endl;
	
	JacobiSVD<MatrixXf> svd(M, Eigen::ComputeFullV);
	std::cerr << "done!" << std::endl;
	MatrixXf V = svd.matrixV();
	
	int numNulls = 0;
	for( int i=0; i < svd.singularValues().size(); ++i )
	{
		std::cerr << svd.singularValues()[i] << std::endl;
		if( svd.singularValues()[i] == 0 )
		{
			++numNulls;
		}
	}



	return matrixTexture;
}

void Grid::testForces( const ParticleData& d )
{
	// save the state so we don't screw the sim up:
	VectorXf originalGridVelocities = m_gridVelocities;

	// calculate da forces brah!
	VectorXf forces( m_gridVelocities.size() );
	calculateForces( d, forces );
	
	// calculate unperturbed energy:
	float e0 = calculateEnergy( d );
	
	// now we're gonna calculate energy derivatives... the stupid way!
	// we're gonna do this component by component, and we're gonna do it
	// by zeroing out the grid velocities, setting the component we're gonna
	// test to delta/m_timeStep, advancing bits of the sim with that velocity field,
	// calculating the energy in the final state (in which one of the grid nodes
	// will have moved a distance delta along one of the axes), and using the result
	// to calculate a finite difference derivative!
	float delta = 0.01f;
	for( int idx = 0; idx < m_gridMasses.size(); ++idx )
	{
		for( size_t dim = 0; dim < 3; ++dim )
		{
			ParticleData dTest = d;
			m_gridVelocities.setZero();
			
			// perturb current grid point a distance delta along the current axis,
			// and calculate the resulting deformation gradients:
			m_gridVelocities( 3 * idx + dim ) = delta / m_timeStep;
			updateDeformationGradients( dTest );
			
			// calculate the resulting energy:
			float e = calculateEnergy( dTest );
			
			// so force = -dE/dX = ( e0 - e ) / delta
			float f = ( e0 - e ) / delta;
			std::cerr << f << " == " << forces( 3 * idx + dim ) << "?  " << (3 * idx + dim) << " of " << forces.size() << std::endl;
		}
	}

	m_gridVelocities = originalGridVelocities;

}

void Grid::testForceDifferentials( const ParticleData& d )
{
	// calculate da forces brah!
	VectorXf forces( m_gridVelocities.size() );
	calculateForces( d, forces );

	// small random perturbation on the grid nodes:
	VectorXf dx( m_gridVelocities.size() );
	dx.setRandom();
	dx = dx * 0.01f;
	
	// calculate force differentials resulting from this perturbation:
	VectorXf forceDifferentials( m_gridVelocities.size() );
	calculateForceDifferentials( d, dx, forceDifferentials );
	
	// save the state so we don't screw the sim up:
	VectorXf originalGridVelocities = m_gridVelocities;
	ParticleData dTest = d;
	
	m_gridVelocities = dx / m_timeStep;
	updateDeformationGradients( dTest );
	VectorXf perturbedForces( m_gridVelocities.size() );
	calculateForces( dTest, perturbedForces );
	
	VectorXf actualForceDifferentials = perturbedForces - forces;
	
	for( int i=0; i <forceDifferentials.size(); ++i )
	{
		#ifdef WIN32
		Sleep(100);
		#endif
		std::cerr << forceDifferentials[i] << " == " << actualForceDifferentials[i] << "? " << i << " of " << forceDifferentials.size() << std::endl;
	}
	
	m_gridVelocities = originalGridVelocities;
}

void Grid::updateGridVelocities( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects, const Solver& implicitSolver )
{
	m_prevGridVelocities = m_gridVelocities;

	// work out forces on grid points:
	VectorXf forces( m_gridVelocities.size() );
	calculateForces( d, forces );
	
	// work out forward velocity update - that's equation 10 in the paper:
	VectorXf forwardMomenta( m_gridVelocities.size() );
	m_nodeCollided.resize( m_gridMasses.size() );
	for( int i=0; i < m_nx; ++i )
	{
		for( int j=0; j < m_ny; ++j )
		{
			for( int k=0; k < m_nz; ++k )
			{
				int idx = coordsToIndex( i, j, k );
				m_nodeCollided[idx] = false;

				if( m_gridMasses[idx] == 0 )
				{
					forwardMomenta.segment<3>( 3 * idx ).setZero();
					m_gridVelocities.segment<3>( 3 * idx ).setZero();
				}
				else
				{
					Vector3f force = forces.segment<3>( 3 * idx );
					Vector3f velocity = m_gridVelocities.segment<3>( 3 * idx );
					Vector3f forwardVelocity = velocity + m_timeStep * force / m_gridMasses[idx];

					// apply collisions:
					Vector3f x( m_gridH * i + m_xmin, m_gridH * j + m_ymin, m_gridH * k + m_zmin );
					for( size_t objIdx = 0; objIdx < collisionObjects.size(); ++objIdx )
					{
						float phi = collisionObjects[objIdx]->phi( x );
						if( phi <= 0 )
						{
							// intersecting the object
							Vector3f n;
							collisionObjects[objIdx]->grad( x, n );
							n.normalize();
							float nDotV = n.dot( forwardVelocity );
							if( nDotV < 0 )
							{
								// trying to move into the object:
								m_nodeCollided[idx] = true;

								// velocity perpendicular to the object
								Vector3f vPerp = nDotV * n;

								// remaining component is velocity paralell to the object:
								Vector3f vTangent = forwardVelocity - vPerp;

								forwardVelocity = vTangent * ( 1 + COULOMBFRICTION * nDotV / vTangent.norm() );
							}
						}
					}
					
					m_gridVelocities.segment<3>( 3 * idx ) = forwardVelocity;
					forwardMomenta.segment<3>( 3 * idx ) = forwardVelocity * m_gridMasses[idx];
				}
			}
		}
	}
	
	// So we want to solve 
	// m * v^(n+1) - m_timeStep * dF(v^(n+1) * m_timeStep) = forwardMomenta
	
	// these crazy "dividing by the square roots of the masses" shenanigans
	// don't really correspond to anything physical - they're just meant to
	// give the matrix better eigenvalues so the solver's happy with it...
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		m_gridVelocities.segment<3>( 3 * i ) *= sqrt( m_gridMasses[i] );
		if( m_gridMasses[i] != 0 )
		{
			forwardMomenta.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
	}

	implicitSolver(
		this,
		d,
		collisionObjects,
		forwardMomenta,
		m_gridVelocities );
	
	// more funky remapping to make the solver happy:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		if( m_gridMasses[i] != 0 )
		{
			m_gridVelocities.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
	}
}


void Grid::updateDeformationGradients( ParticleData& d )
{
	// little array with indexes going from -1 to store shape function derivative weights
	// on each dimension:
	DECLARE_WEIGHTARRAY( w );
	DECLARE_WEIGHTARRAY( dw );
	Vector3i particleCell;
	
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		cellAndWeights( d.particleX[p], particleCell, w, dw );
		
		Vector3f v = Vector3f::Zero();
		Matrix3f delV = Matrix3f::Zero();
		for( int i=-1; i < 3; ++i )
		{
			for( int j=-1; j < 3; ++j )
			{
				for( int k=-1; k < 3; ++k )
				{
					int idx = coordsToIndex( particleCell[0] + i, particleCell[1] + j, particleCell[2] + k );
					Vector3f weightGrad( dw[0][i] * w[1][j] * w[2][k], w[0][i] * dw[1][j] * w[2][k], w[0][i] * w[1][j] * dw[2][k] );
					Vector3f vSample = m_gridVelocities.segment<3>( 3 * idx );
					delV += vSample * weightGrad.transpose();
					v += vSample * w[0][i] * w[1][j] * w[2][k];
				}
			}
		}
		Matrix3f newParticleF = ( Matrix3f::Identity() + m_timeStep * delV ) * d.particleF[p];
		d.particleF[p] = newParticleF;
	}

	m_constituativeModel.updateDeformation( d );
}


inline float Grid::N( float x )
{
	float ax = fabs(x);
	if( ax < 1 )
	{
		return 0.5f * ax * ax * ax - ax * ax + 2.0f/3;
	}
	else if( ax < 2 )
	{
		return -ax * ax * ax / 6 + ax * ax - 2 * ax + 4.0f/3;
	}
	else
	{
		return 0;
	}
}

inline float Grid::DN( float x )
{
	if( x < 0 )
	{
		return -DN( -x );
	}
	
	if( x < 1 )
	{
		return x * ( 1.5f * x - 2 );
	}
	else if( x < 2 )
	{
		x -= 2;
		return -0.5f * x * x;
	}
	else
	{
		return 0;
	}
}

inline int Grid::coordsToIndex( int x, int y, int z ) const
{
	return x + m_nx * ( y + z * m_ny );
}

void Grid::cellAndWeights( const Vector3f& particleX, Vector3i& particleCell, float *w[], float** dw ) const
{
	Vector3f positionInCell;
	positionInCell[0] = ( particleX[0] - m_xmin ) / m_gridH;
	positionInCell[1] = ( particleX[1] - m_ymin ) / m_gridH;
	positionInCell[2] = ( particleX[2] - m_zmin ) / m_gridH;
	
	particleCell[0] = (int)floor( positionInCell[0] );
	particleCell[1] = (int)floor( positionInCell[1] );
	particleCell[2] = (int)floor( positionInCell[2] );
	
	positionInCell -= particleCell.cast<float>();
	if( dw )
	{
		for( int i=0; i < 3; ++i )
		{
			dw[i][-1] = DN( positionInCell[i] + 1 ) / m_gridH;
			dw[i][0] = DN( positionInCell[i] ) / m_gridH;
			dw[i][1] = DN( positionInCell[i] - 1 ) / m_gridH;
			dw[i][2] = DN( positionInCell[i] - 2 ) / m_gridH;
		}
	}
	
	for( int i=0; i < 3; ++i )
	{
		w[i][-1] = N( positionInCell[i] + 1 );
		w[i][0] = N( positionInCell[i] );
		w[i][1] = N( positionInCell[i] - 1 );
		w[i][2] = N( positionInCell[i] - 2 );
	}
	
}

inline void Grid::minMax( float x, float& min, float& max )
{
	if( x < min )
	{
		min = x;
	}
	if( x > max )
	{
		max = x;
	}
}

inline int Grid::fixDim( float& min, float& max ) const
{
	float minPadded = min - 1.5f * m_gridH;
	float maxPadded = max + 1.5f * m_gridH;
	int n = int( ceil( ( maxPadded - minPadded ) / m_gridH ) ) + 1;
	min = minPadded;
	max = min + n * m_gridH;
	return n;
}
