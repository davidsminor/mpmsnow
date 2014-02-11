#include "MpmSim/Grid.h"
#include "MpmSim/ImplicitUpdateMatrix.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>
#include <fstream>

#include <GL/gl.h>

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>


#ifdef HAVE_CORTEX
#include "IECore/ImagePrimitive.h"
#include "IECore/Writer.h"
#endif

using namespace Eigen;
using namespace MpmSim;

#define GRAVITY -9.8f
#define COULOMBFRICTION 0.5f

Grid::GridSplatter::GridSplatter(
	const Grid& g,
	const ParticleData& d,
	const ParticleData::PartitionList& partition,
	Eigen::VectorXf& result,
	const void* args
) :
	m_g( g ),
	m_d( d ),
	m_partition( partition ),
	m_result( result ),
	m_args( args )
{
}

void Grid::GridSplatter::operator()(const tbb::blocked_range<int> &r) const
{
	// iterate over a range of voxels assigned by tbb:
	for (int i = r.begin(); i != r.end(); ++i)
	{
		// splat all the particles in this voxel:
		splat( m_partition[i].first, m_partition[i].second, m_g, m_d, m_result, m_args );
	}
}

class MassSplatter : public Grid::GridSplatter
{
public:
	DECLARE_GRIDSPLATTER_CONSTRUCTOR( MassSplatter )
	
	virtual void splat(
		ParticleData::IndexIterator begin,
		ParticleData::IndexIterator end,
		const Grid& g,
		const ParticleData& d,
		Eigen::VectorXf& gridMasses,
		const void* args
	) const
	{
		ShapeFunction::PointToGridIterator& shIt = g.pointIterator();
		Vector3i particleCell;
		for( ParticleData::IndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			shIt.initialize( d.particleX[p] );
			do
			{
				shIt.gridPos( particleCell );
				int idx = g.coordsToIndex( particleCell );
				gridMasses[ idx ] += d.particleM[p] * shIt.w();
			} while( shIt.next() );
		}
	}
			
};


class VelocitySplatter : public Grid::GridSplatter
{
public:
	
	DECLARE_GRIDSPLATTER_CONSTRUCTOR( VelocitySplatter )
	
	virtual void splat(
		ParticleData::IndexIterator begin,
		ParticleData::IndexIterator end,
		const Grid& g,
		const ParticleData& d,
		Eigen::VectorXf& gridVelocities,
		const void* args
	) const
	{
		ShapeFunction::PointToGridIterator& shIt = g.pointIterator();
		Vector3i particleCell;
		const Eigen::VectorXf& masses = g.masses();
		for( ParticleData::IndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			shIt.initialize( d.particleX[p] );
			
			// splat the velocities:
			do
			{
				shIt.gridPos( particleCell );
				int idx = g.coordsToIndex( particleCell );
				float gridCellMass = masses[idx];
				if( gridCellMass > 0 )
				{
					gridVelocities.segment<3>( 3 * idx ) += shIt.w() * ( d.particleM[p] / gridCellMass ) * d.particleV[p];
				}
			} while( shIt.next() );
		}
	}
			
};



Grid::Grid( const ParticleData& d, float timeStep, const ShapeFunction& shapeFunction, const ConstitutiveModel& model )
	: m_gridH( d.gridSize ), m_timeStep( timeStep ), m_shapeFunction( shapeFunction ), m_constitutiveModel( model )
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
	
	m_gridOrigin[0] = m_xmin;
	m_gridOrigin[1] = m_ymin;
	m_gridOrigin[2] = m_zmin;
	
	// calculate masses:
	m_gridMasses.resize( m_nx * m_ny * m_nz );
	m_gridMasses.setZero();
	
	splat< MassSplatter >( d, m_gridMasses );

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
	
	splat< VelocitySplatter >( d, m_gridVelocities );
	
}

const Eigen::VectorXf& Grid::masses() const
{
	return m_gridMasses;
}

const Eigen::VectorXf& Grid::velocities() const
{
	return m_gridVelocities;
}

const ConstitutiveModel& Grid::constitutiveModel() const
{
	return m_constitutiveModel;
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
	
	ShapeFunction::PointToGridIterator& shIt = pointIterator();
	Vector3i particleCell;
	
	float cellVolume = m_gridH * m_gridH * m_gridH;
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		shIt.initialize( d.particleX[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell );
			
			// accumulate the particle's density:
			d.particleDensities[p] += shIt.w() * m_gridMasses[ idx ] / cellVolume;
					
		} while( shIt.next() );
	}
}

void Grid::updateParticleVelocities( ParticleData& d )
{
	ShapeFunction::PointToGridIterator& shIt = pointIterator();
	Vector3i particleCell;

	// this is, like, totally doing things FLIP style. The paper recommends a combination of FLIP and PIC...
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		shIt.initialize( d.particleX[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell );
			d.particleV[p] += shIt.w() * ( m_gridVelocities.segment<3>( 3 * idx ) - m_prevGridVelocities.segment<3>( 3 * idx ) );
		} while( shIt.next() );
	}
}

float Grid::gridH() const
{
	return m_gridH;
}

void Grid::origin( Eigen::Vector3f& o ) const
{
	o[0] = m_xmin;
	o[1] = m_ymin;
	o[2] = m_zmin;
}

class ForceDifferentialSplatter : public Grid::GridSplatter
{
public:
	
	DECLARE_GRIDSPLATTER_CONSTRUCTOR( ForceDifferentialSplatter )
	
	virtual void splat(
		ParticleData::IndexIterator begin,
		ParticleData::IndexIterator end,
		const Grid& g,
		const ParticleData& d,
		Eigen::VectorXf& df,
		const void* args
	) const
	{
		const Eigen::VectorXf& dx = *(Eigen::VectorXf*)args;
		ShapeFunction::PointToGridIterator& shIt = g.pointIterator();
		Vector3i particleCell;
		Vector3f weightGrad;
		for( ParticleData::IndexIterator it = begin; it != end; ++it )
		{
			int p = *it;

			// work out deformation gradient differential for this particle when grid nodes are
			// moved by their respective v * Dt
			Matrix3f dFp = Matrix3f::Zero();
			shIt.initialize( d.particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = g.coordsToIndex( particleCell );
				dFp += dx.segment<3>( 3 * idx ) * weightGrad.transpose() * d.particleF[p];
			} while( shIt.next() );
			
			Matrix3f forceMatrix;
			g.constitutiveModel().dEdFDifferential( forceMatrix, dFp, d, p );
			forceMatrix = d.particleVolumes[p] * forceMatrix * d.particleF[p].transpose();

			shIt.initialize( d.particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = g.coordsToIndex( particleCell );
				
				// add on difference in velocity due to this force:
				df.segment<3>( 3 * idx ) -= forceMatrix * weightGrad;
						
			} while( shIt.next() );
		}
	}
			
};


void Grid::calculateForceDifferentials( const ParticleData& d, const VectorXf& dx, VectorXf& df ) const
{
	df.setZero();
	splat< ForceDifferentialSplatter >( d, df, &dx );
}


class ForceSplatter : public Grid::GridSplatter
{
public:
	
	DECLARE_GRIDSPLATTER_CONSTRUCTOR( ForceSplatter )
	
	virtual void splat(
		ParticleData::IndexIterator begin,
		ParticleData::IndexIterator end,
		const Grid& g,
		const ParticleData& d,
		Eigen::VectorXf& forces,
		const void* args
	) const
	{
		ShapeFunction::PointToGridIterator& shIt = g.pointIterator();
		Vector3i particleCell;
		Vector3f weightGrad;

		for( ParticleData::IndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			Matrix3f dEdF;
			g.constitutiveModel().dEnergyDensitydF( dEdF, d, p );
			
			shIt.initialize( d.particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = g.coordsToIndex( particleCell );
				forces.segment<3>( 3 * idx ) -= d.particleVolumes[p] * dEdF * d.particleF[p].transpose() * weightGrad;
			} while( shIt.next() );
		}
	}
};



void Grid::calculateForces( const ParticleData& d, VectorXf& forces ) const
{

	// start with gravity:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		forces.segment<3>(3 * i) = m_gridMasses[i] * Vector3f( 0, GRAVITY, 0 );
	}
	
	// add on internal forces:
	splat< ForceSplatter >( d, forces );
}

float Grid::calculateEnergy( const ParticleData& d ) const
{
	float e = 0;
	for( size_t p=0; p < d.particleF.size(); ++p )
	{
		e += d.particleVolumes[p] * m_constitutiveModel.energyDensity( d, p );
	}
	return e;
}


void Grid::outputDiagnostics( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects ) const
{
	std::ofstream f( "/tmp/diagnostics.dat", std::ofstream::out );
	
	
	f << "grid masses: " << std::endl;
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		f << m_gridMasses[i] << std::endl;
	}
	
	f << "particle data: " << std::endl;
	float maxv = 0;
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		f << "p: " << d.particleX[p].transpose() << " v: " << d.particleV[p].transpose() << " m: " << d.particleM[p] << " J: " << d.particleJ[p] << std::endl;
		float v = d.particleV[p].norm();
		if( v > maxv )
		{
			v = maxv;
		}
	}
	
	f << "max v " << maxv << std::endl;
	
	VectorXf x( m_gridVelocities.size() );
	x.setZero();
	
	VectorXf b( m_gridVelocities.size() );
	MatrixXf M( m_gridVelocities.size(), m_gridVelocities.size() );
	
	float maxDiagonal( -1000000000 );
	float minDiagonal( 1000000000 );
	
	std::cerr << "evaluate matrix:" << std::endl;
	
	f << "diagonals:" << std::endl;
	ImplicitUpdateMatrix mat( d, *this, collisionObjects );
	for( int i=0; i < x.size(); ++i )
	{
		std::cerr << i << " of " << x.size() << std::endl;

		x[i] = 1;
		mat.multVector( x, b );
		x[i] = 0;

		M.block( 0, i, x.size(), 1 ) = b;
		f << M(i,i) << std::endl;
		if( M(i,i) > maxDiagonal )
		{
			maxDiagonal = M(i,i);
		}
		if( M(i,i) < minDiagonal )
		{
			minDiagonal = M(i,i);
		}
	}
	
	f << "diagonal check: " << minDiagonal << " - " << maxDiagonal << std::endl;

	MatrixXf shouldBeZero = M.transpose() - M;
	f << "symmetry check: " << shouldBeZero.maxCoeff() << " - " << shouldBeZero.minCoeff() << std::endl;
	 
#ifdef HAVE_CORTEX
	
	Imath::Box2i dataWindow;
	dataWindow.min = Imath::V2i( 1 );
	dataWindow.max = Imath::V2i( x.size() );
	Imath::Box2i displayWindow = dataWindow;
	
	IECore::ImagePrimitivePtr image = new IECore::ImagePrimitive( dataWindow, displayWindow );
	IECore::FloatVectorDataPtr matrixData = new IECore::FloatVectorData;
	matrixData->writable().resize( x.size() * x.size() );
	
	std::vector<float>::iterator pixel = matrixData->writable().begin();
	for( int i=0; i < x.size(); ++i )
	{
		for( int j=0; j < x.size(); ++j )
		{
			*pixel++ = M(i,j);
		}
	}
	image->variables["R"] = IECore::PrimitiveVariable( IECore::PrimitiveVariable::Vertex, matrixData );
	
	IECore::Writer::create( image, "/tmp/matrix.exr" )->write();
	
#endif

	/*
	std::cerr << "compute eigenstuffs!!" << std::endl;

	EigenSolver<MatrixXf> es(M, false);
	
	std::cerr << "done!" << std::endl;
	f << "eigenvalues: " << std::endl;
	for( int i=0; i < es.eigenvalues().size(); ++i )
	{
		f << es.eigenvalues()[i] << std::endl;
	}
	*/
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

void Grid::updateGridVelocities( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects, const LinearSolver& implicitSolver )
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
				int idx = coordsToIndex( Vector3i( i, j, k ) );
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
					
					forwardMomenta.segment<3>( 3 * idx ) = forwardVelocity * m_gridMasses[idx];
				}
			}
		}
	}
	
	m_gridVelocities.setZero();
	
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
		ImplicitUpdateMatrix( d, *this, collisionObjects ),
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

	ShapeFunction::PointToGridIterator& shIt = pointIterator();
	Vector3f weightGrad;
	Vector3i particleCell;
	Matrix3f delV;

	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		delV.setZero();
		shIt.initialize( d.particleX[p], true );
		do
		{
			shIt.gridPos( particleCell );
			shIt.dw( weightGrad );
			int idx = coordsToIndex( particleCell );
			delV += m_gridVelocities.segment<3>( 3 * idx ) * weightGrad.transpose();
		} while( shIt.next() );
		
		Matrix3f newParticleF = ( Matrix3f::Identity() + m_timeStep * delV ) * d.particleF[p];
		d.particleF[p] = newParticleF;
	}

	m_constitutiveModel.updateDeformation( d );
}

int Grid::coordsToIndex( const Eigen::Vector3i& p ) const
{
	return p[0] + m_nx * ( p[1] + m_ny * p[2] );
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
	int cellMin = int( floor( min / m_gridH ) ) - 1;
	int cellMax = int( ceil( max / m_gridH ) ) + 3;
	min = cellMin * m_gridH;
	max = cellMax * m_gridH;
	return cellMax - cellMin;
}

ShapeFunction::PointToGridIterator& Grid::pointIterator() const
{
	std::auto_ptr< ShapeFunction::PointToGridIterator >& pIt = m_pointIterators.local();
	if( pIt.get() == 0 )
	{
		pIt.reset( new ShapeFunction::PointToGridIterator( m_shapeFunction, m_gridH, m_gridOrigin ) );
	}
	return *pIt.get();
}
