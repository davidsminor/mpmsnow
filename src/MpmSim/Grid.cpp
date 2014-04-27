#include "MpmSim/Grid.h"
#include "MpmSim/ImplicitUpdateMatrix.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>

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
		Grid::PointToGridIterator& shIt = g.pointIterator();
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
		Grid::PointToGridIterator& shIt = g.pointIterator();
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



Grid::Grid(
		ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		float timeStep,
		const ShapeFunction& shapeFunction,
		const ConstitutiveModel& model,
		int dimension
) :
	m_d( d ),
	m_collisionObjects( collisionObjects ),
	m_gridH( d.gridSize ),
	m_timeStep( timeStep ),
	m_dimension( dimension ),
	m_shapeFunction( shapeFunction ),
	m_constitutiveModel( model )
{
	// work out the physical size of the grid:
	m_min.setConstant( 1.e10 );
	m_max.setConstant( -1.e10 );

	for( size_t i=0; i < d.particleX.size(); ++i )
	{
		float prod = d.particleX[i][0] * d.particleX[i][1] * d.particleX[i][2];
#ifdef WIN32
		if( !_finite(prod) )
#else

		if( isinff(prod) || isnanf(prod) )
#endif
		{
			throw std::runtime_error( "grid has non finite dimensions!" );
		}
		for( int j=0; j < m_dimension; ++j )
		{
			minMax( d.particleX[i][j], m_min[j], m_max[j] );
		}
	}
	
	Vector3f averageV = Vector3f::Zero();
	for( size_t i=0; i < d.particleV.size(); ++i )
	{
		averageV += d.particleV[i];
	}
	averageV[0] /= (int)d.particleV.size();
	averageV[1] /= (int)d.particleV.size();
	averageV[2] /= (int)d.particleV.size();

	std::cerr << "average v: " << averageV.transpose() << std::endl;
	std::cerr << "grid bbox: " << m_min.transpose() << " --> " << m_max.transpose() << std::endl;

	// calculate grid dimensions and quantize bounding box:
	m_n.setConstant(1);
	for( int j=0; j < m_dimension; ++j )
	{
		m_n[j] = fixDim( m_min[j], m_max[j] );
	}
	
	long long ncells = (long long)m_n[0] * (long long)m_n[1] * (long long)m_n[2];
	if( ncells > 4000000000 || ncells <= 0 )
	{
		throw std::runtime_error( "grid is too big" );
	}
	
	// calculate masses:
	m_gridMasses.resize( ncells );
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
	m_gridVelocities.resize( ncells * 3 );
	m_gridVelocities.setZero();
	
	splat< VelocitySplatter >( d, m_gridVelocities );
	
}

const Eigen::VectorXf& Grid::masses() const
{
	return m_gridMasses;
}

const Eigen::VectorXf& Grid::getVelocities() const
{
	return m_gridVelocities;
}
void Grid::setVelocities( const Eigen::VectorXf& velocities )
{
	m_gridVelocities = velocities;
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
	for( int i=0; i <= m_n[0]; ++i )
	{
		for( int j=0; j <= m_n[1]; ++j )
		{
			glVertex3f( m_min[0] + i * m_gridH, m_min[1] + j * m_gridH, m_min[2] );
			glVertex3f( m_min[0] + i * m_gridH, m_min[1] + j * m_gridH, m_max[2] );
		}
	}
	// zy
	for( int i=0; i <= m_n[2]; ++i )
	{
		for( int j=0; j <= m_n[1]; ++j )
		{
			glVertex3f( m_min[0], m_min[1] + j * m_gridH, m_min[2] + i * m_gridH );
			glVertex3f( m_max[0], m_min[1] + j * m_gridH, m_min[2] + i * m_gridH );
		}
	}

	// xz
	for( int i=0; i <= m_n[0]; ++i )
	{
		for( int j=0; j <= m_n[2]; ++j )
		{
			glVertex3f( m_min[0] + i * m_gridH, m_min[1], m_min[2] + j * m_gridH );
			glVertex3f( m_min[0] + i * m_gridH, m_max[1], m_min[2] + j * m_gridH );
		}
	}
	glEnd();
	
	glEnable( GL_LIGHTING );
	for( size_t i=0; i < m_collisionObjects.size(); ++i )
	{
		m_collisionObjects[i]->draw();
	}
	glDisable( GL_LIGHTING );
	
	for( size_t p = 0; p < m_d.particleX.size(); ++p )
	{
		float r = 2 * pow( m_d.particleVolumes[p], 1.0f/3 ) / ( 4 * 3.1415926 / 3 );
		Eigen::Vector3f x = m_d.particleF[p] * Eigen::Vector3f(1,0,0);
		Eigen::Vector3f y = m_d.particleF[p] * Eigen::Vector3f(0,1,0);
		Eigen::Vector3f z = m_d.particleF[p] * Eigen::Vector3f(0,0,1);
		
		glBegin( GL_QUADS );
		glColor3f( 1,0,1 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
		glEnd();
		
		glBegin( GL_LINE_LOOP );
		glColor3f( 1,1,1 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
		glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
		glEnd();

	}
}

void Grid::computeDensities() const
{
	m_d.particleDensities.resize( m_d.particleX.size(), 0 );
	
	Grid::PointToGridIterator& shIt = pointIterator();
	Vector3i particleCell;
	
	float cellVolume = m_gridH * m_gridH * m_gridH;
	for( size_t p = 0; p < m_d.particleX.size(); ++p )
	{
		shIt.initialize( m_d.particleX[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell );
			
			// accumulate the particle's density:
			m_d.particleDensities[p] += shIt.w() * m_gridMasses[ idx ] / cellVolume;
					
		} while( shIt.next() );
	}
}

void Grid::updateParticleVelocities()
{
	Grid::PointToGridIterator& shIt = pointIterator();
	Vector3i particleCell;
	
	const float alpha = 0.95f;
	float maxv = 0;
	for( size_t p = 0; p < m_gridVelocities.size() / 3; ++p )
	{
		float v = m_gridVelocities.segment<3>(p*3).norm();
		if( v > maxv )
		{
			v = maxv;
		}
	}
	std::cerr << "updateParticleVelocities: " << maxv << std::endl;
	
	// blend FLIP and PIC, as pure FLIP allows spurious particle motion
	// inside the cells:
	for( size_t p = 0; p < m_d.particleX.size(); ++p )
	{
		Vector3f vFlip = m_d.particleV[p];
		Vector3f vPic = Vector3f::Zero();
		shIt.initialize( m_d.particleX[p] );
		do
		{
			// soo... should I be interpolating momentum here instead? Need to experiment...
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell );
			float w = shIt.w();
			vFlip += w * ( m_gridVelocities.segment<3>( 3 * idx ) - m_prevGridVelocities.segment<3>( 3 * idx ) );
			vPic += w * m_gridVelocities.segment<3>( 3 * idx );
		} while( shIt.next() );
		m_d.particleV[p] = alpha * vFlip + ( 1.0f - alpha ) * vPic;
		
		collide( m_d.particleV[p], m_d.particleX[p], m_collisionObjects );
	}
}

float Grid::gridH() const
{
	return m_gridH;
}

void Grid::origin( Eigen::Vector3f& o ) const
{
	o = m_min;
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
		Grid::PointToGridIterator& shIt = g.pointIterator();
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
		Grid::PointToGridIterator& shIt = g.pointIterator();
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


bool Grid::collide( Eigen::Vector3f& v, const Eigen::Vector3f& x, const std::vector<CollisionObject*>& collisionObjects )
{
	bool nodeCollided = false;
	for( size_t objIdx = 0; objIdx < collisionObjects.size(); ++objIdx )
	{
		float phi = collisionObjects[objIdx]->phi( x );
		if( phi <= 0 )
		{
			// intersecting the object
			Vector3f vObj;
			collisionObjects[objIdx]->velocity( x, vObj );
			
			Vector3f n;
			collisionObjects[objIdx]->grad( x, n );
			n.normalize();
			
			// subtract off object velocity:
			v -= vObj;
			
			float nDotV = n.dot( v );
			if( nDotV < 0 )
			{
				// trying to move into the object:
				nodeCollided = true;

				// velocity perpendicular to the object
				Vector3f vPerp = nDotV * n;

				// remaining component is velocity paralell to the object:
				Vector3f vTangent = v - vPerp;
				float vtNorm = vTangent.norm();
				if( vtNorm >= -nDotV * COULOMBFRICTION )
				{
					v = vTangent * ( 1 + COULOMBFRICTION * nDotV / vTangent.norm() );
				}
				else
				{
					v.setZero();
				}
			}
			
			// add object velocity back on:
			v += vObj;
		}
	}
	
	return nodeCollided;
}


void Grid::updateGridVelocities( const LinearSolver& implicitSolver )
{
	m_prevGridVelocities = m_gridVelocities;

	// work out forces on grid points:
	VectorXf forces( m_gridVelocities.size() );
	calculateForces( m_d, forces );
	
	// work out forward velocity update - that's equation 10 in the paper:
	VectorXf forwardMomenta( m_gridVelocities.size() );
	m_nodeCollided.resize( m_gridMasses.size() );
	for( int i=0; i < m_n[0]; ++i )
	{
		for( int j=0; j < m_n[1]; ++j )
		{
			for( int k=0; k < m_n[2]; ++k )
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
					Vector3f x( m_gridH * i + m_min[0], m_gridH * j + m_min[1], m_gridH * k + m_min[2] );
					m_nodeCollided[ idx ] = collide( forwardVelocity, x, m_collisionObjects );
					
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
		if( m_gridMasses[i] > 1.e-8 )
		{
			forwardMomenta.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
	}
	
	// this is broken for moving collision objects! Currently the implicit update matrix thinks they're not moving...
	implicitSolver(
		ImplicitUpdateMatrix( m_d, *this, m_collisionObjects ),
		forwardMomenta,
		m_gridVelocities );
	
	// more funky remapping to make the solver happy:
	for( int i=0; i < m_gridMasses.size(); ++i )
	{
		if( m_gridMasses[i] > 1.e-8 )
		{
			m_gridVelocities.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
		else
		{
			m_gridVelocities.segment<3>( 3 * i ).setZero();
		}
	}
}


float Grid::updateDeformationGradients()
{

	Grid::PointToGridIterator& shIt = pointIterator();
	Vector3f weightGrad;
	Vector3i particleCell;
	Matrix3f delV;
	float maxMovement(0);

	for( size_t p = 0; p < m_d.particleX.size(); ++p )
	{
		delV.setZero();
		shIt.initialize( m_d.particleX[p], true );
		do
		{
			shIt.gridPos( particleCell );
			shIt.dw( weightGrad );
			int idx = coordsToIndex( particleCell );
			delV += m_gridVelocities.segment<3>( 3 * idx ) * weightGrad.transpose();
		} while( shIt.next() );
		
		float movement = 0.5f * ( delV + delV.transpose() ).norm();
		if( movement > maxMovement )
		{
			maxMovement = movement;
		}

		Matrix3f newParticleF = ( Matrix3f::Identity() + m_timeStep * delV ) * m_d.particleF[p];
		m_d.particleF[p] = newParticleF;
	}

	m_constitutiveModel.updateDeformation( m_d );
	return maxMovement;
}

int Grid::coordsToIndex( const Eigen::Vector3i& p ) const
{
	return p[0] + m_n[0] * ( p[1] + m_n[1] * p[2] );
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
	int cellMin = int( floor( min / m_gridH ) ) - m_shapeFunction.supportRadius() + 1;
	int cellMax = int( ceil( max / m_gridH ) ) + m_shapeFunction.supportRadius() + 1;
	min = cellMin * m_gridH;
	max = cellMax * m_gridH;
	return cellMax - cellMin;
}

Grid::PointToGridIterator& Grid::pointIterator() const
{
	std::auto_ptr< PointToGridIterator >& pIt = m_pointIterators.local();
	if( pIt.get() == 0 )
	{
		pIt.reset( new PointToGridIterator( *this ) );
	}
	return *pIt.get();
}

Grid::PointToGridIterator::PointToGridIterator( const Grid& g ) :
	m_grid( g ),
	m_diameter( 2 * g.m_shapeFunction.supportRadius() )
{
	for( int dim=0; dim < 3; ++dim )
	{
		m_w[dim].resize( m_diameter, 1 );
		m_dw[dim].resize( m_diameter, 0 );
	}
}

void Grid::PointToGridIterator::initialize( const Vector3f& p, bool computeDerivatives )
{
	m_gradients = computeDerivatives;
	int r = (int)m_w[0].size() / 2;
	m_pos.setZero();
	m_base.setZero();

	for( int dim=0; dim < m_grid.m_dimension; ++dim )
	{
		float fracDimPos = ( p[dim] - m_grid.m_min[dim] ) / m_grid.m_gridH;
		int dimPos = (int)floor( fracDimPos );
		fracDimPos -= dimPos;
		int j;
		
		m_base[ dim ] = dimPos + 1 - r;
		
		j = 1-r;
		for( int i = 0; i < m_diameter; ++i, ++j )
		{
			m_w[dim][i] = m_grid.m_shapeFunction.w( j - fracDimPos );
		}

		if( computeDerivatives )
		{
			j = 1-r;
			for( int i = 0; i < m_diameter; ++i, ++j )
			{
				m_dw[dim][i] = m_grid.m_shapeFunction.dw( j - fracDimPos ) / m_grid.m_gridH;
			}
		}
	}
}

bool Grid::PointToGridIterator::next()
{

	++m_pos[0];
	if( m_pos[0] >= m_diameter )
	{
		if( m_grid.m_dimension == 1 )
		{
			return false;
		}

		m_pos[0] = 0;
		++m_pos[1];
		
		if( m_pos[1] >= m_diameter )
		{
			
			if( m_grid.m_dimension == 2 )
			{
				return false;
			}

			m_pos[1] = 0;
			++m_pos[2];
		}
		if( m_pos[2] >= m_diameter )
		{
			return false;
		}
	}

	return true;
}

void Grid::PointToGridIterator::gridPos( Eigen::Vector3i& pos ) const
{
	pos[0] = m_pos[0] + m_base[0];
	pos[1] = m_pos[1] + m_base[1];
	pos[2] = m_pos[2] + m_base[2];
}

void Grid::PointToGridIterator::dw( Eigen::Vector3f& g ) const
{
	if( !m_gradients )
	{
		throw std::runtime_error( "Grid::PointToGridIterator::dw(): derivatives not computed!" );
	}
	g[0] = -m_dw[0][ m_pos[0] ] *  m_w[1][ m_pos[1] ] *  m_w[2][ m_pos[2] ];
	g[1] = - m_w[0][ m_pos[0] ] * m_dw[1][ m_pos[1] ] *  m_w[2][ m_pos[2] ];
	g[2] = - m_w[0][ m_pos[0] ] *  m_w[1][ m_pos[1] ] * m_dw[2][ m_pos[2] ];
}

float Grid::PointToGridIterator::w() const
{
	return m_w[0][ m_pos[0] ] * m_w[1][ m_pos[1] ] * m_w[2][ m_pos[2] ];
}

