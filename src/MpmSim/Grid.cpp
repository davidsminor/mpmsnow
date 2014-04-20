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
		float prod = d.particleX[i][0] * d.particleX[i][1] * d.particleX[i][2];
#ifdef WIN32
		if( !_finite(prod) )
#else
		if( isinff(prod) || isnanf(prod) )
#endif
		{
			throw std::runtime_error( "grid has non finite dimensions!" );
		}
		minMax( d.particleX[i][0], m_xmin, m_xmax );
		minMax( d.particleX[i][1], m_ymin, m_ymax );
		minMax( d.particleX[i][2], m_zmin, m_zmax );
	}
	
	Vector3f averageV = Vector3f::Zero();
	for( int i=0; i < d.particleV.size(); ++i )
	{
		averageV += d.particleV[i];
	}
	averageV /= d.particleV.size();

	std::cerr << "average v: " << averageV.transpose() << std::endl;
	std::cerr << "grid bbox: " << m_xmin << ", " << m_ymin << ", " << m_zmin << " --> " << m_xmax << ", " << m_ymax << ", " << m_zmax << std::endl;

	// calculate grid dimensions and quantize bounding box:
	m_nx = fixDim( m_xmin, m_xmax );
	m_ny = fixDim( m_ymin, m_ymax );
	m_nz = fixDim( m_zmin, m_zmax );
	
	m_gridOrigin[0] = m_xmin;
	m_gridOrigin[1] = m_ymin;
	m_gridOrigin[2] = m_zmin;
	
	long long ncells = (long long)m_nx * (long long)m_ny * (long long)m_nz;
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

void Grid::updateParticleVelocities( ParticleData& d, const std::vector<CollisionObject*>& collisionObjects )
{
	ShapeFunction::PointToGridIterator& shIt = pointIterator();
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
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		Vector3f vFlip = d.particleV[p];
		Vector3f vPic = Vector3f::Zero();
		shIt.initialize( d.particleX[p] );
		do
		{
			// soo... should I be interpolating momentum here instead? Need to experiment...
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell );
			float w = shIt.w();
			vFlip += w * ( m_gridVelocities.segment<3>( 3 * idx ) - m_prevGridVelocities.segment<3>( 3 * idx ) );
			vPic += w * m_gridVelocities.segment<3>( 3 * idx );
		} while( shIt.next() );
		d.particleV[p] = alpha * vFlip + ( 1.0f - alpha ) * vPic;
		
		collide( d.particleV[p], d.particleX[p], collisionObjects );
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
					m_nodeCollided[ idx ] = collide( forwardVelocity, x, collisionObjects );
					
					forwardMomenta.segment<3>( 3 * idx ) = forwardVelocity * m_gridMasses[idx];
					m_gridVelocities.segment<3>( 3 * idx ) = forwardVelocity;
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
		if( m_gridMasses[i] > 1.e-8 )
		{
			forwardMomenta.segment<3>( 3 * i ) /= sqrt( m_gridMasses[i] );
		}
		else
		{
			forwardMomenta.segment<3>( 3 * i ).setZero();
		}
	}
	
	// this is broken for moving collision objects! Currently the implicit update matrix thinks they're not moving...
	implicitSolver(
		ImplicitUpdateMatrix( d, *this, collisionObjects ),
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


float Grid::updateDeformationGradients( ParticleData& d )
{

	ShapeFunction::PointToGridIterator& shIt = pointIterator();
	Vector3f weightGrad;
	Vector3i particleCell;
	Matrix3f delV;
	float maxMovement(0);

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
		
		float movement = 0.5 * ( delV + delV.transpose() ).norm();
		if( movement > maxMovement )
		{
			maxMovement = movement;
		}

		Matrix3f newParticleF = ( Matrix3f::Identity() + m_timeStep * delV ) * d.particleF[p];
		d.particleF[p] = newParticleF;
	}

	m_constitutiveModel.updateDeformation( d );
	return maxMovement;
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
