#include "MpmSim/Sim.h"
#include "MpmSim/CollisionObject.h"
#include "MpmSim/ForceField.h"
#include "MpmSim/ConstitutiveModel.h"
#include "MpmSim/Grid.h"
#include "MpmSim/DiagonalPreconditioner.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/KdTree.h"

#include <iostream>
#include <queue>
#include <stdexcept>

using namespace MpmSim;
using namespace Eigen;


int Sim::CollisionObjectSet::collide(
	Vector3f& v,
	const Vector3f& x,
	const Vector3f& frameVelocity,
	bool addCollisionVelocity
) const
{
	int collisionObject(-1);
	bool intersectedAlready(false);
	for( size_t objIdx = 0; objIdx < objects.size(); ++objIdx )
	{
		float phi = objects[objIdx]->phi( x );
		if( phi <= 0 )
		{
			if( intersectedAlready )
			{
				// colliding with more than one object: set
				// velocity to zero and bail
				collisionObject = -2;
				v.setZero();
				break;
			}
			intersectedAlready = true;

			// intersecting the object
			Vector3f vObj;
			objects[objIdx]->velocity( x, vObj );
			
			// express object velocity relative to moving frame:
			vObj -= frameVelocity;

			// subtract off object velocity:
			v -= vObj;
			
			Vector3f n;
			objects[objIdx]->grad( x, n );
			n.normalize();
			
			float nDotV = n.dot( v );
			if( nDotV < 0 )
			{
				// trying to move into the object:
				collisionObject = (int)objIdx;

				if( objects[objIdx]->sticky() )
				{
					v.setZero();
				}
				else
				{
					
					// velocity perpendicular to the object
					Vector3f vPerp = nDotV * n;
					
					// remaining component is velocity paralell to the object:
					Vector3f vTangent = v - vPerp;
					float vtNorm = vTangent.norm();
					float coulombFriction = objects[objIdx]->coulombFriction();
					if( vtNorm >= -nDotV * coulombFriction )
					{
						v = vTangent * ( 1 + coulombFriction * nDotV / vTangent.norm() );
					}
					else
					{
						v.setZero();
					}
				}
			}
			
			if( addCollisionVelocity )
			{
				v += vObj;
			}
		}
	}
	
	return collisionObject;
}



Sim::Sim(
	const std::vector<Vector3f>& x,
	const std::vector<float>& masses,
	float gridSize,
	const ShapeFunction& shapeFunction,
	const ConstitutiveModel& model,
	const CollisionObjectSet& collisionObjects,
	const ForceFieldSet& forceFields,
	int dimension
) :
	m_gridSize( gridSize ),
	m_shapeFunction( shapeFunction ),
	m_constitutiveModel( model ),
	m_collisionObjects( collisionObjects ),
	m_forceFields( forceFields ),
	m_dimension( dimension )
{	
	particleData.variable<Vector3f>("p") = x;
	particleData.variable<Vector3f>("v").resize( x.size(), Vector3f::Zero() );
	particleData.variable<Matrix3f>("F").resize( x.size(), Matrix3f::Identity() );
	particleData.variable<float>("m") = masses;
	particleData.variable<float>("volume").resize( x.size(), 0.0f );
	
	m_constitutiveModel.createParticleData( particleData );
	m_constitutiveModel.setParticles( particleData );
	
	calculateBodies();
	
	for( std::vector< IndexList >::iterator it = m_bodies.begin(); it != m_bodies.end(); ++it )
	{
		IndexList& b = *it;
		Grid g( particleData, b, gridSize, shapeFunction, Eigen::Vector3f::Zero(),m_dimension );
		g.computeParticleVolumes();
	}
}

void Sim::advance( float timeStep, TerminationCriterion& termination, LinearSolver::Debug* d )
{
	std::vector<Eigen::Vector3f>& particleX = particleData.variable<Vector3f>( "p" );
	std::vector<Eigen::Vector3f>& particleV = particleData.variable<Vector3f>( "v" );
	std::vector<float>& particleMasses = particleData.variable<float>( "m" );
	
	// advance ballistic particle velocities:
	std::cerr << m_ballisticParticles.size() << " ballistic" << std::endl;
	IndexIterator ballisticEnd = m_ballisticParticles.end();
	for( IndexIterator it = m_ballisticParticles.begin(); it != ballisticEnd; ++it )
	{
		int p = *it;
		Eigen::Vector3f accn = Eigen::Vector3f::Zero();
		for( size_t i=0; i < m_forceFields.fields.size(); ++i )
		{
			accn += m_forceFields.fields[i]->force( particleX[p], particleMasses[p] );
		}
		accn /= particleMasses[p];
		particleV[p] += timeStep * accn;
	}
	// \todo: apply collisions. Or can I just apply them all at once at the end?
	
	// update velocities on particles in material point bodies:
	BodyIterator bodyEnd = m_bodies.end();
	std::cerr << m_bodies.size() << " bodies" << std::endl;
	for( BodyIterator bIt = m_bodies.begin(); bIt != bodyEnd; ++bIt )
	{
		// find centre of mass velocity, so we can make a comoving grid:
		// \todo: angular velocity sounds worthwhile (possibly more so than linear)
		// although more fiddly
		Eigen::Vector3f centreOfMassVelocity = Eigen::Vector3f::Zero();
		
		float mass = 0;
		std::cerr << "Body " << ( bIt - m_bodies.begin() ) << ": " << bIt->size() << " particles" << std::endl;
		for( IndexIterator it = bIt->begin(); it != bIt->end(); ++it )
		{
			centreOfMassVelocity += particleV[*it] * particleMasses[*it];
			mass += particleMasses[*it];
		}
		centreOfMassVelocity /= mass;
		
		// construct comoving background grid for this body:
		Grid g( particleData, *bIt, m_gridSize, m_shapeFunction, centreOfMassVelocity, m_dimension );
		
		// update grid velocities using internal stresses...
		DiagonalPreconditioner preconditioner( g, m_constitutiveModel, timeStep );
		ConjugateResiduals solver( termination, &preconditioner );

		//ConjugateResiduals solver( termination );
		g.updateGridVelocities(
			timeStep,
			m_constitutiveModel,
			m_collisionObjects,
			m_forceFields.fields,
			solver,
			d
		);
		
		// transfer the grid velocities back onto the particles:
		g.updateParticleVelocities();
		// \todo: apply collisions. Or can I just apply them all at once at the end?
		
		// update particle deformation gradients:
		g.updateDeformationGradients( timeStep );
		m_constitutiveModel.updateParticleData( particleData );
		
	}
	
	// advance particle positions:
	std::vector<Eigen::Vector3f>::iterator end = particleX.end();
	std::vector<Eigen::Vector3f>::iterator it = particleX.begin();
	std::vector<Eigen::Vector3f>::iterator vIt = particleV.begin();
	for( ; it != end; ++it, ++vIt )
	{
		// resolve object collisions:
		m_collisionObjects.collide( *vIt, *it, Eigen::Vector3f::Zero(), true );
		*it += *vIt * timeStep;
	}

	calculateBodies();
}

size_t Sim::numParticleVariables() const
{
	return particleData.numVariables();
}


// used for sorting a list of points by their positions on the "dim" axis
class SpatialComparator
{
public:
	
	SpatialComparator( int dim, const std::vector< Eigen::Vector3f >& x ) :
		m_dim( dim ), m_x( x )
	{
	}
	
	bool operator()( int a, int b )
	{
		return m_x[a][m_dim] < m_x[b][m_dim];
	}

private:
	int m_dim;
	const std::vector< Eigen::Vector3f >& m_x;
};

void Sim::voxelSort(
	IndexIterator begin,
	IndexIterator end,
	float voxelSize,
	const std::vector<Eigen::Vector3f>& particleX,
	int dim )
{
	// sort along dim:
	std::sort( begin, end, SpatialComparator( dim, particleX ) );

	if( dim == 2 )
	{
		// got to z coord - done!
		return;
	}
	
	// now chop into slices along dim, and sort each slice along dim + 1:
	int currentSlice = int( floor( particleX[*begin][dim] / voxelSize ) );
	std::vector<int>::iterator sliceBegin = begin;
	std::vector<int>::iterator sliceEnd = begin + 1;
	for( ;;++sliceEnd )
	{
		int particleSlice = sliceEnd == end ? currentSlice + 1 : int( floor( particleX[*sliceEnd][dim] / voxelSize ) );
		if( particleSlice != currentSlice )
		{
			voxelSort( sliceBegin, sliceEnd, voxelSize, particleX, dim + 1 );
			if( sliceEnd == end )
			{
				break;
			}
			sliceBegin = sliceEnd;
			currentSlice = particleSlice;
		}
	}
}

void Sim::calculateBodies()
{

	const std::vector<Eigen::Vector3f>& particleV = particleData.variable<Vector3f>("v");
	std::vector<Eigen::Vector3f>& particleX = particleData.variable<Vector3f>("p");
	
	m_bodies.clear();
	m_ballisticParticles.clear();

	// build little grid based acceleration structure for neighbour queries:
	V3fTree n( particleX.begin(), particleX.end() );
	std::vector<bool> processed( particleX.size(), false );
	std::vector< std::vector<Eigen::Vector3f>::iterator > nearNeighbours;
	for( size_t i=0; i < particleX.size(); ++i )
	{
		float prod = particleV[i][0] * particleV[i][1] * particleV[i][2];
		#ifdef WIN32
		if( !_finite(prod) )
		#else
		if( isinff(prod) || isnanf(prod) )
		#endif
		{
			throw std::runtime_error( "nans in particle velocity data!" );
		}
		
		if( processed[i] )
		{
			continue;
		}
		n.nearestNeighbours( particleX[i], m_gridSize, nearNeighbours );
		if( nearNeighbours.size() == 1 )
		{
			// this means the only particle within a radius of m_gridSize of
			// particleX[i] IS particleX[i], so stick it on the ballistic list:
			processed[i] = true;
			m_ballisticParticles.push_back((int)i);
			continue;
		}
		
		// looks like the particle actually has neighbours: create a body
		m_bodies.resize( m_bodies.size() + 1 );
		IndexList& b = m_bodies.back();

		// now push all the neighbours onto a queue and do a floodfill to fill up the body:
		std::queue<int> flood;
		for( size_t j=0; j < nearNeighbours.size(); ++j )
		{
			flood.push((int)(nearNeighbours[j] - particleX.begin()));
		}
		while( flood.size() )
		{
			int current = flood.front();
			flood.pop();
			if( processed[current] )
			{
				continue;
			}

			b.push_back( current );
			processed[current] = true;
			n.nearestNeighbours( particleX[current], m_gridSize, nearNeighbours );
			for( size_t j=0; j < nearNeighbours.size(); ++j )
			{
				flood.push( (int)(nearNeighbours[j] - particleX.begin()));
			}
		}

		// sort the spatial index so particles in the same voxel are adjacent:
		voxelSort( b.begin(), b.end(), 2 * m_shapeFunction.supportRadius() * m_gridSize, particleX );
	}
}

Sim::CollisionObjectSet::~CollisionObjectSet()
{
	for( size_t i=0; i < objects.size(); ++i )
	{
		delete objects[i];
	}
}

void Sim::CollisionObjectSet::add( CollisionObject* o )
{
	objects.push_back( o );
}

Sim::ForceFieldSet::~ForceFieldSet()
{
	for( size_t i=0; i < fields.size(); ++i )
	{
		delete fields[i];
	}
}

void Sim::ForceFieldSet::add( ForceField* f )
{
	fields.push_back( f );
}

size_t Sim::numBodies() const
{
	return m_bodies.size();
}

// particle indices in body n:
const Sim::IndexList& Sim::body( size_t n ) const
{
	if( n >= m_bodies.size() )
	{
		throw std::runtime_error( "Sim::body(): index exceeds number of bodies" );
	}
	return m_bodies[n];
}

// particle indices for the ballistic particles:
const Sim::IndexList& Sim::ballisticParticles() const
{
	return m_ballisticParticles;
}

