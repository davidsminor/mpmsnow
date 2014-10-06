#include "MpmSim/Sim.h"
#include "MpmSim/CollisionObject.h"
#include "MpmSim/ForceField.h"
#include "MpmSim/ConstitutiveModel.h"
#include "MpmSim/Grid.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/KdTree.h"

#include <iostream>
#include <queue>
#include <stdexcept>

using namespace MpmSim;
using namespace Eigen;


Sim::Sim(
	const std::vector<Vector3f>& x,
	const std::vector<float>& masses,
	float gridSize,
	const ShapeFunction& shapeFunction,
	ConstitutiveModel& model,
	const CollisionObject::CollisionObjectSet& collisionObjects,
	const ForceField::ForceFieldSet& forceFields,
	int dimension
) :
	m_gridSize( gridSize ),
	m_shapeFunction( shapeFunction ),
	m_constitutiveModel( model ),
	m_collisionObjects( collisionObjects ),
	m_forceFields( forceFields ),
	m_dimension( dimension )
{	
	m_particleData.variable<Vector3f>("p") = x;
	m_particleData.variable<Vector3f>("v").resize( x.size(), Vector3f::Zero() );
	m_particleData.variable<Matrix3f>("F").resize( x.size(), Matrix3f::Identity() );
	m_particleData.variable<float>("m") = masses;
	m_particleData.variable<float>("volume").resize( x.size(), 0.0f );
	
	m_constitutiveModel.setParticles( m_particleData );
	
	calculateBodies();
	
	for( std::vector< IndexList >::iterator it = m_bodies.begin(); it != m_bodies.end(); ++it )
	{
		IndexList& b = *it;
		Grid g( m_particleData, b, gridSize, shapeFunction, Eigen::Vector3f::Zero(),m_dimension );
		g.computeParticleVolumes();
	}
}

MaterialPointData& Sim::particleData()
{
	return m_particleData;
}

void Sim::advance( float timeStep, TerminationCriterion& termination, LinearSolver::Debug* d )
{
	std::vector<Eigen::Vector3f>& particleX = m_particleData.variable<Vector3f>( "p" );
	std::vector<Eigen::Vector3f>& particleV = m_particleData.variable<Vector3f>( "v" );
	std::vector<float>& particleMasses = m_particleData.variable<float>( "m" );
	
	// advance ballistic particle velocities:
	std::cerr << m_ballisticParticles.size() << " ballistic" << std::endl;
	IndexIterator ballisticEnd = m_ballisticParticles.end();
	for( IndexIterator it = m_ballisticParticles.begin(); it != ballisticEnd; ++it )
	{
		int p = *it;
		Eigen::Vector3f accn = Eigen::Vector3f::Zero();
		m_forceFields.force( accn, particleX[p], particleMasses[p] );
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
		Grid g( m_particleData, *bIt, m_gridSize, m_shapeFunction, centreOfMassVelocity, m_dimension );
		
		// update grid velocities using internal stresses...
		g.updateGridVelocities(
			timeStep,
			m_constitutiveModel,
			m_collisionObjects,
			m_forceFields,
			termination,
			d
		);
		if( termination.cancelled() )
		{
			calculateBodies();
			return;
		}

		// transfer the grid velocities back onto the particles:
		g.updateParticleVelocities();
		
		// update particle deformation gradients:
		g.updateDeformationGradients( timeStep );
		m_constitutiveModel.updateParticleData();
		
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

void Sim::calculateBodies()
{

	const std::vector<Eigen::Vector3f>& particleV = m_particleData.variable<Vector3f>("v");
	std::vector<Eigen::Vector3f>& particleX = m_particleData.variable<Vector3f>("p");
	
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
	}
}
