#include "tests/TestSimClass.h"

#include "MpmSim/Sim.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/ConstitutiveModel.h"
#include "MpmSim/SnowConstitutiveModel.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{


class DummyModel : public ConstitutiveModel
{
public:

	virtual void createParticleData( Sim::MaterialPointDataMap& p ) const
	{}
	
	virtual void updateParticleData( Sim& sim ) const
	{}

	virtual void setParticles( const Sim::MaterialPointDataMap& p ) const
	{}
	
	virtual float energyDensity( size_t p ) const
	{ return 0; }
	
	virtual Eigen::Matrix3f dEnergyDensitydF( size_t p ) const
	{ return Eigen::Matrix3f::Zero(); }
	
	virtual Eigen::Matrix3f dEdFDifferential( const Eigen::Matrix3f& dFp, size_t p ) const
	{ return Eigen::Matrix3f::Zero(); }


};

void testNeighbourQuery()
{
	srand(20);
	std::vector<Eigen::Vector3f> particleX;
	for( int i=0; i < 1000; ++i )
	{
		float x = ((float)rand() / RAND_MAX);
		float y = ((float)rand() / RAND_MAX);
		float z = ((float)rand() / RAND_MAX);
		particleX.push_back( Eigen::Vector3f(x,y,z) );
	}

	Sim::NeighbourQuery n( particleX, 0.1f );
	for( size_t i=0; i < particleX.size(); ++i )
	{
		const Eigen::Vector3f& p = particleX[i];
		std::vector<int> actualNeighbours;
		for( size_t j=0; j < particleX.size(); ++j )
		{
			float squaredDist = ( particleX[j] - p ).squaredNorm();
			if( squaredDist > 0 && squaredDist < ( 0.1f * 0.1f ) )
			{
				actualNeighbours.push_back( (int)j );
			}
		}
		std::vector<int> neighbours;
		
		n.neighbours( p, neighbours );
		assert( neighbours.size() == actualNeighbours.size() );
		
		std::sort( neighbours.begin(), neighbours.end() );
		std::sort( actualNeighbours.begin(), actualNeighbours.end() );
		
		std::vector<int>::iterator it = neighbours.begin();
		std::vector<int>::iterator ait = actualNeighbours.begin();
		for( ;it != neighbours.end(); ++it, ++ait )
		{
			assert( *it == *ait );
		}
	}
}

void testInitialization()
{

	std::vector<Vector3f> positions;
	std::vector<float> masses;
	const float gridSize = 0.1f;
	
	// create two cube shaped clusters of particles.
	// particle spacing = 1/2 of a cell:
	for( int i=0; i < 16; ++i )
	{
		for( int j=0; j < 16; ++j )
		{
			for( int k=0; k < 16; ++k )
			{
				positions.push_back( Vector3f( 0.5f * gridSize * (i+0.5f), 0.5f * gridSize * (j+0.5f), 0.5f * gridSize * (k+0.5f) ) );
				masses.push_back( 1.0f );

				positions.push_back( Vector3f( 0.5f * gridSize * (i+0.5f) + 6, 0.5f * gridSize * (j+0.5f), 0.5f * gridSize * (k+0.5f) ) );
				masses.push_back( 1.0f );

			}
		}
	}

	// create a bunch of ballistic particles:
	for( int i=0; i < 16; ++i )
	{
		for( int j=0; j < 16; ++j )
		{
			for( int k=0; k < 16; ++k )
			{
				positions.push_back( Vector3f( 1.1f * gridSize * i + 20, 1.1f * gridSize * j, 1.1f * gridSize * k ) );
				masses.push_back( 1.0f );

			}
		}
	}
	
	CubicBsplineShapeFunction shapeFunction;
	DummyModel constitutiveModel;
	Sim::CollisionObjectSet collisionObjects;
	Sim::ForceFieldSet forceFields;
	Sim sim( positions, masses, gridSize, shapeFunction, constitutiveModel, collisionObjects, forceFields );
	
	// sensible particle variables?
	assert( sim.numParticleVariables() == 5 );
	assert( sim.particleVariable<MpmSim::ScalarData>("m") );
	assert( sim.particleVariable<MpmSim::ScalarData>("volume") );
	assert( sim.particleVariable<MpmSim::VectorData>("v") );
	assert( sim.particleVariable<MpmSim::VectorData>("p") );
	assert( sim.particleVariable<MpmSim::MatrixData>("F") );
	
	assert( sim.particleVariable<MpmSim::ScalarData>("F") == 0 );
	assert( sim.particleVariable<MpmSim::ScalarData>("Fiddlesticks") == 0 );

	assert( sim.ballisticParticles().size() == 16 * 16 * 16 );
	
	assert( sim.numBodies() == 2 );
	assert( sim.body( 0 ).size() == 16 * 16 * 16 );
	assert( sim.body( 1 ).size() == 16 * 16 * 16 );
}

void testSimClass()
{
	testNeighbourQuery();
	testInitialization();
}

}
