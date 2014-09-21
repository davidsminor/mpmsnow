#include "tests/TestSimClass.h"

#include "MpmSim/Sim.h"
#include "MpmSim/GravityField.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/ConstitutiveModel.h"
#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/SquareMagnitudeTermination.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{


class DummyModel : public ConstitutiveModel
{
public:

	virtual void createParticleData( MaterialPointData& p ) const
	{}
	
	virtual void updateParticleData( MaterialPointData& p ) const
	{}

	virtual void setParticles( MaterialPointData& p ) const
	{}
	
	virtual float energyDensity( size_t p ) const
	{ return 0; }
	
	virtual Eigen::Matrix3f dEnergyDensitydF( size_t p ) const
	{ return Eigen::Matrix3f::Zero(); }
	
	virtual Eigen::Matrix3f dEdFDifferential( const Eigen::Matrix3f& dFp, size_t p ) const
	{ return Eigen::Matrix3f::Zero(); }


};

void testInitialization()
{
	std::cerr << "testInitialization()" << std::endl;

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

	// these throw exceptions if the requested items don't exist:
	sim.particleData.variable<float>("m");
	sim.particleData.variable<float>("volume");
	sim.particleData.variable<Eigen::Vector3f>("v");
	sim.particleData.variable<Eigen::Vector3f>("p");
	sim.particleData.variable<Eigen::Matrix3f>("F");
	
	bool exceptionThrown(false);
	try
	{
		sim.particleData.variable<float>("Fiddlesticks");
	}
	catch( ... )
	{
		exceptionThrown = true;
	}
	assert( exceptionThrown );

	int numBallistics = sim.ballisticParticles().size();
	assert( numBallistics == 16 * 16 * 16 );
	
	assert( sim.numBodies() == 2 );
	assert( sim.body( 0 ).size() == 16 * 16 * 16 );
	assert( sim.body( 1 ).size() == 16 * 16 * 16 );
}

void testTimestepAdvance()
{
	std::cerr << "testTimestepAdvance()" << std::endl;

	std::vector<Vector3f> positions;
	std::vector<float> masses;
	const float gridSize = 0.1f;
	
	// create two cube shaped clusters of particles.
	// particle spacing = 1/2 of a cell:
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				positions.push_back( Vector3f( 0.5f * gridSize * (i+0.5f), 0.5f * gridSize * (j+0.5f), 0.5f * gridSize * (k+0.5f) ) );
				masses.push_back( 1.0f );

				positions.push_back( Vector3f( 0.5f * gridSize * (i+0.5f) + 6, 0.5f * gridSize * (j+0.5f), 0.5f * gridSize * (k+0.5f) ) );
				masses.push_back( 1.0f );

			}
		}
	}

	// create a bunch of ballistic particles:
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				positions.push_back( Vector3f( 1.1f * gridSize * i + 20, 1.1f * gridSize * j, 1.1f * gridSize * k ) );
				masses.push_back( 1.0f );
			}
		}
	}
	
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel constitutiveModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);
	Sim::CollisionObjectSet collisionObjects;
	Sim::ForceFieldSet forceFields;
	forceFields.fields.push_back( new GravityField( Eigen::Vector3f( 0, 1.f, 0 ) ) );
	Sim sim( positions, masses, gridSize, shapeFunction, constitutiveModel, collisionObjects, forceFields );
	
	SquareMagnitudeTermination t( 10, 0.0f );
	sim.advance( 0.01f, t );
	sim.advance( 0.01f, t );
	sim.advance( 0.01f, t );
	
	// average velocity should be about 0.03 now, innit
	const std::vector<Eigen::Vector3f>& velocities = sim.particleData.variable<Eigen::Vector3f>( "v" );
	Eigen::Vector3f v = Eigen::Vector3f::Zero();
	for( std::vector<Eigen::Vector3f>::const_iterator it = velocities.begin(); it != velocities.end(); ++it )
	{
		v += *it;
	}
	v /= velocities.size();
	
	assert( fabs( v[0] ) < 1.e-6 );
	assert( fabs( v[2] ) < 1.e-6 );
	assert( fabs( v[1] - 0.03f ) < 1.e-5 );

}

void testSimClass()
{
	std::cerr << "testSimClass()" << std::endl;
	testInitialization();
	testTimestepAdvance();
}

}
