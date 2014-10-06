#include <Eigen/Geometry>

#include "tests/TestSnowConstitutiveModel.h"

#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/Sim.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

void TestSnowConstitutiveModel::test()
{
	std::cerr << "testSnowConstitutiveModel()" << std::endl;
	
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	
	positions.push_back( Eigen::Vector3f::Zero() );
	masses.push_back( 1.0f );
	
	const float gridSize(1);
	CubicBsplineShapeFunction shapeFunction;
	
	SnowConstitutiveModel snowModel(
		1, // young's modulus
		0.25, // poisson ratio
		10, // hardening
		0.1f, // compressive strength
		0.2f	// tensile strength
	);

	Sim::CollisionObjectSet collisionObjects;
	Sim::ForceFieldSet forceFields;
	Sim sim( positions, masses, gridSize, shapeFunction, snowModel, collisionObjects, forceFields );
	Matrix3f p = Matrix3f::Random();
	p = p * p;
	
	// should have created all this particle data for the sim:
	std::vector<Eigen::Matrix3f>& particleF = sim.particleData.variable<Matrix3f>( "F" );
	std::vector<Eigen::Matrix3f>& particleFplastic = sim.particleData.variable<Matrix3f>( "Fp" );
	std::vector<Eigen::Matrix3f>& particleFinvTrans = sim.particleData.variable<Matrix3f>( "FinvTrans" );
	std::vector<Eigen::Matrix3f>& particleR = sim.particleData.variable<Matrix3f>( "R" );
	std::vector<Eigen::Matrix3f>& particleGinv = sim.particleData.variable<Matrix3f>( "Ginv" );
	std::vector<float>& particleJ = sim.particleData.variable<float>( "J" );
	std::vector<float>& particleMu = sim.particleData.variable<float>( "mu" );
	std::vector<float>& particleLambda = sim.particleData.variable<float>( "lambda" );
	
	particleFplastic[0] = p;
	
	Matrix3f m = Matrix3f::Random();
	m = m * m;
	particleF[0] = m;
	
	snowModel.updateParticleData( sim.particleData );
	
	// test general matrix decomposition properties:
	assert( ( particleF[0] * particleFplastic[0] - m * p ).norm() < 1.e-4 );
	assert( ( particleF[0].inverse().transpose() - particleFinvTrans[0] ).norm() < 1.e-4 );
	assert( fabs( particleJ[0] - particleF[0].determinant() ) < 1.e-4 );
	
	// test stretch clamping:
	particleFplastic[0] = Matrix3f::Identity();
	particleF[0] = Matrix3f::Identity();
	particleF[0](0,0) = 1.3f;
	particleF[0](1,1) = 0.8f;
	
	snowModel.updateParticleData( sim.particleData );
	
	assert( abs( particleF[0](0,0) - 1.2f ) < 1.e-4 );
	assert( abs( particleF[0](1,1) - 0.9f ) < 1.e-4 );
	
	// identity particle should have zero enegy density and forces:
	particleF[0] = Matrix3f::Identity();
	snowModel.updateParticleData( sim.particleData );
	assert( snowModel.energyDensity( 0 ) == 0 );
	assert( snowModel.dEnergyDensitydF( 0 ) == Matrix3f::Zero() );
	
	// now test a purely rotated particle:
	particleF[0] = AngleAxisf(float(0.25*M_PI), Vector3f::UnitX())
	  * AngleAxisf(float(0.5*M_PI),  Vector3f::UnitY())
	  * AngleAxisf(float(0.33*M_PI), Vector3f::UnitZ());
	snowModel.updateParticleData( sim.particleData );

	// should also have zero energy density and forces, innit!
	assert( fabs( snowModel.energyDensity( 0 ) ) < 1.e-6 );
	assert( snowModel.dEnergyDensitydF( 0 ).norm() < 1.e-6 );

	
	// set deformation within yeild limits:
	particleF[0] = Matrix3f::Identity();
	particleF[0](0,0) = 1.02f;
	particleF[0](1,1) = 0.95f;

	// rotate its ass:
	particleF[0] = particleF[0] * AngleAxisf(float(0.25*M_PI), Vector3f::UnitX())
	  * AngleAxisf(float(0.5*M_PI),  Vector3f::UnitY())
	  * AngleAxisf(float(0.33*M_PI), Vector3f::UnitZ());

	Matrix3f originalF = particleF[0];
	snowModel.updateParticleData( sim.particleData );
	
	// work out dedf with model:
	snowModel.setParticles( sim.particleData );
	Matrix3f dEdF = snowModel.dEnergyDensitydF( 0 );

	// compare to finite difference version:
	const float delta = 0.01f;
	Matrix3f dEdF_fd;
	for( int i=0; i < 3; ++i )
	{
		for( int j=0; j < 3; ++j )
		{
			particleF[0] = originalF;
			particleF[0](i,j) += delta;
			snowModel.updateParticleData( sim.particleData );
			float ePlusDeltaE = snowModel.energyDensity( 0 );
			
			particleF[0] = originalF;
			particleF[0](i,j) -= delta;
			snowModel.updateParticleData( sim.particleData );
			float eMinusDeltaE = snowModel.energyDensity( 0 );

			dEdF_fd( i, j ) = ( ePlusDeltaE - eMinusDeltaE ) / ( 2 * delta );
		}
	}
	float diffMaxCoeff = fabs( ( dEdF - dEdF_fd ).maxCoeff() );
	float diffMinCoeff = fabs( ( dEdF - dEdF_fd ).minCoeff() );
	assert( diffMaxCoeff < 1.e-5 );
	assert( diffMinCoeff < 1.e-5 );
	
	// now work out "force differential density" (maybe there's a better name for that...)
	particleF[0] = originalF;
	snowModel.updateParticleData( sim.particleData );
	
	Matrix3f dF = Matrix3f::Random() * 0.001f;
	Matrix3f dEdFDifferential = snowModel.dEdFDifferential( dF, 0 );
	
	// compare with finite difference:
	particleF[0] += dF;
	snowModel.updateParticleData( sim.particleData );
	Matrix3f plusDelta = snowModel.dEnergyDensitydF( 0 );
	
	Matrix3f dEdFDifferential_fd = plusDelta - dEdF;
	diffMaxCoeff = fabs( ( dEdFDifferential_fd - dEdFDifferential ).maxCoeff() );
	diffMinCoeff = fabs( ( dEdFDifferential_fd - dEdFDifferential ).minCoeff() );
	
	std::cerr << dEdFDifferential << std::endl;
	std::cerr << dEdFDifferential_fd << std::endl;

	assert( diffMaxCoeff < 1.e-5 );
	assert( diffMinCoeff < 1.e-5 );

}

}
