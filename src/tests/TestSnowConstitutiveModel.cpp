#include "tests/TestSnowConstitutiveModel.h"

#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/ParticleData.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

void testSnowConstitutiveModel()
{
	SnowConstitutiveModel snowModel(
		1, // young's modulus
		0.25, // poisson ratio
		10, // hardening
		0.1f, // compressive strength
		0.2f	// tensile strength
	);
	
	std::vector< Eigen::Vector3f > positions;
	std::vector< float > masses;
	positions.push_back( Eigen::Vector3f::Zero() );
	masses.push_back( 1.0f );

	ParticleData d( positions, masses, 1 );
	snowModel.initParticles( d );
	Matrix3f p = Matrix3f::Random();
	p = p * p;
	d.particleFplastic[0] = p;
	
	Matrix3f m = Matrix3f::Random();
	m = m * m;
	d.particleF[0] = m;
	
	snowModel.updateDeformation(d);
	
	// test general matrix decomposition properties:
	assert( ( d.particleF[0] * d.particleFplastic[0] - m * p ).norm() < 1.e-4 );
	assert( ( d.particleR[0] * d.particleS[0] - d.particleF[0] ).norm() < 1.e-4 );
	assert( ( d.particleF[0].inverse().transpose() - d.particleFinvTrans[0] ).norm() < 1.e-4 );
	assert( fabs( d.particleJ[0] - d.particleF[0].determinant() ) < 1.e-4 );
	
	// test stretch clamping:
	d.particleFplastic[0] = Matrix3f::Identity();
	d.particleF[0] = Matrix3f::Identity();
	d.particleF[0](0,0) = 1.3f;
	d.particleF[0](1,1) = 0.8f;
	
	snowModel.updateDeformation(d);
	
	assert( abs( d.particleF[0](0,0) - 1.2f ) < 1.e-4 );
	assert( abs( d.particleF[0](1,1) - 0.9f ) < 1.e-4 );
	
}

}
