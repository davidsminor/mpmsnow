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
	
	// set deformation within yeild limits:
	d.particleF[0](0,0) = 1.02f;
	d.particleF[0](1,1) = 0.95f;
	Matrix3f originalF = d.particleF[0];
	snowModel.updateDeformation(d);
	
	// work out dedf with model:
	Matrix3f dEdF;
	snowModel.dEnergyDensitydF( dEdF, d, 0 );

	// compare to finite difference version:
	const float delta = 0.01f;
	Matrix3f dEdF_fd;
	for( int i=0; i < 3; ++i )
	{
		for( int j=0; j < 3; ++j )
		{
			d.particleF[0] = originalF;
			d.particleF[0](i,j) += delta;
			snowModel.updateDeformation(d);
			float ePlusDeltaE = snowModel.energyDensity( d, 0 );
			
			d.particleF[0] = originalF;
			d.particleF[0](i,j) -= delta;
			snowModel.updateDeformation(d);
			float eMinusDeltaE = snowModel.energyDensity( d, 0 );

			dEdF_fd( i, j ) = ( ePlusDeltaE - eMinusDeltaE ) / ( 2 * delta );
		}
	}
	assert( ( dEdF - dEdF_fd ).maxCoeff() - ( dEdF - dEdF_fd ).minCoeff() < 1.e-6 );
	
	// now work out "force differential density" (maybe there's a better name for that...)
	d.particleF[0] = originalF;
	snowModel.updateDeformation(d);
	
	Matrix3f dF = Matrix3f::Random() * 0.001f;
	Matrix3f dEdFDifferential;
	snowModel.dEdFDifferential( dEdFDifferential, dF, d, 0 );
	
	// compare with finite difference:
	Matrix3f plusDelta;
	d.particleF[0] += dF;
	snowModel.updateDeformation(d);
	snowModel.dEnergyDensitydF( plusDelta, d, 0 );
	
	assert( ( plusDelta - dEdF - dEdFDifferential ).maxCoeff() - ( plusDelta - dEdF - dEdFDifferential ).minCoeff() < 1.e-5 );
}

}
