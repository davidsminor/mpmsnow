#include "tests/TestGrid.h"

#include "MpmSim/Grid.h"
#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/ConjugateGradients.h"
#include "MpmSim/CubicBsplineShapeFunction.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

class TestConstitutiveModel : public ConstitutiveModel
{
public:
	
	TestConstitutiveModel() {}

	virtual void initParticles( ParticleData& p ) const {};

	// update deformation at particle p:
	virtual void updateDeformation( ParticleData& d ) const {};

	// energy density for particle p:
	virtual float energyDensity( const ParticleData& d, size_t p ) const
	{
		return 0;
	}
	
	// derivative of energy density with respect to the deformation gradient at particle p:
	virtual void dEnergyDensitydF( Eigen::Matrix3f& result, const ParticleData& d, size_t p ) const
	{
	}
	
	// energy differentials, using second derivatives of energy function
	virtual void dEdFDifferential( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const
	{
	}

};

static void testMassSplatting()
{
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	positions.push_back( Vector3f( 0, 0, 0 ) );
	masses.push_back( 1.0f );
	for( int i=0; i < 5000; ++i )
	{
		float xr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float yr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float zr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
		positions.push_back( Vector3f( xr, yr, zr ) );
		masses.push_back( 1.0f );
	}
	
	const float gridSize = 0.01f;
	ParticleData d( positions, masses, gridSize );
	
	// create a grid, which will immediately splat the particle masses onto itself,
	// multi threaded tbb stylee:
	CubicBsplineShapeFunction shapeFunction;
	TestConstitutiveModel testConstitutiveModel;
	Grid g( d, 0.1f, shapeFunction, testConstitutiveModel );
	
	// now manually splat the masses, serial stylee:
	Eigen::VectorXf massVector( g.masses().size() );
	massVector.setZero();
	
	ShapeFunction::PointToGridIterator& shIt = g.pointIterator();
	Vector3i particleCell;
	for( size_t p = 0; p < d.particleX.size(); ++p )
	{
		shIt.initialize( d.particleX[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = g.coordsToIndex( particleCell );
			massVector[ idx ] += d.particleM[p] * shIt.w();
		} while( shIt.next() );
	}	
	
	massVector -= g.masses();
	assert( massVector.maxCoeff() - massVector.minCoeff() < 1.e-6 );

}

static void testDeformationGradients()
{
	// unit cube full of particles centered at the origin:
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	positions.push_back( Vector3f( 0, 0, 0 ) );
	masses.push_back( 1.0f );
	for( int i=0; i < 10; ++i )
	{
		for( int j=0; j < 10; ++j )
		{
			for( int k=0; k < 10; ++k )
			{
				positions.push_back( Vector3f( float(i) / 9 - 0.5, float(j) / 9 - 0.5, float(k) / 9 - 0.5 ) );
				masses.push_back( 1.0f );
			}
		}
	}
	
	// now set the velocities so it grows in the xz plane and shrinks in the y axis:
	const float gridSize = 0.2f;
	ParticleData d( positions, masses, gridSize );
	
	for( size_t p=0; p < d.particleX.size(); ++p )
	{
		d.particleV[p][0] = d.particleX[p][0];
		d.particleV[p][1] = -d.particleX[p][1];
		d.particleV[p][2] = 0.5f * d.particleX[p][2];
	}
	
	// advance simulation 10 steps without internal forces, and compare the computed deformation
	// gradients with their analytic values:
	CubicBsplineShapeFunction shapeFunction;
	TestConstitutiveModel testConstitutiveModel;
	
	for( int i=1; i < 8; ++i )
	{
		Grid g( d, 0.01f, shapeFunction, testConstitutiveModel );
		g.updateDeformationGradients( d );
		d.advance( 0.01f );
		assert( fabs( d.particleF[0]( 0,0 ) - ( 1 + 0.01f * i ) ) < 0.001f );
		assert( fabs( d.particleF[0]( 1,1 ) - ( 1 - 0.01f * i ) ) < 0.001f );
		assert( fabs( d.particleF[0]( 2,2 ) - ( 1 + 0.005f * i ) ) < 0.001f );
	}

	// NB: significant errors get introduced to the deformation gradient round the edges of
	// the object it seems, which is probably responsible for some visual artefacts... Maybe
	// it's worth fixing these some time?

}



static void testForces()
{
	// create some particles:
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	positions.push_back( Vector3f( 0.1f,0.2f,0.f ) );
	masses.push_back( 1.0f );
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				positions.push_back( Vector3f( float( i ) / 3 - 0.5f, float( j ) / 3 - 0.5f, float( k ) / 3 - 0.5f ) );
				masses.push_back( 1.0f );
			}
		}
	}

	// create particle data:
	const float gridSize = 0.5f;
	ParticleData d( positions, masses, gridSize );
	
	// give it a sinusoidal velocity field and displace it a bit:
	for( size_t p=0; p < d.particleX.size(); ++p )
	{
		d.particleV[p][0] = 100 * sin( 6 * d.particleX[p][0] );
		d.particleV[p][1] = 100 * sin( 6 * d.particleX[p][1] );
		d.particleV[p][2] = 100 * sin( 6 * d.particleX[p][2] );
	}
	
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel snowModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);
	snowModel.initParticles( d );
	
	const float timeStep = 0.001f;
	Grid g( d, timeStep, shapeFunction, snowModel );

	g.computeDensities( d );
	d.particleVolumes.resize( d.particleX.size() );
	for( size_t p = 0; p < d.particleDensities.size(); ++p )
	{
		d.particleVolumes[p] = d.particleM[p] / d.particleDensities[p];
	}

	g.updateDeformationGradients( d );
	VectorXf gridVelocities = g.getVelocities();

	// calculate da forces brah!
	VectorXf forces = VectorXf::Zero( gridVelocities.size() );
	g.calculateForces( d, forces );
	

	// calculate unperturbed energy:
	float e0 = g.calculateEnergy( d );
	
	// sanity check: this shouldn't affect ANYTHING:
	ParticleData dTest = d;

	gridVelocities.setZero();
	g.setVelocities( gridVelocities );
	g.updateDeformationGradients( dTest );
	
	for( size_t p=0; p < d.particleF.size(); ++p )
	{
		assert( d.particleF[p] == dTest.particleF[p] );
		assert( d.particleR[p] == dTest.particleR[p] );
		assert( d.particleJ[p] == dTest.particleJ[p] );
		assert( d.particleMu[p] == dTest.particleMu[p] );
		assert( d.particleLambda[p] == dTest.particleLambda[p] );
		assert( d.particleVolumes[p] == dTest.particleVolumes[p] );
	}
	
	assert( e0 == g.calculateEnergy( dTest ) );
	
	const VectorXf& m = g.masses();

	// now we're gonna calculate energy derivatives... the stupid way!
	// we're gonna do this component by component, and we're gonna do it
	// by zeroing out the grid velocities, setting the component we're gonna
	// test to delta/m_timeStep, advancing bits of the sim with that velocity field,
	// calculating the energy in the final state (in which one of the grid nodes
	// will have moved a distance delta along one of the axes), and using the result
	// to calculate a finite difference derivative!
	float delta = 0.05f;
	float squaredError(0);
	float squaredForce(0);
	int count(0);
	for( int idx = 0; idx < m.size(); ++idx )
	{
		if( m[idx] < 1 )
		{
			continue;
		}
		for( size_t dim = 0; dim < 3; ++dim )
		{
			gridVelocities.setZero();
			
			// perturb current grid point a distance delta along the current axis,
			// and calculate the resulting deformation gradients:
			dTest = d;
			gridVelocities( 3 * idx + dim ) = delta / timeStep;
			g.setVelocities( gridVelocities );
			g.updateDeformationGradients( dTest );
			float ePlus = g.calculateEnergy( dTest );

			dTest = d;
			gridVelocities( 3 * idx + dim ) = -delta / timeStep;
			g.setVelocities( gridVelocities );
			g.updateDeformationGradients( dTest );
			float eMinus = g.calculateEnergy( dTest );
			
			// so force = -dE/dX = ( e0 - e ) / delta
			float f = ( eMinus - ePlus ) / ( 2 * delta );
			squaredError += ( f - forces( 3 * idx + dim ) ) * ( f - forces( 3 * idx + dim ) );
			squaredForce += forces( 3 * idx + dim ) * forces( 3 * idx + dim );
			++count;
			std::cerr << f << " == " << forces( 3 * idx + dim ) << "?  " << (3 * idx + dim) << " of " << forces.size() << std::endl;
		}
	}
	
	assert( squaredError / squaredForce < 1.e-4 );
	
}

void testImplicitUpdate()
{
	// create some particles:
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	positions.push_back( Vector3f( 0.1f,0.2f,0.f ) );
	masses.push_back( 1.0f );
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				positions.push_back( Vector3f( float( i ) / 3 - 0.5f, float( j ) / 3 - 0.5f, float( k ) / 3 - 0.5f ) );
				masses.push_back( 1.0f );
			}
		}
	}

	// create particle data:
	const float gridSize = 0.5f;
	ParticleData d( positions, masses, gridSize );
	
	// give it a sinusoidal velocity field and displace it a bit:
	for( size_t p=0; p < d.particleX.size(); ++p )
	{
		d.particleV[p][0] = 0.01f*sin( 6 * d.particleX[p][0] );
		d.particleV[p][1] = 0.01f*sin( 6 * d.particleX[p][1] );
		d.particleV[p][2] = 0.01f*sin( 6 * d.particleX[p][2] );
	}
	
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel snowModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);
	snowModel.initParticles( d );
	
	const float timeStep = 0.005f;
	Grid g( d, timeStep, shapeFunction, snowModel );

	g.computeDensities( d );
	d.particleVolumes.resize( d.particleX.size() );
	for( size_t p = 0; p < d.particleDensities.size(); ++p )
	{
		d.particleVolumes[p] = d.particleM[p] / d.particleDensities[p];
	}
	
	// save grid velocities:
	VectorXf initialGridVelocities( g.getVelocities().size() );
	for( int i=0; i < g.getVelocities().size(); ++i )
	{
		initialGridVelocities[i] = g.getVelocities()[i];
	}

	// do an implicit update:
	std::vector<CollisionObject*> collisionObjects;
	g.updateGridVelocities( d, collisionObjects, ConjugateResiduals( 500, 1.e-7 ) );

	// transfer the grid velocities back onto the particles:
	g.updateParticleVelocities( d, collisionObjects );
	
	// update particle deformation gradients:
	g.updateDeformationGradients( d );
	
	// calculate da forces brah, and do a backwards explicit update!
	VectorXf forces = VectorXf::Zero( initialGridVelocities.size() );
	g.calculateForces( d, forces );
	
	const VectorXf& gridVelocities = g.getVelocities();
	const VectorXf& gridMasses = g.masses();
	for( int i=0; i < gridMasses.size(); ++i )
	{
		if( gridMasses[i] > 0 )
		{
			std::cerr << initialGridVelocities.segment<3>( 3*i ).transpose() << " ==  ";
			std::cerr << ( gridVelocities.segment<3>( 3*i ) - ( timeStep / gridMasses[i] ) * forces.segment<3>( 3*i ) ).transpose() << "?" << std::endl;
		}
	}
	std::cerr << "dun" << std::endl;
}

void testGrid()
{
	testMassSplatting();
	testDeformationGradients();
	testForces();
	testImplicitUpdate();
}

}
