#include "tests/TestGrid.h"

#include "MpmSim/Grid.h"
#include "MpmSim/ConstitutiveModel.h"
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
	virtual void forceDifferentialDensity( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const
	{
	}

};

void testGrid()
{
	std::vector<Vector3f> positions;
	std::vector<float> masses;
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

}
