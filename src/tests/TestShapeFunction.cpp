
#include "MpmSim/ParticleData.h"
#include "MpmSim/Grid.h"

#include "tests/TestShapeFunction.h"

using namespace MpmSim; 

static float evaluateShapeFunction( const ShapeFunction& shapeFunction, const Eigen::Vector3f& origin, float h, const Eigen::Vector3f& evalPos )
{
	Eigen::Vector3f relativePos = ( evalPos - origin ) / h;
	return shapeFunction.w( relativePos[0] ) * shapeFunction.w( relativePos[1] ) * shapeFunction.w( relativePos[2] );
}

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

void testShapeFunction( const ShapeFunction& shapeFunction )
{
	assert( shapeFunction.supportRadius() == 2 );

	assert( shapeFunction.w( (float)shapeFunction.supportRadius() ) == 0. );
	assert( shapeFunction.w( -(float)shapeFunction.supportRadius() ) == 0.0f );

	assert( shapeFunction.dw( (float)shapeFunction.supportRadius() ) == 0.0f );
	assert( shapeFunction.dw( -(float)shapeFunction.supportRadius() ) == 0.0f );
	
	const float dx = 0.001f;
	for( int i=0; i <= 20; ++i )
	{
		float x = 0.1f * i;

		// test symmetry:
		assert( shapeFunction.w( x ) == shapeFunction.w( -x ) );
		assert( shapeFunction.dw( x ) == -shapeFunction.dw( -x ) );
		
		// test mass conservation:
		float sum = shapeFunction.w( x-4 ) + shapeFunction.w( x-3 ) + shapeFunction.w( x-2 ) + shapeFunction.w( x-1 ) + shapeFunction.w( x ) + shapeFunction.w( x+1 ) + shapeFunction.w( x+2 );
		assert( fabs( sum - 1 ) < 1.e-6f );

		// test derivative:
		float dwFiniteDiff = ( shapeFunction.w( x + dx ) - shapeFunction.w( x - dx ) ) / ( 2 * dx );
		assert( fabs( dwFiniteDiff - shapeFunction.dw( x ) ) < 1.e-4f );
	}
	
	// test iterator for splatting to grid:
	Eigen::Vector3f gridOrigin( -1.7638f, -2.85363f, 2.25722f );
	Eigen::Vector3f particlePos( 0.15321f, 2.29587345f, -3.897576f );
	Eigen::Vector3i gridPos;
	const float gridH = 0.314857f;
	
	std::vector<Eigen::Vector3f> positions;
	positions.push_back( gridOrigin );

	std::vector<float> masses;
	masses.push_back( 1 );

	ParticleData d( positions, masses, gridH );
	TestConstitutiveModel model;
	std::vector<CollisionObject*> collisionObjects;
	Grid g(d, collisionObjects, 0.01f, shapeFunction, model );
	Grid::PointToGridIterator it( g );
	it.initialize( particlePos, true );
	int pointsVisited = 0;
	float totalWeight = 0;

	Eigen::Vector3f minPos(1000000,1000000,1000000);
	Eigen::Vector3f maxPos(-1000000,-1000000,-1000000);

	Eigen::Vector3i minGridPos(1000000,1000000,1000000);
	Eigen::Vector3i maxGridPos(-1000000,-1000000,-1000000);
	
	do
	{
		it.gridPos( gridPos );
		if( gridPos[0] < minGridPos[0] ) minGridPos[0] = gridPos[0];
		if( gridPos[1] < minGridPos[1] ) minGridPos[1] = gridPos[1];
		if( gridPos[2] < minGridPos[2] ) minGridPos[2] = gridPos[2];
		
		if( gridPos[0] > maxGridPos[0] ) maxGridPos[0] = gridPos[0];
		if( gridPos[1] > maxGridPos[1] ) maxGridPos[1] = gridPos[1];
		if( gridPos[2] > maxGridPos[2] ) maxGridPos[2] = gridPos[2];
		
		Eigen::Vector3f worldPos;
		worldPos[0] = gridPos[0] * gridH + gridOrigin[0];
		worldPos[1] = gridPos[1] * gridH + gridOrigin[1];
		worldPos[2] = gridPos[2] * gridH + gridOrigin[2];

		if( worldPos[0] < minPos[0] ) minPos[0] = worldPos[0];
		if( worldPos[1] < minPos[1] ) minPos[1] = worldPos[1];
		if( worldPos[2] < minPos[2] ) minPos[2] = worldPos[2];
		
		if( worldPos[0] > maxPos[0] ) maxPos[0] = worldPos[0];
		if( worldPos[1] > maxPos[1] ) maxPos[1] = worldPos[1];
		if( worldPos[2] > maxPos[2] ) maxPos[2] = worldPos[2];
		
		float w = evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos );
		assert( fabs( w - it.w() ) < 1.e-4f );
		
		Eigen::Vector3f weightGrad;
		it.dw( weightGrad );

		float dwx = -(
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( dx, 0, 0 ) ) -
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( -dx, 0, 0 ) )
		) / ( 2 * dx );

		float dwy = -(
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( 0, dx, 0 ) ) -
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( 0, -dx, 0 ) )
		) / ( 2 * dx );

		float dwz = -(
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( 0, 0, dx ) ) -
			evaluateShapeFunction( shapeFunction, particlePos, gridH, worldPos + Eigen::Vector3f( 0, 0, -dx ) )
		) / ( 2 * dx );
		
		assert( fabs( weightGrad[0] - dwx ) < 1.e-3 );
		assert( fabs( weightGrad[1] - dwy ) < 1.e-3 );
		assert( fabs( weightGrad[2] - dwz ) < 1.e-3 );

		totalWeight += it.w();
		++pointsVisited;
	} while( it.next() );
	
	// have we iterated over the right number of points?
	assert( pointsVisited == 8 * shapeFunction.supportRadius() * shapeFunction.supportRadius() * shapeFunction.supportRadius() );

	// is mass conserved?
	assert( fabs( totalWeight - 1 ) < 1.0e-6f );

	// correct bounds?
	Eigen::Vector3f gridRelativeParticle = ( particlePos - gridOrigin )/gridH;
	Eigen::Vector3f expectedMin(
		gridOrigin[0] + ( floor( gridRelativeParticle[0] ) - 1 ) * gridH,
		gridOrigin[1] + ( floor( gridRelativeParticle[1] ) - 1 ) * gridH,
		gridOrigin[2] + ( floor( gridRelativeParticle[2] ) - 1 ) * gridH
	);
	Eigen::Vector3f expectedMax = expectedMin + ( 2.0f * shapeFunction.supportRadius() - 1.0f ) * Eigen::Vector3f( gridH, gridH, gridH );
	
	for( int i=0; i < 3; ++i )
	{
		assert( fabs( minPos[i] - expectedMin[i] ) < 1.e-4 );
		assert( fabs( maxPos[i] - expectedMax[i] ) < 1.e-4 );
	}

	Eigen::Vector3i expectedGridMin(
		(int)floor( gridRelativeParticle[0] ) - 1,
		(int)floor( gridRelativeParticle[1] ) - 1,
		(int)floor( gridRelativeParticle[2] ) - 1
	);
	Eigen::Vector3i expectedGridMax = expectedGridMin + ( 2 * shapeFunction.supportRadius() - 1 ) * Eigen::Vector3i( 1, 1, 1 );
	
	assert( minGridPos == expectedGridMin );
	assert( maxGridPos == expectedGridMax );
	
}

}
