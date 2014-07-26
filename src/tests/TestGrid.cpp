#include "tests/TestGrid.h"

#include "MpmSim/Grid.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/SnowConstitutiveModel.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

static void testProcessingPartitions()
{
	VectorData* positionData = new VectorData;
	VectorData* velocityData = new VectorData;
	ScalarData* massData = new ScalarData;

	// make some initial particles:
	std::vector<Vector3f>& positions = positionData->m_data;
	std::vector<Vector3f>& velocities = velocityData->m_data;
	std::vector<float>& masses = massData->m_data;
	
	Sim::IndexList particleInds;

	float gridSize(0.01f);
	const int n = 16;
	for( int i=0; i < n; ++i )
	{
		for( int j=0; j < n; ++j )
		{
			for( int k=0; k < n; ++k )
			{
				velocities.push_back( Vector3f::Zero() );
				positions.push_back( Vector3f( 0.5f * gridSize * ( i + 0.5f ), 0.5f * gridSize * ( j + 0.5f ), 0.5f * gridSize * ( k + 0.5f ) ) );
				masses.push_back( 1.0f );
				particleInds.push_back( (int)particleInds.size() );
			}
		}
	}
	
	CubicBsplineShapeFunction shapeFunction;
	float voxelSize( 2 * shapeFunction.supportRadius() * gridSize );
	Sim::voxelSort(
		particleInds.begin(),
		particleInds.end(),
		voxelSize,
		positions );
	
	Sim::MaterialPointDataMap d;
	d["p"] = positionData;
	d["v"] = velocityData;
	d["m"] = massData;
	
	// check the indices have the right spatial properties:
	Vector3i currentVoxel;
	int numVoxels(0);
	for( size_t p=0; p<particleInds.size(); ++p )
	{
		Vector3f x = positions[ particleInds[p] ] / ( voxelSize );
		Vector3i voxel( (int)floor( x[0] ), (int)floor( x[1] ), (int)floor( x[2] ) );
		if( p == 0 || currentVoxel != voxel )
		{
			++numVoxels;
			// ok, we're gonna test the lexicographical ordering, meaning we demand
			// that the first non zero component of the voxel difference is positive:
			Vector3i voxelDiff = ( voxel - currentVoxel );
			for( int i=0; i < 3; ++i )
			{
				if( voxelDiff[i] != 0 )
				{
					assert( voxelDiff[i] > 0 );
					break;
				}
			}
			currentVoxel = voxel;
		}
	}
	
	Grid g( d, particleInds, gridSize, shapeFunction );
	
	// now check the processing partitioning
	int numPartitionedVoxels(0);
	int numPartitionedParticles(0);
	for( int i=0; i < 2; ++i )
	{
		for( int j=0; j < 2; ++j )
		{
			for( int k=0; k < 2; ++k )
			{
				const Grid::PartitionList& partition = g.partition( i, j, k );
				numPartitionedVoxels += (int)partition.size();
				for( size_t v=0; v < partition.size(); ++v )
				{
					Sim::ConstIndexIterator begin = partition[v].first;
					Sim::ConstIndexIterator end = partition[v].second;
					numPartitionedParticles += (int)(end - begin);
					Vector3i currentVoxel;
					for( Sim::ConstIndexIterator it = begin; it != end; ++it )
					{
						// all particles must be in the same voxel...
						// and the voxel must be in the right partition
						Vector3f x = positions[ *it ] / ( 4 * gridSize );
						Vector3i voxel( (int)floor( x[0] ), (int)floor( x[1] ), (int)floor( x[2] ) );
						if( it == begin )
						{
							currentVoxel = voxel;
						}
						assert( voxel == currentVoxel );
						assert( ( voxel[0] & 1 ) == i );
						assert( ( voxel[1] & 1 ) == j );
						assert( ( voxel[2] & 1 ) == k );
					}
				}
			}
		}
	}
	assert( numPartitionedVoxels == numVoxels );
	assert( numPartitionedParticles == positions.size() );
	
	for( Sim::MaterialPointDataMap::iterator it = d.begin(); it != d.end(); ++it )
	{
		delete it->second;
	}
}

static void testSplatting()
{

	VectorData* positionData = new VectorData;
	VectorData* velocityData = new VectorData;
	ScalarData* massData = new ScalarData;

	// make some initial particles:
	std::vector<Vector3f>& positions = positionData->m_data;
	std::vector<Vector3f>& velocities = velocityData->m_data;
	std::vector<float>& masses = massData->m_data;
	
	Sim::IndexList particleInds;
	
	positions.push_back( Vector3f::Zero() );
	velocities.push_back( Vector3f(1,0,0) + Vector3f::Random() );
	masses.push_back( 1.0f );
	particleInds.push_back( (int)particleInds.size() );
	for( int i=0; i < 500; ++i )
	{
		float xr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float yr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float zr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
		positions.push_back( Vector3f( xr, yr, zr ) );
		velocities.push_back( Vector3f(1,0,0) + Vector3f::Random() );
		masses.push_back( 1.0f );
		particleInds.push_back( (int)particleInds.size() );
	}
	
	CubicBsplineShapeFunction shapeFunction;
	const float gridSize = 0.01f;
	
	float voxelSize( 2 * shapeFunction.supportRadius() * gridSize );
	Sim::voxelSort(
		particleInds.begin(),
		particleInds.end(),
		voxelSize,
		positions );
	
	Sim::MaterialPointDataMap d;
	d["p"] = positionData;
	d["v"] = velocityData;
	d["m"] = massData;
	
	// create a grid, which will immediately splat the particle masses onto itself,
	// multi threaded tbb stylee:
	Grid g( d, particleInds, gridSize, shapeFunction );
	
	// now manually splat the masses, serial stylee:
	Eigen::VectorXf massVector( g.masses.size() );
	massVector.setZero();
	
	Grid::ShapeFunctionIterator& shIt = g.shapeFunctionIterator();
	Vector3i particleCell;
	for( size_t p = 0; p < positions.size(); ++p )
	{
		shIt.initialize( positions[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
			massVector[ idx ] += masses[p] * shIt.w();
		} while( shIt.next() );
	}	
	
	massVector -= g.masses;
	assert( fabs( massVector.maxCoeff() ) < 1.e-6 );
	assert( fabs( massVector.minCoeff() ) < 1.e-6 );
	
	// check conservation of mass:
	float totalGridMass = g.masses.sum();
	float totalParticleMass = 0;
	for( size_t p = 0; p < masses.size(); ++p )
	{
		totalParticleMass += masses[p];
	}
	
	// ideally these should be equal, but they're gonna differ for numerical reasons:
	assert( fabs( totalParticleMass - totalGridMass ) / totalParticleMass < 1.e-4 );
	
	// check conservation of momentum:
	Eigen::Vector3f gridMomentum = Eigen::Vector3f::Zero();
	for( int i=0; i < g.masses.size(); ++i )
	{
		gridMomentum += g.masses[i] * g.velocities.segment<3>( 3 * i );
	}
	Eigen::Vector3f particleMomentum = Eigen::Vector3f::Zero();
	for( size_t p = 0; p < masses.size(); ++p )
	{
		particleMomentum += masses[p] * velocities[p];
	}
	// again: should be equal but for numerical shiz:
	assert( ( particleMomentum - gridMomentum ).norm() / particleMomentum.norm() < 1.e-4 );

	for( Sim::MaterialPointDataMap::iterator it = d.begin(); it != d.end(); ++it )
	{
		delete it->second;
	}
}


static void testDeformationGradients()
{
	// unit cube full of particles centered at the origin:
	Sim::MaterialPointDataMap particleData;

	VectorData* p = new VectorData;
	particleData["p"] = p;
	
	VectorData* v = new VectorData;
	particleData["v"] = v;

	MatrixData* f = new MatrixData;
	particleData["F"] = f;
	
	ScalarData* m = new ScalarData;
	particleData["m"] = m;

	ScalarData* volume = new ScalarData;
	particleData["volume"] = volume;
	
	std::vector<Matrix3f>& F = f->m_data;
	std::vector<Vector3f>& velocities = v->m_data;
	std::vector<Vector3f>& positions = p->m_data;
	std::vector<float>& masses = m->m_data;
	std::vector<float>& volumes = volume->m_data;
	Sim::IndexList inds;
	
	// set the velocities up so they expand in the x/z directions and shrink in y:

	// add a particle at the origin, which is the one we be testin'
	inds.push_back( (int)inds.size() );
	positions.push_back( Vector3f::Zero() );
	masses.push_back( 1.0f );
	velocities.push_back( Eigen::Vector3f::Zero() );
	F.push_back( Eigen::Matrix3f::Identity() );
	volumes.push_back( 1.0f );

	for( int i=0; i < 10; ++i )
	{
		for( int j=0; j < 10; ++j )
		{
			for( int k=0; k < 10; ++k )
			{
				inds.push_back( (int)inds.size() );
				positions.push_back( Vector3f( float(i) / 9 - 0.5f, float(j) / 9 - 0.5f, float(k) / 9 - 0.5f ) );
				masses.push_back( 1.0f );
				velocities.push_back( Eigen::Vector3f( positions.back()[0], -positions.back()[1], 0.5f * positions.back()[2] ) );
				F.push_back( Eigen::Matrix3f::Identity() );
				volumes.push_back( 1.0f );
			}
		}
	}
	
	CubicBsplineShapeFunction shapeFunction;
	const float gridSize = 0.2f;
	const float timeStep = 0.01f;
	for( int i=1; i < 8; ++i )
	{
		Grid g( particleData, inds, gridSize, shapeFunction );
		g.updateDeformationGradients( timeStep );
		for( size_t p=0; p < velocities.size(); ++p )
		{
			positions[p] += velocities[p] * timeStep;
		}

		assert( fabs( F[0]( 0,0 ) - ( 1 + 0.01f * i ) ) < 0.001f );
		assert( fabs( F[0]( 1,1 ) - ( 1 - 0.01f * i ) ) < 0.001f );
		assert( fabs( F[0]( 2,2 ) - ( 1 + 0.005f * i ) ) < 0.001f );
	}

	// NB: significant errors get introduced to the deformation gradient round the edges of
	// the object it seems, which is probably responsible for some visual artefacts... Maybe
	// it's worth fixing these some time?

	std::cerr << std::endl;

}


static void testForces()
{
	// create a single particle:
	Sim::MaterialPointDataMap particleData;

	VectorData* p = new VectorData( 1, Eigen::Vector3f::Zero() );
	particleData["p"] = p;
	
	VectorData* v = new VectorData( 1, Eigen::Vector3f::Zero() );
	particleData["v"] = v;

	MatrixData* f = new MatrixData( 1, Eigen::Matrix3f::Identity() );
	particleData["F"] = f;
	
	ScalarData* m = new ScalarData( 1, 1 );
	particleData["m"] = m;

	ScalarData* volume = new ScalarData( 1, 1.0f );
	particleData["volume"] = volume;
	
	const float gridSize = 1.0f;
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel constitutiveModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);

	constitutiveModel.createParticleData( particleData );
	constitutiveModel.setParticles( particleData );
	
	Sim::IndexList inds;
	inds.push_back(0);
	Grid g( particleData, inds, gridSize, shapeFunction );

	// at equilibrium, so no forces on the grid nodes:
	Sim::ForceFieldSet forceFields;
	Eigen::VectorXf forces( g.velocities.size() );
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	assert( forces.minCoeff() == forces.maxCoeff() );
	assert( forces.maxCoeff() == 0 );
	
	
	// fiddle with F for one of the particles to get it out of equilibrium:
	std::vector<Eigen::Matrix3f>& F = dynamic_cast<MatrixData*>( particleData["F"] )->m_data;
	F[0] += 0.1f * Eigen::Matrix3f::Random();
	
	// recompute forces on grid nodes:
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	
	Eigen::Vector3f fracDimPos = ( Eigen::Vector3f::Zero() - g.minCoord() ) / gridSize;
	const Eigen::Vector3i& n = g.n();
	int centerIdx = int(fracDimPos[0]) + n[0] *( int(fracDimPos[1]) + n[1] * int(fracDimPos[2]) );
	
	Eigen::Matrix3f Forig = F[0];
	
	// now work out the forces on the grid nodes by finite differences!

	// initial energy of the particle:
	float e0 = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
	const float dx = 0.1f;
	for( int i=0; i <g.velocities.size(); ++i )
	{
		if( forces[i] == 0 )
		{
			continue;
		}

		// use the deformation gradient update move one of the nodes
		// a small distance along one of the axes, so it squishes the
		// particle around a bit:
		g.velocities[i] = dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		
		// work out the energy in this perturbed configuration:
		float ePlus = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
		
		// now deform it the other way and work out the energy:
		g.velocities[i] = -dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		
		float eMinus = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
		
		// x component of force on this node is -dE/dx, same for y and z. Make fd and analytic values
		// agree within a percent:
		assert( fabs( ( forces[i] + ( ePlus - eMinus ) / ( 2 * dx ) ) / forces[i] ) < 0.01f );
		
		// put things back in their place
		g.velocities[i] = 0;
	}
	
}

/*
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
	std::vector<CollisionObject*> collisionObjects;
	Grid g( d, collisionObjects, timeStep, shapeFunction, snowModel );

	g.computeDensities();
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
	g.updateGridVelocities( ConjugateResiduals( 500, 1.e-7f ) );

	// transfer the grid velocities back onto the particles:
	g.updateParticleVelocities();
	
	// update particle deformation gradients:
	g.updateDeformationGradients();
	
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
*/

void testGrid()
{
	testProcessingPartitions();
	testSplatting();
	testDeformationGradients();
	testForces();
	//testImplicitUpdate();
}

}
