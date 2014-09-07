#include <Eigen/Geometry>

#include "tests/TestGrid.h"

#include "MpmSim/Grid.h"
#include "MpmSim/MassMatrix.h"
#include "MpmSim/CollisionPlane.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/SnowConstitutiveModel.h"

#include <iostream>
#include <fstream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

static void testProcessingPartitions()
{
	std::cerr << "testProcessingPartitions()" << std::endl;
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

	std::cerr << "testSplatting()" << std::endl;
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
	std::cerr << "testDeformationGradients()" << std::endl;
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
	std::cerr << "testForces()" << std::endl;
	// create a single particle:
	Sim::MaterialPointDataMap particleData;

	VectorData* p = new VectorData( 1, Eigen::Vector3f::Zero() );
	particleData["p"] = p;
	
	VectorData* v = new VectorData( 1, Eigen::Vector3f::Zero() );
	particleData["v"] = v;
	
	// initialize deformation gradient with an arbitrary rotation:
	Matrix3f rotato = Matrix3f::Identity();
	rotato = AngleAxisf(float( 0.25*M_PI ), Vector3f::UnitX())
	  * AngleAxisf(float( 0.5*M_PI ),  Vector3f::UnitY())
	  * AngleAxisf(float( 0.33*M_PI ), Vector3f::UnitZ());

	MatrixData* f = new MatrixData( 1, rotato );
	particleData["F"] = f;
	
	ScalarData* m = new ScalarData( 1, 1 );
	particleData["m"] = m;

	ScalarData* volume = new ScalarData( 1, 0.5f );
	particleData["volume"] = volume;
	
	const float gridSize = 1.0f;
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel constitutiveModel(
		1.4e3f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);

	constitutiveModel.createParticleData( particleData );
	constitutiveModel.setParticles( particleData );
	constitutiveModel.updateParticleData( particleData );

	Sim::IndexList inds;
	inds.push_back(0);
	Grid g( particleData, inds, gridSize, shapeFunction );

	// at equilibrium, so no forces on the grid nodes:
	Sim::ForceFieldSet forceFields;
	Eigen::VectorXf forces( g.velocities.size() );
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	float minF = fabs( forces.minCoeff() );
	float maxF = fabs( forces.maxCoeff() );
	assert( minF < 1.e-4 && maxF < 1.e-4 );
	
	
	// fiddle with F for one of the particles to get it out of equilibrium:
	std::vector<Eigen::Matrix3f>& F = dynamic_cast<MatrixData*>( particleData["F"] )->m_data;
	F[0] += 0.01f * Eigen::Matrix3f::Random();
	constitutiveModel.updateParticleData( particleData );
	
	// recompute forces on grid nodes:
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	
	Eigen::Vector3f fracDimPos = ( Eigen::Vector3f::Zero() - g.minCoord() ) / gridSize;
	const Eigen::Vector3i& n = g.n();
	int centerIdx = int(fracDimPos[0]) + n[0] *( int(fracDimPos[1]) + n[1] * int(fracDimPos[2]) );
	
	Eigen::Matrix3f Forig = F[0];
	
	Eigen::VectorXf forcesPerturbed( g.velocities.size() );
	Eigen::VectorXf forcesPerturbedNegative( g.velocities.size() );
	
	// now work out the forces on the grid nodes by finite differences!
	float e0 = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
	const float dx = 0.01f;
	for( int i=0; i <g.velocities.size(); ++i )
	{
		if( forces[i] == 0 )
		{
			continue;
		}
		
		F[0] = Forig;
		constitutiveModel.updateParticleData( particleData );
		VectorXf df(g.velocities.size());
		g.velocities[i] = dx;
		g.calculateForceDifferentials(
			df,
			g.velocities,
			constitutiveModel,
			forceFields.fields );

		// use the deformation gradient update to move one of the nodes
		// a small distance along one of the axes, so it squishes the
		// particle around a bit and work out the energy in this perturbed
		// configuration:
		g.velocities[i] = dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		constitutiveModel.updateParticleData( particleData );
		float ePlus = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
		g.calculateForces( forcesPerturbed, constitutiveModel, forceFields.fields );
		
		// now deform it the other way and work out the energy:
		g.velocities[i] = -dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		constitutiveModel.updateParticleData( particleData );
		float eMinus = constitutiveModel.energyDensity( 0 ) * volume->m_data[0];
		g.calculateForces( forcesPerturbedNegative, constitutiveModel, forceFields.fields );
		
		// x component of force on this node is -dE/dx, same for y and z. Make fd and analytic values
		// agree within a percent:
		assert( fabs( ( forces[i] + ( ePlus - eMinus ) / ( 2 * dx ) ) / forces[i] ) < 0.01f );
		
		// I wish I knew why this stuff converged so badly... looks like I can only get the force
		// differentials to agree with the finite difference approximation to within 2.5% or something.
		std::cerr << ( 0.5f * ( forcesPerturbed - forcesPerturbedNegative ) - df ).norm() / df.norm() << " " << df.norm() << std::endl;
		assert( ( 0.5f * ( forcesPerturbed - forcesPerturbedNegative ) - df ).norm() / df.norm() < 0.025f );

		// put things back in their place
		g.velocities[i] = 0;
	}
}

class ImplicitUpdateRecord : public LinearSolver::Debug
{
public:

	ImplicitUpdateRecord(
		const Grid& g,
		const Eigen::VectorXf& explicitMomenta,
		const std::vector<const CollisionObject*>& collisionObjects,
		const std::vector<char>& nodeCollided,
		const std::string& fileName
	) :
		outFile( fileName.c_str(), std::ofstream::binary )
	{
		g.collisionVelocities( vc, collisionObjects, nodeCollided );

		Eigen::VectorXf explicitVelocities( explicitMomenta.size() );
		for( int i=0; i < g.masses.size(); ++i )
		{
			if( g.masses[i] != 0 )
			{
				explicitVelocities.segment<3>( 3 * i ) = explicitMomenta.segment<3>( 3 * i ) / g.masses[i];
			}
			else
			{
				explicitVelocities.segment<3>( 3 * i ).setZero();
			}
		}
		explicitVelocities += vc;

		float gridH = g.gridSize();
		outFile.write( (const char*)&gridH, sizeof( float ) );

		outFile.write( (const char*)&(g.minCoord()[0]), sizeof( float ) );
		outFile.write( (const char*)&(g.minCoord()[1]), sizeof( float ) );
		outFile.write( (const char*)&(g.minCoord()[2]), sizeof( float ) );

		outFile.write( (const char*)&(g.n()[0]), sizeof( int ) );
		outFile.write( (const char*)&(g.n()[1]), sizeof( int ) );
		outFile.write( (const char*)&(g.n()[2]), sizeof( int ) );

		outFile.write( (const char*)&(g.masses[0]), g.n()[0] * g.n()[1] * g.n()[2] * sizeof( float ) );
		
		outFile.write( (const char*)&(explicitVelocities[0]), 3 * g.n()[0] * g.n()[1] * g.n()[2] * sizeof( float ) );
	}
	
	virtual void operator()( Eigen::VectorXf& x )
	{
		Eigen::VectorXf iterate = x + vc;
		outFile.write( (const char*)&(iterate[0]), iterate.size() * sizeof( float ) );
	}
	
	Eigen::VectorXf vc;
	std::ofstream outFile;
	std::vector< Eigen::VectorXf > record;

};


void testImplicitUpdate()
{
	std::cerr << "testImplicitUpdate()" << std::endl;

	// create some particules:
	Sim::MaterialPointDataMap particleData;

	VectorData* p = new VectorData;
	particleData["p"] = p;
	
	VectorData* v = new VectorData;
	particleData["v"] = v;

	MatrixData* fData = new MatrixData;
	particleData["F"] = fData;
	
	ScalarData* m = new ScalarData;
	particleData["m"] = m;

	ScalarData* volume = new ScalarData;
	particleData["volume"] = volume;

	std::vector<Vector3f>& velocities = v->m_data;
	std::vector<Vector3f>& positions = p->m_data;
	std::vector<Matrix3f>& F = fData->m_data;
	std::vector<float>& masses = m->m_data;
	std::vector<float>& volumes = volume->m_data;
	
	Matrix3f distortion = Matrix3f::Random() * 0.001f;
	
	const float gridSize = 0.5f;
	Sim::IndexList inds;
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				inds.push_back( (int)positions.size() );
				positions.push_back(
					Vector3f(
						float( i -0.5f ) * gridSize,
						float( j - 0.5f ) * gridSize,
						float( k - 0.5f ) * gridSize
					)
				);
				
				masses.push_back( 1.0f );
				volumes.push_back( 1.0f );
				
				velocities.push_back(
					Vector3f(
						0.1f + 0.01f*sin( 2 * positions.back()[0] ),
						0.1f + 0.01f*sin( 2 * positions.back()[1] ),
						0.1f + 0.01f*sin( 2 * positions.back()[2] )
					)
				);
				
				F.push_back(
					Matrix3f::Identity() + distortion * cos( 2 * positions.back()[0] )
				);
				
			}
		}
	}
	
	// create particle data:
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel snowModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);
	float voxelSize( 2 * shapeFunction.supportRadius() * gridSize );
	Sim::voxelSort(
		inds.begin(),
		inds.end(),
		voxelSize,
		positions );
	
	snowModel.createParticleData( particleData );
	snowModel.setParticles( particleData );
	Grid g( particleData, inds, gridSize, shapeFunction );
	g.computeParticleVolumes();
	
	float timeStep = 0.002f;
	Sim::CollisionObjectSet collisionObjects;
	CollisionPlane* plane1 = new CollisionPlane( Eigen::Vector4f( -1,-1,-1,1 ) );
	plane1->setV( Eigen::Vector3f( -.1f, -.2f, -.3f ) );
	collisionObjects.objects.push_back( plane1 );

	CollisionPlane* plane2 = new CollisionPlane( Eigen::Vector4f( -1,0,0,1.5 ) );
	plane2->setV( Eigen::Vector3f( -.2f, -.1f, -.4f ) );
	collisionObjects.objects.push_back( plane2 );

	std::vector<const ForceField*> fields;
	ConjugateResiduals solver( 400, 0 );
	
	VectorXf explicitMomenta;
	std::vector<char> nodeCollided;
	g.calculateExplicitMomenta(
		explicitMomenta,
		nodeCollided,
		timeStep,
		snowModel,
		collisionObjects,
		fields
	);

	{
		ImplicitUpdateRecord d( g, explicitMomenta, collisionObjects.objects, nodeCollided, "debug.dat" );
		g.updateGridVelocities(
			timeStep, 
			snowModel,
			collisionObjects,
			fields,
			solver,
			&d
		);
	}

	//system( "C:\\Users\\david\\Documents\\GitHub\\mpmsnow\\Debug\\viewer.exe" );

	// check nothing's moving in or out of the collision objects:
	VectorXf vc( g.velocities.size() );
	g.collisionVelocities( vc, collisionObjects.objects, nodeCollided );
	for( int i=0; i < g.n()[0]; ++i )
	{
		for( int j=0; j < g.n()[1]; ++j )
		{
			for( int k=0; k < g.n()[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				Vector3f velocity = g.velocities.segment<3>( 3 * idx );
				Vector3f x(
					g.gridSize() * i + g.minCoord()[0],
					g.gridSize() * j + g.minCoord()[1],
					g.gridSize() * k + g.minCoord()[2]
				);

				int objIdx = nodeCollided[idx];
				if( objIdx >= 0 )
				{
					const CollisionObject* obj = collisionObjects.objects[ objIdx ];
					
					// express velocity and momentum relative to this object:
					Vector3f vObj;
					obj->velocity( x, vObj );
					assert( vc.segment<3>( 3 * idx ) == vObj );
					velocity -= vObj;
					
					// find object normal:
					Vector3f n;
					obj->grad( x, n );
					n.normalize();
					float nDotV = fabs( n.dot( velocity ) );
					assert( nDotV < 1.e-5 );
				}
				else if( objIdx == -2 )
				{
					assert( velocity.norm() == 0 );
					assert( vc.segment<3>( 3 * idx ).norm() == 0 );
				}
			}
		}
	}
	
	// so.... did it work?
	
	// ok - so it should have solved this equation:
	// P * ( M * ( g.velocities - vc ) - DF * g.velocities * dt * dt ) = explicitMomenta
	// (here, P means we project out the collided velocity components)
	
	for( int i=0; i < g.n()[0]; ++i )
	{
		for( int j=0; j < g.n()[1]; ++j )
		{
			for( int k=0; k < g.n()[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				
				// apply P to velocities:
				if( nodeCollided[idx] >= 0 )
				{
					const CollisionObject* obj = collisionObjects.objects[ nodeCollided[idx] ];
					
					// find object normal:
					Vector3f x(
						g.gridSize() * i + g.minCoord()[0],
						g.gridSize() * j + g.minCoord()[1],
						g.gridSize() * k + g.minCoord()[2]
					);
					
					Vector3f vObj;
					obj->velocity( x, vObj );
					assert( vc.segment<3>( 3 * idx ) == vObj );	
				}
				else
				{
					assert( vc.segment<3>( 3 * idx ).norm() == 0 );
				}
			}
		}
	}
	
	// work out DF * g.velocities:
	VectorXf df( g.velocities.size() );
	g.calculateForceDifferentials( df, g.velocities, snowModel, fields );

	//  M * ( g.velocities - vc ) - DF * g.velocities * dt * dt 
	VectorXf lhs( g.velocities.size() );
	for( int idx=0; idx < g.masses.size(); ++idx )
	{
		lhs.segment<3>( 3 * idx ) =
			g.masses[idx] * ( g.velocities.segment<3>( 3 * idx ) - vc.segment<3>( 3 * idx ) ) - df.segment<3>( 3 * idx ) * timeStep * timeStep;
	}
	
	// run the projection again:
	for( int i=0; i < g.n()[0]; ++i )
	{
		for( int j=0; j < g.n()[1]; ++j )
		{
			for( int k=0; k < g.n()[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				
				// apply P:
				if( nodeCollided[idx] >= 0 )
				{
					const CollisionObject* obj = collisionObjects.objects[ nodeCollided[idx] ];
					
					// find object normal:
					Vector3f x(
						g.gridSize() * i + g.minCoord()[0],
						g.gridSize() * j + g.minCoord()[1],
						g.gridSize() * k + g.minCoord()[2]
					);
					
					Vector3f momentum = lhs.segment<3>( 3 * idx );
					Vector3f n;
					obj->grad( x, n );
					n.normalize();
					float nDotP = n.dot( momentum );
					
					// project out component perpendicular to the object
					lhs.segment<3>( 3 * idx ) = momentum - nDotP * n;
				}
				else if( nodeCollided[idx] == -2 )
				{
					lhs.segment<3>( 3 * idx ).setZero();
				}
			}
		}
	}
	
	float maxNorm(0);
	for( int idx=0; idx < g.masses.size(); ++idx )
	{
		float norm = ( lhs.segment<3>( 3 * idx ) - explicitMomenta.segment<3>( 3 * idx ) ).norm();
		if( norm > maxNorm )
		{
			maxNorm = norm;
		}
	}
	assert( maxNorm < 1.e-6 );
}

void testMovingGrid()
{
	std::cerr << "testMovingGrid" << std::endl;
	
	// create some particules:
	Sim::MaterialPointDataMap particleData;

	VectorData* p = new VectorData;
	particleData["p"] = p;
	
	VectorData* v = new VectorData;
	particleData["v"] = v;

	MatrixData* fData = new MatrixData;
	particleData["F"] = fData;
	
	ScalarData* m = new ScalarData;
	particleData["m"] = m;

	ScalarData* volume = new ScalarData;
	particleData["volume"] = volume;


	std::vector<Vector3f>& velocities = v->m_data;
	std::vector<Vector3f>& positions = p->m_data;
	std::vector<Matrix3f>& F = fData->m_data;
	std::vector<float>& masses = m->m_data;
	std::vector<float>& volumes = volume->m_data;
	
	const float gridSize = 0.5f;
	Sim::IndexList inds;
	for( int i=0; i < 4; ++i )
	{
		for( int j=0; j < 4; ++j )
		{
			for( int k=0; k < 4; ++k )
			{
				inds.push_back( (int)positions.size() );
				positions.push_back(
					Vector3f(
						float( i -0.5f ) * gridSize,
						float( j - 0.5f ) * gridSize,
						float( k - 0.5f ) * gridSize
					)
				);
				
				masses.push_back( 1.0f );
				volumes.push_back( 1.0f );
				velocities.push_back( Vector3f( 0, -1, 0 ) );
				F.push_back( Matrix3f::Identity() );
				
			}
		}
	}
	
	// create particle data:
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel snowModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);
	float voxelSize( 2 * shapeFunction.supportRadius() * gridSize );
	Sim::voxelSort(
		inds.begin(),
		inds.end(),
		voxelSize,
		positions );
	
	snowModel.createParticleData( particleData );
	snowModel.setParticles( particleData );
	Grid g( particleData, inds, gridSize, shapeFunction, Vector3f( 0, -1, 0 ) );
	g.computeParticleVolumes();
	
	float timeStep = 0.002f;
	Sim::CollisionObjectSet collisionObjects;
	CollisionPlane* plane1 = new CollisionPlane( Eigen::Vector4f( 0,1,0,1 ) );
	collisionObjects.objects.push_back( plane1 );
	
	std::vector<const ForceField*> fields;
	ConjugateResiduals solver( 400, 0 );
	
	VectorXf explicitMomenta;
	std::vector<char> nodeCollided;
	g.calculateExplicitMomenta(
		explicitMomenta,
		nodeCollided,
		timeStep,
		snowModel,
		collisionObjects,
		fields
	);
	
	// ok, at the moment I'm just doing something and eyeballing the result - how about
	// I do a proper unit test for this some time?
	{
		ImplicitUpdateRecord d( g, explicitMomenta, collisionObjects.objects, nodeCollided, "debug.dat" );
		g.updateGridVelocities(
			timeStep, 
			snowModel,
			collisionObjects,
			fields,
			solver,
			&d
		);
	}
	system( "C:\\Users\\david\\Documents\\GitHub\\mpmsnow\\Debug\\viewer.exe" );


}


void testGrid()
{
	std::cerr << "testGrid()" << std::endl;
	//testProcessingPartitions();
	//testSplatting();
	//testDeformationGradients();
	//testForces();
	//testImplicitUpdate();
	testMovingGrid();
}

}
