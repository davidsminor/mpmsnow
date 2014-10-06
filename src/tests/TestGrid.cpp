#include <Eigen/Geometry>

#include "tests/TestGrid.h"

#include "MpmSim/Grid.h"
#include "MpmSim/CollisionPlane.h"
#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/SquareMagnitudeTermination.h"

#include <iostream>
#include <fstream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

void TestGrid::testProcessingPartitions()
{
	std::cerr << "testProcessingPartitions()" << std::endl;
	MaterialPointData d;
	
	// make some initial particles:
	std::vector<Vector3f>& positions = d.variable<Vector3f>("p");
	std::vector<Vector3f>& velocities = d.variable<Vector3f>("v");
	std::vector<float>& masses = d.variable<float>("m");
	
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
	Grid::voxelSort(
		particleInds.begin(),
		particleInds.end(),
		voxelSize,
		positions );
	
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
				const Grid::ParticlesInVoxelList& partition = g.m_processingPartitions[i][j][k];
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
}

void TestGrid::testSplatting()
{

	std::cerr << "testSplatting()" << std::endl;
	MaterialPointData d;
	
	// make some initial particles:
	std::vector<Vector3f>& positions = d.variable<Vector3f>("p");
	std::vector<Vector3f>& velocities = d.variable<Vector3f>("v");
	std::vector<float>& masses = d.variable<float>("m");
	
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
	Grid::voxelSort(
		particleInds.begin(),
		particleInds.end(),
		voxelSize,
		positions );
	
	// create a grid, which will immediately splat the particle masses onto itself,
	// multi threaded tbb stylee:
	Grid g( d, particleInds, gridSize, shapeFunction );
	
	// now manually splat the masses, serial stylee:
	Eigen::VectorXf massVector( g.m_masses.size() );
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
	
	massVector -= g.m_masses;
	assert( fabs( massVector.maxCoeff() ) < 1.e-6 );
	assert( fabs( massVector.minCoeff() ) < 1.e-6 );
	
	// check conservation of mass:
	float totalGridMass = g.m_masses.sum();
	float totalParticleMass = 0;
	for( size_t p = 0; p < masses.size(); ++p )
	{
		totalParticleMass += masses[p];
	}
	
	// ideally these should be equal, but they're gonna differ for numerical reasons:
	assert( fabs( totalParticleMass - totalGridMass ) / totalParticleMass < 1.e-4 );
	
	// check conservation of momentum:
	Eigen::Vector3f gridMomentum = Eigen::Vector3f::Zero();
	for( int i=0; i < g.m_masses.size(); ++i )
	{
		gridMomentum += g.m_masses[i] * g.m_velocities.segment<3>( 3 * i );
	}
	Eigen::Vector3f particleMomentum = Eigen::Vector3f::Zero();
	for( size_t p = 0; p < masses.size(); ++p )
	{
		particleMomentum += masses[p] * velocities[p];
	}
	// again: should be equal but for numerical shiz:
	assert( ( particleMomentum - gridMomentum ).norm() / particleMomentum.norm() < 1.e-4 );

}


void TestGrid::testDeformationGradients()
{
	std::cerr << "testDeformationGradients()" << std::endl;
	// unit cube full of particles centered at the origin:
	MaterialPointData particleData;
	
	std::vector<Matrix3f>& F = particleData.variable<Matrix3f>( "F" );
	std::vector<Vector3f>& velocities = particleData.variable<Vector3f>( "v" );
	std::vector<Vector3f>& positions = particleData.variable<Vector3f>( "p" );
	std::vector<float>& masses = particleData.variable<float>( "m" );
	std::vector<float>& volumes = particleData.variable<float>( "volume" );
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


void TestGrid::testForces()
{
	std::cerr << "testForces()" << std::endl;
	// create a single particle:
	MaterialPointData particleData;

	particleData.variable<Vector3f>("p").resize( 1, Eigen::Vector3f::Zero() );
	particleData.variable<Vector3f>("v").resize( 1, Eigen::Vector3f::Zero() );
	
	// initialize deformation gradient with an arbitrary rotation:
	Matrix3f rotato = Matrix3f::Identity();
	rotato = AngleAxisf(float( 0.25*M_PI ), Vector3f::UnitX())
	  * AngleAxisf(float( 0.5*M_PI ),  Vector3f::UnitY())
	  * AngleAxisf(float( 0.33*M_PI ), Vector3f::UnitZ());

	particleData.variable<Matrix3f>("F").resize( 1, rotato );

	particleData.variable<float>("m").resize( 1, 1 );	
	particleData.variable<float>("volume").resize( 1, 0.5f );
		
	const float gridSize = 1.0f;
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel constitutiveModel(
		1.4e3f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);

	constitutiveModel.setParticles( particleData );
	constitutiveModel.updateParticleData();

	Sim::IndexList inds;
	inds.push_back(0);
	Grid g( particleData, inds, gridSize, shapeFunction );

	// at equilibrium, so no forces on the grid nodes:
	Sim::ForceFieldSet forceFields;
	Eigen::VectorXf forces( g.m_velocities.size() );
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	float minF = fabs( forces.minCoeff() );
	float maxF = fabs( forces.maxCoeff() );
	assert( minF < 1.e-4 && maxF < 1.e-4 );
	
	
	// fiddle with F for one of the particles to get it out of equilibrium:
	std::vector<Eigen::Matrix3f>& F = particleData.variable<Matrix3f>("F");
	F[0] += 0.01f * Eigen::Matrix3f::Random();
	constitutiveModel.updateParticleData();
	
	// recompute forces on grid nodes:
	g.calculateForces( forces, constitutiveModel, forceFields.fields );
	
	Eigen::Vector3f fracDimPos = ( Eigen::Vector3f::Zero() - g.m_min ) / gridSize;
	const Eigen::Vector3i& n = g.m_n;
	int centerIdx = int(fracDimPos[0]) + n[0] *( int(fracDimPos[1]) + n[1] * int(fracDimPos[2]) );
	
	Eigen::Matrix3f Forig = F[0];
	
	Eigen::VectorXf forcesPerturbed( g.m_velocities.size() );
	Eigen::VectorXf forcesPerturbedNegative( g.m_velocities.size() );
	
	// now work out the forces on the grid nodes by finite differences!
	std::vector<float>& volume = particleData.variable<float>("volume");
	float e0 = constitutiveModel.energyDensity( 0 ) * volume[0];
	const float dx = 0.01f;
	for( int i=0; i <g.m_velocities.size(); ++i )
	{
		if( forces[i] == 0 )
		{
			continue;
		}
		
		F[0] = Forig;
		constitutiveModel.updateParticleData();
		VectorXf df(g.m_velocities.size());
		g.m_velocities[i] = dx;
		g.calculateForceDifferentials(
			df,
			g.m_velocities,
			constitutiveModel,
			forceFields.fields );

		// use the deformation gradient update to move one of the nodes
		// a small distance along one of the axes, so it squishes the
		// particle around a bit and work out the energy in this perturbed
		// configuration:
		g.m_velocities[i] = dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		constitutiveModel.updateParticleData();
		float ePlus = constitutiveModel.energyDensity( 0 ) * volume[0];
		g.calculateForces( forcesPerturbed, constitutiveModel, forceFields.fields );
		
		// now deform it the other way and work out the energy:
		g.m_velocities[i] = -dx;
		F[0] = Forig;
		g.updateDeformationGradients( 1.0f );
		constitutiveModel.updateParticleData();
		float eMinus = constitutiveModel.energyDensity( 0 ) * volume[0];
		g.calculateForces( forcesPerturbedNegative, constitutiveModel, forceFields.fields );
		
		// x component of force on this node is -dE/dx, same for y and z. Make fd and analytic values
		// agree within a percent:
		assert( fabs( ( forces[i] + ( ePlus - eMinus ) / ( 2 * dx ) ) / forces[i] ) < 0.01f );
		
		// I wish I knew why this stuff converged so badly... looks like I can only get the force
		// differentials to agree with the finite difference approximation to within 2.5% or something.
		std::cerr << ( 0.5f * ( forcesPerturbed - forcesPerturbedNegative ) - df ).norm() / df.norm() << " " << df.norm() << std::endl;
		assert( ( 0.5f * ( forcesPerturbed - forcesPerturbedNegative ) - df ).norm() / df.norm() < 0.025f );

		// put things back in their place
		g.m_velocities[i] = 0;
	}
}

class TestGrid::ImplicitUpdateRecord : public LinearSolver::Debug
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
		for( int i=0; i < g.m_masses.size(); ++i )
		{
			if( g.m_masses[i] != 0 )
			{
				explicitVelocities.segment<3>( 3 * i ) = explicitMomenta.segment<3>( 3 * i ) / g.m_masses[i];
			}
			else
			{
				explicitVelocities.segment<3>( 3 * i ).setZero();
			}
		}
		explicitVelocities += vc;

		float gridH = g.m_gridSize;
		outFile.write( (const char*)&gridH, sizeof( float ) );

		outFile.write( (const char*)&(g.m_min[0]), sizeof( float ) );
		outFile.write( (const char*)&(g.m_min[1]), sizeof( float ) );
		outFile.write( (const char*)&(g.m_min[2]), sizeof( float ) );

		outFile.write( (const char*)&(g.m_n[0]), sizeof( int ) );
		outFile.write( (const char*)&(g.m_n[1]), sizeof( int ) );
		outFile.write( (const char*)&(g.m_n[2]), sizeof( int ) );

		outFile.write( (const char*)&(g.m_masses[0]), g.m_n[0] * g.m_n[1] * g.m_n[2] * sizeof( float ) );
		
		outFile.write( (const char*)&(explicitVelocities[0]), 3 * g.m_n[0] * g.m_n[1] * g.m_n[2] * sizeof( float ) );
	}
	
	virtual void operator()( Eigen::VectorXf& x )
	{
		Eigen::VectorXf iterate = x + vc;
		outFile.write( (const char*)&(iterate[0]), (std::streamsize)(iterate.size() * sizeof( float )) );
	}
	
	Eigen::VectorXf vc;
	std::ofstream outFile;
	std::vector< Eigen::VectorXf > record;

};


void TestGrid::testImplicitUpdate()
{
	std::cerr << "testImplicitUpdate()" << std::endl;

	// create some particules:
	MaterialPointData particleData;
	
	std::vector<Vector3f>& velocities = particleData.variable<Vector3f>( "v" );
	std::vector<Vector3f>& positions = particleData.variable<Vector3f>( "p" );
	std::vector<Matrix3f>& F = particleData.variable<Matrix3f>( "F" );
	std::vector<float>& masses = particleData.variable<float>( "m" );
	std::vector<float>& volumes = particleData.variable<float>( "volume" );
	
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
	Grid::voxelSort(
		inds.begin(),
		inds.end(),
		voxelSize,
		positions );
	
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
	
	SquareMagnitudeTermination t( 40, 0.0f );
	
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
			t,
			&d
		);
	}

	//system( "C:\\Users\\david\\Documents\\GitHub\\mpmsnow\\Debug\\viewer.exe" );

	// check nothing's moving in or out of the collision objects:
	VectorXf vc( g.m_velocities.size() );
	g.collisionVelocities( vc, collisionObjects.objects, nodeCollided );
	for( int i=0; i < g.m_n[0]; ++i )
	{
		for( int j=0; j < g.m_n[1]; ++j )
		{
			for( int k=0; k < g.m_n[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				Vector3f velocity = g.m_velocities.segment<3>( 3 * idx );
				Vector3f x(
					g.m_gridSize * i + g.m_min[0],
					g.m_gridSize * j + g.m_min[1],
					g.m_gridSize * k + g.m_min[2]
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
	// P * ( M * ( g.m_velocities - vc ) - DF * g.m_velocities * dt * dt ) = explicitMomenta
	// (here, P means we project out the collided velocity components)
	
	for( int i=0; i < g.m_n[0]; ++i )
	{
		for( int j=0; j < g.m_n[1]; ++j )
		{
			for( int k=0; k < g.m_n[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				
				// apply P to velocities:
				if( nodeCollided[idx] >= 0 )
				{
					const CollisionObject* obj = collisionObjects.objects[ nodeCollided[idx] ];
					
					// find object normal:
					Vector3f x(
						g.m_gridSize * i + g.m_min[0],
						g.m_gridSize * j + g.m_min[1],
						g.m_gridSize * k + g.m_min[2]
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
	
	// work out DF * g.m_velocities:
	VectorXf df( g.m_velocities.size() );
	g.calculateForceDifferentials( df, g.m_velocities, snowModel, fields );

	//  M * ( g.m_velocities - vc ) - DF * g.m_velocities * dt * dt 
	VectorXf lhs( g.m_velocities.size() );
	for( int idx=0; idx < g.m_masses.size(); ++idx )
	{
		lhs.segment<3>( 3 * idx ) =
			g.m_masses[idx] * ( g.m_velocities.segment<3>( 3 * idx ) - vc.segment<3>( 3 * idx ) ) - df.segment<3>( 3 * idx ) * timeStep * timeStep;
	}
	
	// run the projection again:
	for( int i=0; i < g.m_n[0]; ++i )
	{
		for( int j=0; j < g.m_n[1]; ++j )
		{
			for( int k=0; k < g.m_n[2]; ++k )
			{
				int idx = g.coordsToIndex( i, j, k );
				
				// apply P:
				if( nodeCollided[idx] >= 0 )
				{
					const CollisionObject* obj = collisionObjects.objects[ nodeCollided[idx] ];
					
					// find object normal:
					Vector3f x(
						g.m_gridSize * i + g.m_min[0],
						g.m_gridSize * j + g.m_min[1],
						g.m_gridSize * k + g.m_min[2]
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
	for( int idx=0; idx < g.m_masses.size(); ++idx )
	{
		float norm = ( lhs.segment<3>( 3 * idx ) - explicitMomenta.segment<3>( 3 * idx ) ).norm();
		if( norm > maxNorm )
		{
			maxNorm = norm;
		}
	}
	assert( maxNorm < 1.e-6 );
}

void TestGrid::testMovingGrid()
{
	std::cerr << "testMovingGrid" << std::endl;
	
	// create some particules:
	MaterialPointData particleData;

	std::vector<Vector3f>& velocities = particleData.variable<Vector3f>( "v" );
	std::vector<Vector3f>& positions = particleData.variable<Vector3f>( "p" );
	std::vector<Matrix3f>& F = particleData.variable<Matrix3f>( "F" );
	std::vector<float>& masses = particleData.variable<float>( "m" );
	std::vector<float>& volumes = particleData.variable<float>( "volume" );
	
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
	Grid::voxelSort(
		inds.begin(),
		inds.end(),
		voxelSize,
		positions );
	
	snowModel.setParticles( particleData );
	Grid g( particleData, inds, gridSize, shapeFunction, Vector3f( 0, -1, 0 ) );
	g.computeParticleVolumes();
	
	float timeStep = 0.002f;
	Sim::CollisionObjectSet collisionObjects;
	CollisionPlane* plane1 = new CollisionPlane( Eigen::Vector4f( 0,1,0,1 ) );
	collisionObjects.objects.push_back( plane1 );
	
	std::vector<const ForceField*> fields;
	SquareMagnitudeTermination t( 40, 0.0f );
	
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
			t,
			&d
		);
	}
	system( "C:\\Users\\david\\Documents\\GitHub\\mpmsnow\\Debug\\viewer.exe" );


}

void TestGrid::testDfiDxi()
{

	std::cerr << "testForces()" << std::endl;
	// create a single particle:
	MaterialPointData particleData;

	particleData.variable<Vector3f>( "v" ).push_back( Vector3f::Zero() );
	particleData.variable<Vector3f>( "p" ).push_back( Vector3f::Zero() );
	particleData.variable<Matrix3f>( "F" ).push_back( Matrix3f::Identity() + 0.1f * Matrix3f::Random() );
	particleData.variable<float>( "m" ).push_back( 1.0f );
	particleData.variable<float>( "volume" ).push_back( 0.5f );;
	
	const float gridSize = 1.0f;
	CubicBsplineShapeFunction shapeFunction;
	SnowConstitutiveModel constitutiveModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		0, // hardening
		100000.0f, // compressive strength
		100000.0f	// tensile strength
	);

	constitutiveModel.setParticles( particleData );
	constitutiveModel.updateParticleData();

	Sim::IndexList inds;
	inds.push_back(0);
	Grid g( particleData, inds, gridSize, shapeFunction );
	Eigen::VectorXf dfidxi;
	g.dForceidXi( dfidxi, constitutiveModel );

	std::vector<const ForceField*> fields;
	Eigen::VectorXf dx = Eigen::VectorXf::Constant(dfidxi.size(),0);
	Eigen::VectorXf df( dfidxi.size() );
	for( int i=0; i < dfidxi.size(); ++i )
	{
		dx[i] = 1;
		g.calculateForceDifferentials(
			df,
			dx,
			constitutiveModel,
			fields );
		assert( fabs( df[i] - dfidxi[i] ) < 1.e-8 );
		dx[i] = 0;
	}
}

void TestGrid::test()
{
	std::cerr << "testGrid()" << std::endl;
	testProcessingPartitions();
	testSplatting();
	testDeformationGradients();
	testForces();
	testImplicitUpdate();
	testMovingGrid();
	testDfiDxi();
}

}
