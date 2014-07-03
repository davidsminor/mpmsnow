#include "tests/TestSimClass.h"

#include "MpmSim/Sim.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

void testNeighbourQuery()
{
	srand(20);
	std::vector<Eigen::Vector3f> particleX;
	for( int i=0; i < 1000; ++i )
	{
		float x = ((float)rand() / RAND_MAX);
		float y = ((float)rand() / RAND_MAX);
		float z = ((float)rand() / RAND_MAX);
		particleX.push_back( Eigen::Vector3f(x,y,z) );
	}

	Sim::NeighbourQuery n( particleX, 0.1f );
	for( size_t i=0; i < particleX.size(); ++i )
	{
		const Eigen::Vector3f& p = particleX[i];
		std::vector<int> actualNeighbours;
		for( size_t j=0; j < particleX.size(); ++j )
		{
			float squaredDist = ( particleX[j] - p ).squaredNorm();
			if( squaredDist > 0 && squaredDist < ( 0.1f * 0.1f ) )
			{
				actualNeighbours.push_back( (int)j );
			}
		}
		std::vector<int> neighbours;
		
		n.neighbours( p, neighbours );
		assert( neighbours.size() == actualNeighbours.size() );
		
		std::sort( neighbours.begin(), neighbours.end() );
		std::sort( actualNeighbours.begin(), actualNeighbours.end() );
		
		std::vector<int>::iterator it = neighbours.begin();
		std::vector<int>::iterator ait = actualNeighbours.begin();
		for( ;it != neighbours.end(); ++it, ++ait )
		{
			assert( *it == *ait );
		}
	}
}

void testBodies()
{

	std::vector<Vector3f> positions;
	std::vector<float> masses;
	const float gridSize = 0.1f;
	
	// create two cube shaped clusters of particles.
	// particle spacing = 1/4 of a cell:
	for( int i=0; i < 16; ++i )
	{
		for( int j=0; j < 16; ++j )
		{
			for( int k=0; k < 16; ++k )
			{
				positions.push_back( Vector3f( 0.25f *gridSize* i, 0.25f *gridSize* j, 0.25f *gridSize* k ) );
				masses.push_back( 1.0f );

				positions.push_back( Vector3f( 0.25f *gridSize* i + 6, 0.25f *gridSize* j, 0.25f *gridSize* k ) );
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
	
	// 4 particles per cell
	Sim sim( positions, masses, gridSize );
	
	// sensible particle variables?
	assert( sim.numParticleVariables() == 4 );
	assert( sim.particleVariable<MpmSim::ScalarData>("m") );
	assert( sim.particleVariable<MpmSim::VectorData>("v") );
	assert( sim.particleVariable<MpmSim::VectorData>("p") );
	assert( sim.particleVariable<MpmSim::MatrixData>("F") );
	
	assert( sim.particleVariable<MpmSim::ScalarData>("F") == 0 );
	assert( sim.particleVariable<MpmSim::ScalarData>("Fiddlesticks") == 0 );
	
	sim.calculateBodies();
	assert( sim.ballisticParticles().size() == 16 * 16 * 16 );
	
	assert( sim.numBodies() == 2 );
	assert( sim.body( 0 ).size() == 16 * 16 * 16 );
	assert( sim.body( 1 ).size() == 16 * 16 * 16 );

}

void testProcessingPartitions()
{
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	const float gridSize = 0.1f;

	const int n = 16;
	for( int i=0; i < n; ++i )
	{
		for( int j=0; j < n; ++j )
		{
			for( int k=0; k < n; ++k )
			{
				positions.push_back( Vector3f( 0.5f * gridSize * ( i + 0.5f ), 0.5f * gridSize * ( j + 0.5f ), 0.5f * gridSize * ( k + 0.5f ) ) );
				masses.push_back( 1.0f );
			}
		}
	}
		
	Sim sim( positions, masses, gridSize );
	sim.calculateBodies();
	
	assert( sim.numBodies() == 1 );
	assert( sim.ballisticParticles().empty() );
	
	// get the indices of body 0:
	const Sim::IndexList& spatialIndex = sim.body( 0 );
	
	// only one body, so this should contain all the particles:
	assert( spatialIndex.size() == positions.size() );
	
	// check the indices have the right spatial properties:
	Vector3i currentVoxel;
	int numVoxels(0);
	for( size_t p=0; p<spatialIndex.size(); ++p )
	{
		Vector3f x = positions[ spatialIndex[p] ] / ( 4 * gridSize );
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
	
	// now check the processing partitioning
	int numPartitionedVoxels(0);
	int numPartitionedParticles(0);
	for( int i=0; i < 2; ++i )
	{
		for( int j=0; j < 2; ++j )
		{
			for( int k=0; k < 2; ++k )
			{
				const Sim::PartitionList& partition = sim.bodyPartitionList( 0, i, j, k );
				numPartitionedVoxels += partition.size();
				for( size_t v=0; v < partition.size(); ++v )
				{
					Sim::IndexIterator begin = partition[v].first;
					Sim::IndexIterator end = partition[v].second;
					numPartitionedParticles += end - begin;
					Vector3i currentVoxel;
					for( Sim::IndexIterator it = begin; it != end; ++it )
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

void testSimClass()
{
	testNeighbourQuery();
	testBodies();
	testProcessingPartitions();
}

}
