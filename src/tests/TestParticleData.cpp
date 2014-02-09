#include "tests/TestParticleData.h"

#include "MpmSim/ParticleData.h"

#include <iostream>

using namespace MpmSim;
using namespace Eigen;

namespace MpmSimTest
{

void testParticleData()
{
	std::vector<Vector3f> positions;
	std::vector<float> masses;
	for( int i=0; i < 10000; ++i )
	{
		float xr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float yr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float zr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
		positions.push_back( Vector3f( xr, yr, zr ) );
		masses.push_back( 1.0f );
	}
	
	const float gridSize = 0.01f;
	ParticleData d( positions, masses, gridSize );
	assert( d.spatialIndex.size() == positions.size() );
	Vector3i currentVoxel;
	int numVoxels(0);
	for( size_t p=0; p<d.spatialIndex.size(); ++p )
	{
		Vector3f x = d.particleX[ d.spatialIndex[p] ] / ( 4 * gridSize );
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
	
	int numPartitionedVoxels(0);
	int numPartitionedParticles(0);
	for( int i=0; i < 2; ++i )
	{
		for( int j=0; j < 2; ++j )
		{
			for( int k=0; k < 2; ++k )
			{
				const std::vector< std::pair< ParticleData::IndexIterator, ParticleData::IndexIterator > >& partition = d.processingPartitions[i][j][k];
				numPartitionedVoxels += partition.size();
				for( size_t v=0; v < partition.size(); ++v )
				{
					ParticleData::IndexIterator begin = partition[v].first;
					ParticleData::IndexIterator end = partition[v].second;
					numPartitionedParticles += end - begin;
					Vector3i currentVoxel;
					for( ParticleData::IndexIterator it = begin; it != end; ++it )
					{
						// all particles must be in the same voxel...
						// and the voxel must be in the right partition
						Vector3f x = d.particleX[ *it ] / ( 4 * gridSize );
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

}
