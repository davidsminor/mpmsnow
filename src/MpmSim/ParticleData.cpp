#include "MpmSim/ParticleData.h"
#include <algorithm>
#include <stdexcept>

using namespace MpmSim;
using namespace Eigen;

ParticleData::ParticleData( const std::vector<Eigen::Vector3f>& x, const std::vector<float>& m, float g )
{
	
	if( x.size() != m.size() )
	{
		throw std::runtime_error( "ParticleData::ParticleData(): positions and masses have different sizes" );
	}
	
	gridSize = g;
	particleX = x;
	particleM = m;
	
	// an index table which permutes the particles. All particles in the same
	// voxel will end up adjacent in this table.
	spatialIndex.resize( x.size() );
	for( size_t i=0; i < x.size(); ++i )
	{
		spatialIndex[i] = (int)i;
	}
	computeProcessingPartitions();

	particleV.resize( x.size(), Eigen::Vector3f::Zero() );
	particleF.resize( x.size(), Eigen::Matrix3f::Identity() );
	
}

// used for sorting a list of points by their positions on the "dim" axis
class SpatialComparator
{
public:
	
	SpatialComparator( int dim, const std::vector< Eigen::Vector3f >& x ) :
		m_dim( dim ), m_x( x )
	{
	}
	
	bool operator()( int a, int b )
	{
		return m_x[a][m_dim] < m_x[b][m_dim];
	}

private:
	int m_dim;
	const std::vector< Eigen::Vector3f >& m_x;
};

void ParticleData::voxelSort( std::vector<int>::iterator begin, std::vector<int>::iterator end, int dim )
{
	// sort along dim:
	std::sort( begin, end, SpatialComparator( dim, particleX ) );

	if( dim == 2 )
	{
		// got to z coord - done!
		return;
	}
	
	// now chop into slices along dim, and sort each slice along dim + 1:
	int currentSlice = int( floor( particleX[*begin][dim] / ( 4 * gridSize ) ) );
	std::vector<int>::iterator sliceBegin = begin;
	std::vector<int>::iterator sliceEnd = begin + 1;
	for( ;;++sliceEnd )
	{
		int particleSlice = sliceEnd == end ? currentSlice + 1 : int( floor( particleX[*sliceEnd][dim] / ( 4 * gridSize ) ) );
		if( particleSlice != currentSlice )
		{
			voxelSort( sliceBegin, sliceEnd, dim + 1 );
			if( sliceEnd == end )
			{
				break;
			}
			sliceBegin = sliceEnd;
			currentSlice = particleSlice;
		}
	}
}


void ParticleData::computeProcessingPartitions()
{
	// sort the spatial index so particles in the same voxel are adjacent:
	voxelSort( spatialIndex.begin(), spatialIndex.end() );

	for(int i=0;i<8;++i)
	{
		processingPartitions[ i&1 ][ (i&2) / 2][ (i&4) / 4].clear();
	}

	// now imagine chopping space up into little 2x2x2 voxel blocks. All
	// the voxels in the (0,0,0) corners go in processingPartitions[0][0][0],
	// all the voxels in the (1,0,0) corners go in processingPartitions[1][0][0],
	// etc etc.
	Eigen::Vector3i currentVoxel;
	IndexIterator begin = spatialIndex.begin();
	IndexIterator end = spatialIndex.end();
	IndexIterator* partitionEnd = 0;
	for( IndexIterator it = begin; it != end; ++it )
	{
		Eigen::Vector3f x = particleX[ *it ] / ( 4 * gridSize );
		Eigen::Vector3i voxel( (int)floor( x[0] ), (int)floor( x[1] ), (int)floor( x[2] ) );
		if( voxel != currentVoxel || it == begin )
		{
			currentVoxel = voxel;
			PartitionList& partition = processingPartitions[ voxel[0]&1 ][ voxel[1]&1 ][ voxel[2]&1 ];
			partition.push_back( std::make_pair( it, it ) );
			partitionEnd = &partition.back().second;
		}
		++*partitionEnd;
	}
}

void ParticleData::advance( float timeStep )
{
	for( size_t p = 0; p < particleX.size(); ++p )
	{
		particleX[p] += particleV[p] * timeStep;
	}
	computeProcessingPartitions();
}
