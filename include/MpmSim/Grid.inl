#ifndef MPMSIM_GRID_INL
#define MPMSIM_GRID_INL

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

namespace MpmSim
{

template <class Splatter>
void Grid::splat( const ParticleData& d, Eigen::VectorXf& result, const void* args ) const
{
	for( int i=0; i < 8; ++i )
	{
		const ParticleData::PartitionList& partition = d.processingPartitions[ i&1 ][ (i&2) / 2][ (i&4) / 4];
		Splatter s( *this, d, partition, result, args );
		tbb::parallel_for( tbb::blocked_range<int>(0, (int)partition.size()), s );
	}
}


} //namespace MpmSim

#endif // MPMSIM_GRID_INL
