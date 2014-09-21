#ifndef MPMSIM_GRID_INL
#define MPMSIM_GRID_INL

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

namespace MpmSim
{

template <class Splatter>
void Grid::splat( Splatter& s ) const
{
	for( int i=0; i < 8; ++i )
	{
		s.setPartition( i&1, (i&2) / 2, (i&4) / 4 );
		tbb::parallel_for( tbb::blocked_range<int>(0, (int)m_processingPartitions[i&1][(i&2) / 2][(i&4) / 4].size()), s );
	}
}


template< class T >
T* Grid::GridSplatter::particleVariable( const std::string& name )
{
	MaterialPointDataMap::iterator it = m_g.m_d.find( name );
	if( it == m_g.m_d.end() )
	{
		return 0;
	}
	return dynamic_cast<T*>( it->second );
}

} //namespace MpmSim

#endif // MPMSIM_GRID_INL
