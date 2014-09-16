#include <algorithm>
#include <limits>

namespace MpmSim
{

template<class BaseType>
inline bool KDTree<BaseType>::Node::isLeaf() const
{
	return m_cutAxisAndLeaf==255;
}

template<class BaseType>
inline typename KDTree<BaseType>::PointIterator *KDTree<BaseType>::Node::permFirst() const
{
	return m_perm.first;
}

template<class BaseType>
inline typename KDTree<BaseType>::PointIterator *KDTree<BaseType>::Node::permLast() const
{
	return m_perm.last;
}

template<class BaseType>
inline bool KDTree<BaseType>::Node::isBranch() const
{
	return m_cutAxisAndLeaf!=255;
}

template<class BaseType>
inline unsigned char KDTree<BaseType>::Node::cutAxis() const
{
	return m_cutAxisAndLeaf;
}

template<class BaseType>
inline BaseType KDTree<BaseType>::Node::cutValue() const
{
	return m_cutValue;
}

template<class BaseType>
inline void KDTree<BaseType>::Node::makeLeaf( PermutationIterator permFirst, PermutationIterator permLast )
{
	m_cutAxisAndLeaf = 255;
	m_perm.first = &(*permFirst);
	m_perm.last = &(*permLast);
}

template<class BaseType>
inline void KDTree<BaseType>::Node::makeBranch( unsigned char cutAxis, BaseType cutValue )
{
	m_cutAxisAndLeaf = cutAxis;
	m_cutValue = cutValue;
}

template<class BaseType>
class KDTree<BaseType>::AxisSort
{
	public :
		AxisSort( unsigned int axis ) : m_axis( axis )
		{
		}

		bool operator() ( PointIterator i, PointIterator j )
		{
			return (*i)[m_axis] < (*j)[m_axis];
		}

	private :
		const unsigned int m_axis;
};

// initialisation

template<class BaseType>
KDTree<BaseType>::KDTree()
{
}

template<class BaseType>
KDTree<BaseType>::KDTree( PointIterator first, PointIterator last, int maxLeafSize )
{
	init( first, last, maxLeafSize );
}

template<class BaseType>
void KDTree<BaseType>::init( PointIterator first, PointIterator last, int maxLeafSize  )
{
	m_maxLeafSize = maxLeafSize;
	m_lastPoint = last;
	m_perm.resize( last - first + 1 );
	unsigned int i=0;
	for( PointIterator it=first; it!=last; it++ )
	{
		m_perm[i++] = it;
	}

	/// \todo Can we reserve() enough space for m_nodes before doing this?
	build( rootIndex(), m_perm.begin(), m_perm.end() - 1 );
}

template<class BaseType>
unsigned char KDTree<BaseType>::majorAxis( PermutationConstIterator permFirst, PermutationConstIterator permLast )
{
	Point min, max;
	for( unsigned char i=0; i<3; i++ ) {
		min[i] = 1.e10;
		max[i] = -1.e10;
	}
	for( PermutationConstIterator it=permFirst; it!=permLast; it++ )
	{
		for( unsigned char i=0; i<3; i++ )
		{
			if( (**it)[i] < min[i] )
			{
				min[i] = (**it)[i];
			}
			if( (**it)[i] > max[i] )
			{
				max[i] = (**it)[i];
			}
		}
	}
	unsigned char major = 0;
	Point size = max - min;
	for( unsigned char i=1; i<3; i++ )
	{
		if( size[i] > size[major] )
		{
			major = i;
		}
	}
	return major;
}

template<class BaseType>
void KDTree<BaseType>::build( NodeIndex nodeIndex, PermutationIterator permFirst, PermutationIterator permLast )
{
	// make room for the new node
	if( nodeIndex>=m_nodes.size() )
	{
		m_nodes.resize( nodeIndex+1 );
	}

	if( permLast - permFirst > m_maxLeafSize )
	{
		unsigned int cutAxis = majorAxis( permFirst, permLast );
		PermutationIterator permMid = permFirst  + (permLast - permFirst)/2;
		std::nth_element( permFirst, permMid, permLast, AxisSort( cutAxis ) );
		BaseType cutValue = (**permMid)[cutAxis];
		// insert node
		m_nodes[nodeIndex].makeBranch( cutAxis, cutValue );

		build( lowChildIndex( nodeIndex ), permFirst, permMid );
		build( highChildIndex( nodeIndex ), permMid, permLast );
	}
	else
	{
		// leaf node
		m_nodes[nodeIndex].makeLeaf( permFirst, permLast );
	}
}

template<class BaseType>
unsigned int KDTree<BaseType>::nearestNeighbours( const Point &p, BaseType r, std::vector<PointIterator> &nearNeighbours ) const
{
	nearNeighbours.clear();

	nearestNeighboursWalk(rootIndex(), p, r*r, nearNeighbours );

	return (unsigned int)nearNeighbours.size();
}

template<class BaseType>
void KDTree<BaseType>::nearestNeighboursWalk( NodeIndex nodeIndex, const Point &p, BaseType r2, std::vector<PointIterator> &nearNeighbours ) const
{
	const Node &node = m_nodes[nodeIndex];
	if( node.isLeaf() )
	{
		PointIterator *permLast = node.permLast();
		for( PointIterator *perm = node.permFirst(); perm!=permLast; perm++ )
		{
			const Point &pp = **perm;
			BaseType dist2 = ( p - pp ).squaredNorm();

			if (dist2 < r2 )
			{
				nearNeighbours.push_back( *perm );
			}
		}
	}
	else
	{
		// node is a branch 
		BaseType d = p[node.cutAxis()] - node.cutValue();
		NodeIndex firstChild, secondChild;
		if( d>0.0 )
		{
			firstChild = highChildIndex( nodeIndex );
			secondChild = lowChildIndex( nodeIndex );
		}
		else
		{
			firstChild = lowChildIndex( nodeIndex );
			secondChild = highChildIndex( nodeIndex );
		}

		nearestNeighboursWalk( firstChild, p, r2, nearNeighbours );
		if( d*d < r2 )
		{
			nearestNeighboursWalk( secondChild, p, r2, nearNeighbours );
		}
	}
}

template<class BaseType>
inline typename KDTree<BaseType>::NodeIndex KDTree<BaseType>::numNodes() const
{
	return m_nodes.size();
}

template<class BaseType>
inline const typename KDTree<BaseType>::Node &KDTree<BaseType>::node( NodeIndex index ) const
{
	return m_nodes[index];
}

template<class BaseType>
inline typename KDTree<BaseType>::NodeIndex KDTree<BaseType>::rootIndex() const
{
	return 1;
}

template<class BaseType>
inline typename KDTree<BaseType>::NodeIndex KDTree<BaseType>::lowChildIndex( NodeIndex parentIndex ) const
{
	return parentIndex * 2;
}

template<class BaseType>
inline typename KDTree<BaseType>::NodeIndex KDTree<BaseType>::highChildIndex( NodeIndex parentIndex ) const
{
	return parentIndex * 2 + 1;
}

} // namespace IECore