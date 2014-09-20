#ifndef MPMSIM_KDTREE_H
#define MPMSIM_KDTREE_H

#include <set>
#include <vector>

#include <Eigen/Dense>

namespace MpmSim
{

template< class BaseType >
class KDTree
{
	public :
	
		typedef Eigen::Matrix<BaseType,3,1> Point;
		typedef typename std::vector<Point>::iterator PointIterator;
		typedef PointIterator Iterator;

		class Node;
		typedef std::vector<Node> NodeVector;
		typedef typename NodeVector::size_type NodeIndex;

		/// Constructs an unititialised tree - you must call init()
		/// before using it.
		KDTree();

		/// Creates a tree for the fast searching of points.
		/// Note that the tree does not own the passed points -
		/// it is up to you to ensure that they remain valid and
		/// unchanged as long as the KDTree is in use.
		KDTree( PointIterator first, PointIterator last, int maxLeafSize=4 );

		/// Builds the tree for the specified points - the iterator range
		/// must remain valid and unchanged as long as the tree is in use.
		/// This method can be called again to rebuild the tree at any time.
		/// \threading This can't be called while other threads are
		/// making queries.
		void init( PointIterator first, PointIterator last, int maxLeafSize=4 );

		/// Populates the passed vector of iterators with the neighbours of point p which are closer than radius r. Returns the number of points found.
		/// \todo There should be a form where nearNeighbours is an output iterator, to allow any container to be filled.
		/// See enclosedPoints for an example of this form.
		/// \threading May be called by multiple concurrent threads provided they are each using a different vector for the result.
		unsigned int nearestNeighbours( const Point &p, BaseType r, std::vector<PointIterator> &nearNeighbours ) const;
		
		/// Returns the number of nodes in the tree.
		inline NodeIndex numNodes() const;
		/// Returns the specified Node of the tree. See rootIndex(), lowChildIndex() and highChildIndex() for
		/// means of getting appropriate indices. This can be used to implement algorithms not provided as
		/// member functions.
		inline const Node &node( NodeIndex index ) const;
		/// Returns the index for the root node.
		inline NodeIndex rootIndex() const;
		/// Returns the index for the "low" child of the specified Node. This will only
		/// be valid if parent.isBranch() is true.
		inline NodeIndex lowChildIndex( NodeIndex parentIndex ) const;
		/// Returns the index for the "high" child of the specified Node. This will only
		/// be valid if parent.isBranch() is true.
		inline NodeIndex highChildIndex( NodeIndex parentIndex ) const;

	private :

		typedef std::vector<PointIterator> Permutation;
		typedef typename Permutation::iterator PermutationIterator;
		typedef typename Permutation::const_iterator PermutationConstIterator;
		
		class Neighbour;
		class AxisSort;

		unsigned char majorAxis( PermutationConstIterator permFirst, PermutationConstIterator permLast );
		void build( NodeIndex nodeIndex, PermutationIterator permFirst, PermutationIterator permLast );

		void nearestNeighboursWalk( NodeIndex nodeIndex, const Point &p, BaseType r2, std::vector<PointIterator> &nearNeighbours ) const;

		Permutation m_perm;
		NodeVector m_nodes;
		int m_maxLeafSize;
		PointIterator m_lastPoint;

};

/// The Node class which is used to implement the branching structure in the KDTree.
template<class BaseType>
class KDTree<BaseType>::Node
{
	public :

		/// Returns true if this is a leaf node of the tree.
		inline bool isLeaf() const;
		/// Returns a pointer to an iterator referencing the first
		/// child of this Node. Only valid if isLeaf() is true.
		inline PointIterator *permFirst() const;
		/// Returns a pointer to an iterator referencing the last
		/// child of this Node. Only valid if isLeaf() is true.
		inline PointIterator *permLast() const;
		/// Returns true if this is a branch node of the tree;
		inline bool isBranch() const;
		/// Returns the axis in which this node cuts the space. Only
		/// valid if isBranch() is true.
		inline unsigned char cutAxis() const;
		/// Returns the point within cutAxis() at which the node
		/// cuts the space.
		inline BaseType cutValue() const;

	private :

		friend class KDTree<BaseType>;

		inline void makeLeaf( PermutationIterator permFirst, PermutationIterator permLast );
		inline void makeBranch( unsigned char cutAxis, BaseType cutValue );

		unsigned char m_cutAxisAndLeaf;
		union {
			BaseType m_cutValue;
			struct {
				PointIterator *first;
				/// \todo Could we just store an offset instead of last? If we limit the maxLeafSize to 255
				/// then this could be a single byte instead of 8 bytes.
				PointIterator *last;
			} m_perm;
		};

};

/// The Neighbour class is used to return information from the KDTree::nearestNNeighbours() query.
template<class BaseType>
class KDTree<BaseType>::Neighbour
{
	public :

		Neighbour( Iterator p, BaseType d2 )
			:       point( p ), distSquared( d2 )
		{
		}
		
		Iterator point;
		BaseType distSquared;
	
		bool operator < ( const Neighbour &other ) const
		{
			return distSquared < other.distSquared;
		}
		
};

typedef KDTree<float> V3fTree;
typedef KDTree<double> V3dTree;

} // namespace MpmSim

#include "KdTree.inl"

#endif // MPMSIM_KDTREE_H
