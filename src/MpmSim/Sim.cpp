#include "MpmSim/Sim.h"

#include <iostream>
#include <queue>

using namespace MpmSim;

Sim::Sim( const std::vector<Eigen::Vector3f>& x, const std::vector<float>& masses, float gridSize ) :
	m_gridSize( gridSize )
{
	VectorData* p = new VectorData;
	p->m_data = x;
	m_pointData["p"] = p;
	
	VectorData* v = new VectorData;
	v->m_data.resize( x.size(), Eigen::Vector3f::Zero() );
	m_pointData["v"] = v;

	MatrixData* f = new MatrixData;
	f->m_data.resize( x.size(), Eigen::Matrix3f::Identity() );
	m_pointData["F"] = f;
	
	ScalarData* m = new ScalarData;
	m->m_data = masses;
	m_pointData["m"] = m;
}

Sim::~Sim()
{
	for( MaterialPointDataMap::iterator it = m_pointData.begin(); it != m_pointData.end(); ++it )
	{
		delete it->second;
	}
}

size_t Sim::numParticleVariables() const
{
	return m_pointData.size();
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

void Sim::voxelSort(
	IndexIterator begin,
	IndexIterator end,
	float voxelSize,
	const std::vector<Eigen::Vector3f>& particleX,
	int dim )
{
	// sort along dim:
	std::sort( begin, end, SpatialComparator( dim, particleX ) );

	if( dim == 2 )
	{
		// got to z coord - done!
		return;
	}
	
	// now chop into slices along dim, and sort each slice along dim + 1:
	int currentSlice = int( floor( particleX[*begin][dim] / voxelSize ) );
	std::vector<int>::iterator sliceBegin = begin;
	std::vector<int>::iterator sliceEnd = begin + 1;
	for( ;;++sliceEnd )
	{
		int particleSlice = sliceEnd == end ? currentSlice + 1 : int( floor( particleX[*sliceEnd][dim] / voxelSize ) );
		if( particleSlice != currentSlice )
		{
			voxelSort( sliceBegin, sliceEnd, voxelSize, particleX, dim + 1 );
			if( sliceEnd == end )
			{
				break;
			}
			sliceBegin = sliceEnd;
			currentSlice = particleSlice;
		}
	}
}

static void minMax( float x, float& min, float& max )
{
	if( x < min )
	{
		min = x;
	}
	if( x > max )
	{
		max = x;
	}
}
void Sim::Body::computeProcessingPartitions( const std::vector<Eigen::Vector3f>& particleX, float voxelSize )
{
	// sort the spatial index so particles in the same voxel are adjacent:
	voxelSort( particleInds.begin(), particleInds.end(), voxelSize, particleX );
	
	for(int i=0;i<8;++i)
	{
		processingPartitions[ i&1 ][ (i&2) / 2][ (i&4) / 4].clear();
	}

	// now imagine chopping space up into little 2x2x2 voxel blocks. All
	// the voxels in the (0,0,0) corners go in processingPartitions[0][0][0],
	// all the voxels in the (1,0,0) corners go in processingPartitions[1][0][0],
	// etc etc.
	Eigen::Vector3i currentVoxel;
	IndexIterator begin = particleInds.begin();
	IndexIterator end = particleInds.end();
	IndexIterator* partitionEnd = 0;
	for( IndexIterator it = begin; it != end; ++it )
	{
		Eigen::Vector3f x = particleX[ *it ] / voxelSize;
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

void Sim::calculateBodies()
{
	const VectorData* p = particleVariable<VectorData>("p");
	if( !p )
	{
		throw std::runtime_error( "Sim::calculateBodies(): couldn't find 'p' data" );
	}
	const std::vector<Eigen::Vector3f>& particleX = p->m_data;
	
	m_bodies.clear();
	m_ballisticParticles.clear();

	// build little grid based acceleration structure for neighbour queries:
	NeighbourQuery n( particleX, m_gridSize );
	
	std::vector<bool> processed( particleX.size(), false );
	std::vector<int> neighbourInds;
	for( size_t i=0; i < particleX.size(); ++i )
	{
		if( processed[i] )
		{
			continue;
		}
		processed[i] = true;
		n.neighbours( particleX[i], neighbourInds );
		if( neighbourInds.size() == 0 )
		{
			m_ballisticParticles.push_back((int)i);
			continue;
		}
		
		m_bodies.resize( m_bodies.size() + 1 );
		Body& b = m_bodies.back();
		std::queue<int> flood;
		b.particleInds.push_back( i );
		for( size_t j=0; j < neighbourInds.size(); ++j )
		{
			flood.push(neighbourInds[j]);
		}
		while( flood.size() )
		{
			int current = flood.front();
			flood.pop();
			if( processed[current] )
			{
				continue;
			}

			b.particleInds.push_back( current );
			processed[current] = true;
			n.neighbours( particleX[current], neighbourInds );
			for( size_t j=0; j < neighbourInds.size(); ++j )
			{
				flood.push(neighbourInds[j]);
			}
		}

		// partition the particles in this body in an 8 color checkerboard
		// pattern so they can be processed in paralell
		b.computeProcessingPartitions( particleX, 4 * m_gridSize );
	}
}

Sim::NeighbourQuery::NeighbourQuery( const std::vector<Eigen::Vector3f>& particleX, float r ) :
	m_particleX( particleX ), m_r(r)
{

	// we're gonna use a background grid to accelerate neighbour lookups, with a voxel
	// size of twice the grid size:
	float voxelSize( 2 * r );
	
	// sort particles into voxels:
	m_spatialSorting.resize( particleX.size() );
	for( IndexIterator it = m_spatialSorting.begin(); it != m_spatialSorting.end(); ++it )
	{
		*it = (int)( it - m_spatialSorting.begin() );
	}
	voxelSort( m_spatialSorting.begin(), m_spatialSorting.end(), voxelSize, particleX );

	// find bounding box:
	Eigen::Vector3f minCoord(1000000000, 1000000000, 1000000000 );
	Eigen::Vector3f maxCoord(-1000000000, -1000000000, -1000000000 );
	for( std::vector<Eigen::Vector3f>::const_iterator it = particleX.begin(); it != particleX.end(); ++it )
	{
		float prod = (*it)[0] * (*it)[1] * (*it)[2];
		#ifdef WIN32
		if( !_finite(prod) )
		#else
		if( isinff(prod) || isnanf(prod) )
		#endif
		{
			throw std::runtime_error( "nans in particle data!" );
		}
		minMax( (*it)[0], minCoord[0], maxCoord[0] );
		minMax( (*it)[1], minCoord[1], maxCoord[1] );
		minMax( (*it)[2], minCoord[2], maxCoord[2] );
	}
	
	// wot's that in grid coordinates?
	m_cellMin = Eigen::Vector3i(
		int( floor( minCoord[0] / voxelSize ) ),
		int( floor( minCoord[1] / voxelSize ) ),
		int( floor( minCoord[2] / voxelSize ) )
	);
	Eigen::Vector3i cellMax(
		int( floor( maxCoord[0] / voxelSize ) ),
		int( floor( maxCoord[1] / voxelSize ) ),
		int( floor( maxCoord[2] / voxelSize ) )
	);
	m_gridDims = cellMax - m_cellMin + Eigen::Vector3i( 1, 1, 1 );
	
	// make a voxel grid containing offsets into "spatialSorting":
	m_voxels.resize( m_gridDims[0] * m_gridDims[1] * m_gridDims[2] );
	std::fill( m_voxels.begin(), m_voxels.end(), -1 );
	for( IndexIterator it = m_spatialSorting.begin(); it != m_spatialSorting.end(); ++it )
	{
		const Eigen::Vector3f& pos = particleX[*it];
		Eigen::Vector3i cell(
			int( floor( pos[0] / (2 * m_r) ) - m_cellMin[0] ),
			int( floor( pos[1] / (2 * m_r) ) - m_cellMin[1] ),
			int( floor( pos[2] / (2 * m_r) ) - m_cellMin[2] )
		);
		size_t offset = voxelOffset( cell );
		if( m_voxels[ offset ] == -1 )
		{
			m_voxels[ offset ] = (int)( it - m_spatialSorting.begin() );
		}
	}
}

size_t Sim::NeighbourQuery::voxelOffset( const Eigen::Vector3i& cell ) const
{
	return cell[0] + m_gridDims[0] * ( cell[1] + m_gridDims[1] * cell[2]);
}

void Sim::NeighbourQuery::neighboursInCell(
	const Eigen::Vector3f& p,
	const Eigen::Vector3i& cell,
	std::vector<int>& neighbourInds ) const
{
	size_t currentVoxelIndex = voxelOffset( cell );
	int cellOffset = m_voxels[currentVoxelIndex];
	if( cellOffset == -1 )
	{
		return;
	}

	for( int offset = cellOffset; offset < (int)m_spatialSorting.size(); ++offset )
	{
		int particleIdx = m_spatialSorting[offset];
		const Eigen::Vector3f& thisPp = m_particleX[particleIdx];
		Eigen::Vector3i cell(
			int( floor( thisPp[0] / (2 * m_r) ) - m_cellMin[0] ),
			int( floor( thisPp[1] / (2 * m_r) ) - m_cellMin[1] ),
			int( floor( thisPp[2] / (2 * m_r) ) - m_cellMin[2] )
		);
		if( cellOffset != offset && voxelOffset( cell ) != currentVoxelIndex )
		{
			break;
		}
		
		float squaredDist = ( m_particleX[particleIdx] - p ).squaredNorm();
		if( squaredDist > 0 && squaredDist < m_r * m_r )
		{
			neighbourInds.push_back( particleIdx );
		}
	}
}

void Sim::NeighbourQuery::neighbours( const Eigen::Vector3f& p, std::vector<int>& neighbourInds ) const
{
	neighbourInds.clear();

	// search neighbours in this cell:
	float xSubPos( p[0] / (2 * m_r) - m_cellMin[0] );
	float ySubPos( p[1] / (2 * m_r) - m_cellMin[1] );
	float zSubPos( p[2] / (2 * m_r) - m_cellMin[2] );
	Eigen::Vector3i cell(
		int( floor( xSubPos ) ),
		int( floor( ySubPos ) ),
		int( floor( zSubPos ) )
	);
	xSubPos -= cell[0];
	ySubPos -= cell[1];
	zSubPos -= cell[2];
	
	// find cell search range:
	Eigen::Vector3i cellMin( cell - Eigen::Vector3i( xSubPos < 0.5f, ySubPos < 0.5f, zSubPos < 0.5f ) );
	Eigen::Vector3i cellMax( cell + Eigen::Vector3i( xSubPos > 0.5f, ySubPos > 0.5f, zSubPos > 0.5f ) );
	if( cellMin[0] < 0 )
	{
		cellMin[0] = 0;
	}
	if( cellMin[1] < 0 )
	{
		cellMin[1] = 0;
	}
	if( cellMin[2] < 0 )
	{
		cellMin[2] = 0;
	}
	if( cellMax[0] > m_gridDims[0] - 1 )
	{
		cellMax[0] = m_gridDims[0] - 1;
	}
	if( cellMax[1] > m_gridDims[1] - 1 )
	{
		cellMax[1] = m_gridDims[1] - 1;
	}
	if( cellMax[2] > m_gridDims[2] - 1 )
	{
		cellMax[2] = m_gridDims[2] - 1;
	}
	
	for( int i=cellMin[0]; i <= cellMax[0]; ++i )
	{
		for( int j=cellMin[1]; j <= cellMax[1]; ++j )
		{
			for( int k=cellMin[2]; k <= cellMax[2]; ++k )
			{
				neighboursInCell( p, Eigen::Vector3i(i,j,k), neighbourInds );
			}	
		}
	}
}

size_t Sim::numBodies() const
{
	return m_bodies.size();
}

// particle indices in body n:
const Sim::IndexList& Sim::body( size_t n ) const
{
	if( n >= m_bodies.size() )
	{
		throw std::runtime_error( "Sim::body(): index exceeds number of bodies" );
	}
	return m_bodies[n].particleInds;
}

const Sim::PartitionList& Sim::bodyPartitionList( size_t n, unsigned i, unsigned j, unsigned k ) const
{
	if( n >= m_bodies.size() )
	{
		throw std::runtime_error( "Sim::bodyPartitionList(): index exceeds number of bodies" );
	}

	if( i > 1 || j > 1 || k > 1 )
	{
		throw std::runtime_error( "Sim::bodyPartitionList(): partition index out of range" );
	}

	return m_bodies[n].processingPartitions[i][j][k];
}

// particle indices for the ballistic particles:
const Sim::IndexList& Sim::ballisticParticles() const
{
	return m_ballisticParticles;
}
