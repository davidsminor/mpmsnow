#include "MpmSim/Sim.h"
#include "MpmSim/CollisionObject.h"
#include "MpmSim/ForceField.h"
#include "MpmSim/ConstitutiveModel.h"
#include "MpmSim/Grid.h"

#include <iostream>
#include <queue>
#include <stdexcept>

using namespace MpmSim;


int Sim::CollisionObjectSet::collide(
	Eigen::Vector3f& v,
	const Eigen::Vector3f& x,
	const Eigen::Vector3f& frameVelocity,
	bool addCollisionVelocity
) const
{
	int collisionObject(-1);
	bool intersectedAlready(false);
	for( size_t objIdx = 0; objIdx < objects.size(); ++objIdx )
	{
		float phi = objects[objIdx]->phi( x );
		if( phi <= 0 )
		{
			if( intersectedAlready )
			{
				// colliding with more than one object: set
				// velocity to zero and bail
				collisionObject = -2;
				v.setZero();
				break;
			}
			intersectedAlready = true;

			// intersecting the object
			Eigen::Vector3f vObj;
			objects[objIdx]->velocity( x, vObj );
			
			// express object velocity relative to moving frame:
			vObj -= frameVelocity;

			// subtract off object velocity:
			v -= vObj;
			
			Eigen::Vector3f n;
			objects[objIdx]->grad( x, n );
			n.normalize();
			
			float nDotV = n.dot( v );
			if( nDotV < 0 )
			{
				// trying to move into the object:
				collisionObject = (int)objIdx;

				if( objects[objIdx]->sticky() )
				{
					v.setZero();
				}
				else
				{
					
					// velocity perpendicular to the object
					Eigen::Vector3f vPerp = nDotV * n;
					
					// remaining component is velocity paralell to the object:
					Eigen::Vector3f vTangent = v - vPerp;
					float vtNorm = vTangent.norm();
					float coulombFriction = objects[objIdx]->coulombFriction();
					if( vtNorm >= -nDotV * coulombFriction )
					{
						v = vTangent * ( 1 + coulombFriction * nDotV / vTangent.norm() );
					}
					else
					{
						v.setZero();
					}
				}
			}
			
			if( addCollisionVelocity )
			{
				v += vObj;
			}
		}
	}
	
	return collisionObject;
}



Sim::Sim(
	const std::vector<Eigen::Vector3f>& x,
	const std::vector<float>& masses,
	float gridSize,
	const ShapeFunction& shapeFunction,
	const ConstitutiveModel& model,
	const CollisionObjectSet& collisionObjects,
	const ForceFieldSet& forceFields,
	int dimension
) :
	m_gridSize( gridSize ),
	m_shapeFunction( shapeFunction ),
	m_constitutiveModel( model ),
	m_collisionObjects( collisionObjects ),
	m_forceFields( forceFields ),
	m_dimension( dimension )
{
	VectorData* p = new VectorData;
	p->m_data = x;
	particleData["p"] = p;
	
	VectorData* v = new VectorData( x.size(), Eigen::Vector3f::Zero() );
	particleData["v"] = v;

	MatrixData* f = new MatrixData( x.size(), Eigen::Matrix3f::Identity() );
	particleData["F"] = f;
	
	ScalarData* m = new ScalarData;
	m->m_data = masses;
	particleData["m"] = m;

	ScalarData* volume = new ScalarData( x.size(), 0.0f );
	particleData["volume"] = volume;
	
	m_constitutiveModel.createParticleData( particleData );
	m_constitutiveModel.setParticles( particleData );

	calculateBodies();

	for( std::vector< IndexList >::iterator it = m_bodies.begin(); it != m_bodies.end(); ++it )
	{
		IndexList& b = *it;
		Grid g( particleData, b, gridSize, shapeFunction, Eigen::Vector3f::Zero(),m_dimension );
		g.computeParticleVolumes();
	}
}

Sim::~Sim()
{
	for( MaterialPointDataMap::iterator it = particleData.begin(); it != particleData.end(); ++it )
	{
		delete it->second;
	}
}


void Sim::advance( float timeStep, const LinearSolver& solver, LinearSolver::Debug* d )
{
	std::vector<Eigen::Vector3f>& particleX = particleVariable<VectorData>( "p" )->m_data;
	std::vector<Eigen::Vector3f>& particleV = particleVariable<VectorData>( "v" )->m_data;
	std::vector<float>& particleMasses = particleVariable<ScalarData>( "m" )->m_data;
	
	// advance ballistic particle velocities:
	std::cerr << m_ballisticParticles.size() << " ballistic" << std::endl;
	IndexIterator ballisticEnd = m_ballisticParticles.end();
	for( IndexIterator it = m_ballisticParticles.begin(); it != ballisticEnd; ++it )
	{
		int p = *it;
		Eigen::Vector3f accn = Eigen::Vector3f::Zero();
		for( size_t i=0; i < m_forceFields.fields.size(); ++i )
		{
			accn += m_forceFields.fields[i]->force( particleX[p], particleMasses[p] );
		}
		accn /= particleMasses[p];
		particleV[p] += timeStep * accn;
	}
	// \todo: apply collisions. Or can I just apply them all at once at the end?
	
	// update velocities on particles in material point bodies:
	BodyIterator bodyEnd = m_bodies.end();
	std::cerr << m_bodies.size() << " bodies" << std::endl;
	for( BodyIterator bIt = m_bodies.begin(); bIt != bodyEnd; ++bIt )
	{
		// find centre of mass velocity, so we can make a comoving grid:
		// \todo: angular velocity sounds worthwhile (possibly more so than linear)
		// although more fiddly
		Eigen::Vector3f centreOfMassVelocity = Eigen::Vector3f::Zero();
		
		float mass = 0;
		for( IndexIterator it = bIt->begin(); it != bIt->end(); ++it )
		{
			centreOfMassVelocity += particleV[*it] * particleMasses[*it];
			mass += particleMasses[*it];
		}
		centreOfMassVelocity /= mass;
		
		// construct comoving background grid for this body:
		Grid g( particleData, *bIt, m_gridSize, m_shapeFunction, centreOfMassVelocity, m_dimension );
		
		// update grid velocities using internal stresses...
		g.updateGridVelocities(
			timeStep,
			m_constitutiveModel,
			m_collisionObjects,
			m_forceFields.fields,
			solver,
			d
		);
		
		// transfer the grid velocities back onto the particles:
		g.updateParticleVelocities();
		// \todo: apply collisions. Or can I just apply them all at once at the end?
		
		// update particle deformation gradients:
		g.updateDeformationGradients( timeStep );
		m_constitutiveModel.updateParticleData( particleData );
		
	}
	
	// advance particle positions:
	std::vector<Eigen::Vector3f>::iterator end = particleX.end();
	std::vector<Eigen::Vector3f>::iterator it = particleX.begin();
	std::vector<Eigen::Vector3f>::iterator vIt = particleV.begin();
	for( ; it != end; ++it, ++vIt )
	{
		// resolve object collisions:
		m_collisionObjects.collide( *vIt, *it, Eigen::Vector3f::Zero(), true );
		*it += *vIt * timeStep;
	}

	calculateBodies();
}

size_t Sim::numParticleVariables() const
{
	return particleData.size();
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
		IndexList& b = m_bodies.back();
		std::queue<int> flood;
		b.push_back( (int)i );
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

			b.push_back( current );
			processed[current] = true;
			n.neighbours( particleX[current], neighbourInds );
			for( size_t j=0; j < neighbourInds.size(); ++j )
			{
				flood.push(neighbourInds[j]);
			}
		}

		// sort the spatial index so particles in the same voxel are adjacent:
		voxelSort( b.begin(), b.end(), 2 * m_shapeFunction.supportRadius() * m_gridSize, particleX );
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

Sim::CollisionObjectSet::~CollisionObjectSet()
{
	for( size_t i=0; i < objects.size(); ++i )
	{
		delete objects[i];
	}
}

void Sim::CollisionObjectSet::add( CollisionObject* o )
{
	objects.push_back( o );
}

Sim::ForceFieldSet::~ForceFieldSet()
{
	for( size_t i=0; i < fields.size(); ++i )
	{
		delete fields[i];
	}
}

void Sim::ForceFieldSet::add( ForceField* f )
{
	fields.push_back( f );
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
	return m_bodies[n];
}

// particle indices for the ballistic particles:
const Sim::IndexList& Sim::ballisticParticles() const
{
	return m_ballisticParticles;
}

