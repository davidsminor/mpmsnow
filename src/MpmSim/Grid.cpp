#include "MpmSim/Grid.h"
#include "MpmSim/ForceField.h"
#include "MpmSim/ImplicitUpdateMatrix.h"

#include <iostream>
#include <stdexcept>

using namespace MpmSim;
using namespace Eigen;

Grid::GridSplatter::GridSplatter(
	const Grid& g,
	Eigen::VectorXf& result
) :
	m_g( g ),
	m_result( result )
{
}

void Grid::GridSplatter::setPartition( int i, int j, int k )
{
	m_partition = &m_g.m_processingPartitions[i][j][k];
}

void Grid::GridSplatter::operator()(const tbb::blocked_range<int> &r) const
{
	// iterate over a range of voxels assigned by tbb:
	for (int i = r.begin(); i != r.end(); ++i)
	{
		// splat all the particles in this voxel:
		splat( (*m_partition)[i].first, (*m_partition)[i].second, m_result );
	}
}

class MassSplatter : public Grid::GridSplatter
{
public:
	MassSplatter( const Grid& g, Eigen::VectorXf& result )
		:
		Grid::GridSplatter( g, result ),
		m_particleX( particleVariable<VectorData>("p")->m_data ),
		m_particleM( particleVariable<ScalarData>("m")->m_data )
	{
	}
	
	virtual void splat(
		Sim::ConstIndexIterator begin,
		Sim::ConstIndexIterator end,
		Eigen::VectorXf& gridMasses
	) const
	{
		Grid::ShapeFunctionIterator& shIt = m_g.shapeFunctionIterator();
		Vector3i particleCell;
		for( Sim::ConstIndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			shIt.initialize( m_particleX[p] );
			do
			{
				shIt.gridPos( particleCell );
				int idx = m_g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
				gridMasses[ idx ] += m_particleM[p] * shIt.w();
			} while( shIt.next() );
		}
	}

private:

	const std::vector<Eigen::Vector3f>& m_particleX;
	const std::vector<float>& m_particleM;

};


class VelocitySplatter : public Grid::GridSplatter
{

public:

	VelocitySplatter( const Grid& g, Eigen::VectorXf& result )
		:
		Grid::GridSplatter( g, result ),
		m_particleM( particleVariable<ScalarData>("m")->m_data ),
		m_particleX( particleVariable<VectorData>("p")->m_data ),
		m_particleV( particleVariable<VectorData>("v")->m_data )
	{	
	}
	
	virtual void splat(
		Sim::ConstIndexIterator begin,
		Sim::ConstIndexIterator end,
		Eigen::VectorXf& gridVelocities
	) const
	{
		Grid::ShapeFunctionIterator& shIt = m_g.shapeFunctionIterator();
		Vector3i particleCell;
		const Eigen::VectorXf& masses = m_g.masses;
		for( Sim::ConstIndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			shIt.initialize( m_particleX[p] );
			
			// splat the velocities:
			do
			{
				shIt.gridPos( particleCell );
				int idx = m_g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
				float gridCellMass = masses[idx];
				if( gridCellMass > 0 )
				{
					gridVelocities.segment<3>( 3 * idx ) += shIt.w() * ( m_particleM[p] / gridCellMass ) * ( m_particleV[p] - m_g.m_frameVelocity );
				}
			} while( shIt.next() );
		}
	}
	
private:

	const std::vector<float>& m_particleM;
	const std::vector<Eigen::Vector3f>& m_particleX;
	const std::vector<Eigen::Vector3f>& m_particleV;

};

class ForceSplatter : public Grid::GridSplatter
{
public:
	
	ForceSplatter(
		const Grid& g,
		Eigen::VectorXf& result,
		const ConstitutiveModel& constitutiveModel
	)
		:
		Grid::GridSplatter( g, result ),
		m_particleX( particleVariable<VectorData>("p")->m_data ),
		m_particleF( particleVariable<MatrixData>("F")->m_data ),
		m_particleVolumes( particleVariable<ScalarData>("volume")->m_data ),
		m_constitutiveModel( constitutiveModel )
	{
	}

	virtual void splat(
		Sim::ConstIndexIterator begin,
		Sim::ConstIndexIterator end,
		Eigen::VectorXf& forces
	) const
	{
		Grid::ShapeFunctionIterator& shIt = m_g.shapeFunctionIterator();
		Vector3i particleCell;
		Vector3f weightGrad;

		for( Sim::ConstIndexIterator it = begin; it != end; ++it )
		{
			int p = *it;
			Eigen::Matrix3f forceMatrix = m_particleVolumes[p] * m_constitutiveModel.dEnergyDensitydF( p ) * m_particleF[p].transpose();
			shIt.initialize( m_particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = m_g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
				forces.segment<3>( 3 * idx ) -= forceMatrix * weightGrad;
			} while( shIt.next() );
		}
	}	
private:

	const std::vector<float>& m_particleVolumes;
	const std::vector<Eigen::Vector3f>& m_particleX;
	const std::vector<Eigen::Matrix3f>& m_particleF;
	const ConstitutiveModel& m_constitutiveModel;
};


class ForceDifferentialSplatter : public Grid::GridSplatter
{
public:
	
	ForceDifferentialSplatter(
		const Grid& g,
		Eigen::VectorXf& result,
		const ConstitutiveModel& constitutiveModel,
		const Eigen::VectorXf& dx
	)
		:
		Grid::GridSplatter( g, result ),
		m_particleX( particleVariable<VectorData>("p")->m_data ),
		m_particleF( particleVariable<MatrixData>("F")->m_data ),
		m_particleVolumes( particleVariable<ScalarData>("volume")->m_data ),
		m_constitutiveModel( constitutiveModel ),
		m_dx( dx )
	{
	}
	
	virtual void splat(
		Sim::ConstIndexIterator begin,
		Sim::ConstIndexIterator end,
		Eigen::VectorXf& df
	) const
	{
		Grid::ShapeFunctionIterator& shIt = m_g.shapeFunctionIterator();
		Vector3i particleCell;
		Vector3f weightGrad;
		for( Sim::ConstIndexIterator it = begin; it != end; ++it )
		{
			int p = *it;

			// work out deformation gradient differential for this particle when grid nodes are
			// all moved by their respective dx
			Matrix3f dFp = Matrix3f::Zero();
			shIt.initialize( m_particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = m_g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
				dFp += m_dx.segment<3>( 3 * idx ) * weightGrad.transpose() * m_particleF[p];
			} while( shIt.next() );
			
			Matrix3f forceMatrix =
				m_particleVolumes[p] *
				m_constitutiveModel.dEdFDifferential( dFp, p ) *
				m_particleF[p].transpose();
			
			shIt.initialize( m_particleX[p], true );
			do
			{
				shIt.gridPos( particleCell );
				shIt.dw( weightGrad );
				int idx = m_g.coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
				
				// add on difference in velocity due to this force:
				df.segment<3>( 3 * idx ) -= forceMatrix * weightGrad;
						
			} while( shIt.next() );
		}
	}
	
	const std::vector<float>& m_particleVolumes;
	const std::vector<Eigen::Vector3f>& m_particleX;
	const std::vector<Eigen::Matrix3f>& m_particleF;
	const Eigen::VectorXf& m_dx;
	const ConstitutiveModel& m_constitutiveModel;
};



Grid::Grid(
		Sim::MaterialPointDataMap& d,
		const Sim::IndexList& particleInds,
		float gridSize,
		const ShapeFunction& shapeFunction,
		Eigen::Vector3f frameVelocity,
		int dimension
) :
	m_d( d ),
	m_particleInds( particleInds ),
	m_gridSize( gridSize ),
	m_shapeFunction( shapeFunction ),
	m_frameVelocity( frameVelocity ),
	m_dimension( dimension )
{
	// work out the physical size of the grid:
	m_min.setConstant( 1.e10 );
	m_max.setConstant( -1.e10 );

	const std::vector<Eigen::Vector3f>& particleX = dynamic_cast< VectorData* >( d["p"] )->m_data;

	for( Sim::ConstIndexIterator it = particleInds.begin(); it != particleInds.end(); ++it )
	{
		int p = *it;
		for( int j=0; j < m_dimension; ++j )
		{
			minMax( particleX[p][j], m_min[j], m_max[j] );
		}
	}
	
	// calculate grid dimensions and quantize bounding box:
	for( int j=0; j < m_dimension; ++j )
	{
		m_n[j] = fixDim( m_min[j], m_max[j] );
	}
	for( int j=m_dimension; j < 3; ++j )
	{
		m_min[j] = m_max[j] = 0;
		m_n[j] = 1;
	}
	
	long long ncells = (long long)m_n[0] * (long long)m_n[1] * (long long)m_n[2];
	if( ncells > 4000000000 || ncells <= 0 )
	{
		throw std::runtime_error( "grid is too big" );
	}
	
	// partition the particle inds for paralell processing:
	computeProcessingPartitions();

	// calculate masses:
	masses.resize( ncells );
	masses.setZero();
	
	MassSplatter sM( *this, masses );
	splat( sM );
	
	// grid masses can end up less than zero due to numerical issues in the shape functions, so clamp 'em:
	for( int i=0; i < masses.size(); ++i )
	{
		if( masses[i] < 0 )
		{
			masses[i] = 0;
		}
	}

	// calculate velocities:
	velocities.resize( ncells * 3 );
	velocities.setZero();
	
	VelocitySplatter sV( *this, velocities );
	splat( sV );
	
	prevVelocities = velocities;
}

const Eigen::Vector3f& Grid::minCoord() const
{
	return m_min;
}

const Eigen::Vector3f& Grid::maxCoord() const
{
	return m_max;
}

const Eigen::Vector3i& Grid::n() const
{
	return m_n;
}


float Grid::gridSize() const
{
	return m_gridSize;
}

void Grid::computeParticleVolumes() const
{
	const std::vector<Eigen::Vector3f>& particleX = dynamic_cast< VectorData* >( m_d["p"] )->m_data;
	const std::vector<float>& particleM = dynamic_cast< ScalarData* >( m_d["m"] )->m_data;
	std::vector<float>& particleVolumes = dynamic_cast< ScalarData* >( m_d["volume"] )->m_data;
	
	ShapeFunctionIterator& shIt = shapeFunctionIterator();
	Vector3i particleCell;
	
	float cellVolume = m_gridSize * m_gridSize * m_gridSize;
	for( Sim::ConstIndexIterator it = m_particleInds.begin(); it != m_particleInds.end(); ++it )
	{
		float density(0);
		int p = *it;
		shIt.initialize( particleX[p] );
		do
		{
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
			
			// accumulate the particle's density:
			density += shIt.w() * masses[ idx ] / cellVolume;
					
		} while( shIt.next() );
		
		particleVolumes[p] = particleM[p] / density;
	}
}

const Grid::PartitionList& Grid::partition( int i, int j, int k ) const
{
	return m_processingPartitions[i][j][k];
}

void Grid::computeProcessingPartitions()
{
	const std::vector<Eigen::Vector3f>& particleX = dynamic_cast< VectorData* >( m_d["p"] )->m_data;
	
	// now imagine chopping space up into little 2x2x2 voxel blocks. All
	// the voxels in the (0,0,0) corners go in processingPartitions[0][0][0],
	// all the voxels in the (1,0,0) corners go in processingPartitions[1][0][0],
	// etc etc.
	Eigen::Vector3i currentVoxel;
	Sim::ConstIndexIterator begin = m_particleInds.begin();
	Sim::ConstIndexIterator end = m_particleInds.end();
	Sim::ConstIndexIterator* partitionEnd = 0;
	
	float voxelSize = 2 * m_shapeFunction.supportRadius() * m_gridSize;
	for( Sim::ConstIndexIterator it = begin; it != end; ++it )
	{
		Eigen::Vector3f x = particleX[ *it ] / voxelSize;
		Eigen::Vector3i voxel( (int)floor( x[0] ), (int)floor( x[1] ), (int)floor( x[2] ) );
		if( voxel != currentVoxel || it == begin )
		{
			currentVoxel = voxel;
			PartitionList& partition = m_processingPartitions[ voxel[0]&1 ][ voxel[1]&1 ][ voxel[2]&1 ];
			partition.push_back( std::make_pair( it, it ) );
			partitionEnd = &partition.back().second;
		}
		++*partitionEnd;
	}
}



void Grid::calculateForces(
	VectorXf& forces,
	const ConstitutiveModel& constitutiveModel,
	const std::vector<const ForceField*>& fields ) const
{
	forces.setZero();

	// force fields:
	for( size_t n=0; n < fields.size(); ++n )
	{
		for( int i=0; i < m_n[0]; ++i )
		{
			for( int j=0; j < m_n[1]; ++j )
			{
				for( int k=0; k < m_n[2]; ++k )
				{
					int idx = coordsToIndex( i, j, k );
					forces.segment<3>(3 * idx) += fields[n]->force( Eigen::Vector3f( float(i), float(j), float(k) ) * m_gridSize + m_min, masses[idx] );
				}
			}
		}
	}

	// add on internal forces:
	ForceSplatter s( *this, forces, constitutiveModel );
	splat( s );
}



void Grid::calculateForceDifferentials(
		VectorXf& df,
		const VectorXf& dx,
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const ForceField*>& fields
) const
{
	// TODO: this doesn't deal with force fields which vary in space
	df.resize( velocities.size() );
	df.setZero();
	ForceDifferentialSplatter s( *this, df, constitutiveModel, dx );
	splat( s );
}

void Grid::calculateExplicitMomenta(
	VectorXf& explicitMomenta,
	std::vector<char>& nodeCollided,
	float timeStep,
	const ConstitutiveModel& constitutiveModel,
	const Sim::CollisionObjectSet& collisionObjects,
	const std::vector<const ForceField*>& fields
)
{
	// work out forces on grid points:
	VectorXf forces( velocities.size() );
	calculateForces( forces, constitutiveModel, fields );
	
	// work out explicit velocity update, and convert it to momenta:
	explicitMomenta.resize( velocities.size() );
	nodeCollided.resize( masses.size() );
	for( int i=0; i < m_n[0]; ++i )
	{
		for( int j=0; j < m_n[1]; ++j )
		{
			for( int k=0; k < m_n[2]; ++k )
			{
				int idx = coordsToIndex( i, j, k );

				Vector3f force = forces.segment<3>( 3 * idx );
				Vector3f velocity = velocities.segment<3>( 3 * idx );
				Vector3f explicitVelocity = velocity;
				if( masses[idx] > 0 )
				{
					explicitVelocity += timeStep * force / masses[idx];
				}

				// which collision objects affect this node? -1 means none, -2 means more than one, >= 0
				// is the object index:
				Vector3f x( m_gridSize * i + m_min[0], m_gridSize * j + m_min[1], m_gridSize * k + m_min[2] );
				nodeCollided[idx] = collisionObjects.collide( explicitVelocity, x, m_frameVelocity );
				explicitMomenta.segment<3>( 3 * idx ) = explicitVelocity * masses[idx];
			}
		}
	}
	
}

void Grid::collisionVelocities(
	Eigen::VectorXf& vc,
	const std::vector<const CollisionObject*>& collisionObjects,
	const std::vector<char>& nodeCollided
) const
{
	vc.resize( 3 * m_n[0] * m_n[1] * m_n[2] );
	for( int i=0; i < m_n[0]; ++i )
	{
		for( int j=0; j < m_n[1]; ++j )
		{
			for( int k=0; k < m_n[2]; ++k )
			{
				int idx = coordsToIndex( i, j, k );
				if( nodeCollided[idx] < 0 )
				{
					vc.segment<3>( 3 * idx ).setZero();
				}
				else
				{
					Vector3f x( m_gridSize * i + m_min[0], m_gridSize * j + m_min[1], m_gridSize * k + m_min[2] );
					const CollisionObject* obj = collisionObjects[ nodeCollided[idx] ];
					
					Vector3f vObj;
					obj->velocity( x, vObj );
					// express collision velocity relative to moving frame:
					vc.segment<3>( 3 * idx ) = vObj - m_frameVelocity;
				}
			}
		}
	}

}

void Grid::updateGridVelocities(
	float timeStep,
	const ConstitutiveModel& constitutiveModel,
	const Sim::CollisionObjectSet& collisionObjects,
	const std::vector<const ForceField*>& fields,
	const LinearSolver& implicitSolver,
	LinearSolver::Debug* d )
{
	VectorXf explicitMomenta;
	calculateExplicitMomenta(
		explicitMomenta,
		m_nodeCollided,
		timeStep,
		constitutiveModel,
		collisionObjects,
		fields
	);
	
	// so, lets work out vc:
	VectorXf vc;
	collisionVelocities( vc, collisionObjects.objects, m_nodeCollided );
	
	// work out new centre of mass velocity, and set all velocities to that
	// to provide an initial guess for the solver:
	float totalMass = 0;
	float collidedMass = 0;
	Eigen::Vector3f totalMomentum = Eigen::Vector3f::Zero();
	Eigen::Vector3f collidedMomentum = Eigen::Vector3f::Zero();
	for( int idx=0; idx<masses.size(); ++idx )
	{
		totalMass += masses[idx];
		totalMomentum += explicitMomenta.segment<3>( 3*idx );
		if( m_nodeCollided[idx] >= 0 )
		{
			collidedMass += masses[idx];
			collidedMomentum += masses[idx] * vc.segment<3>( 3*idx );
		}
	}
	Eigen::Vector3f initialGuessVelocity;
	if( totalMass == 0 )
	{
		initialGuessVelocity.setZero();
	}
	else if( collidedMass == 0 )
	{
		// use centre of mass velocity as an initial guess:
		initialGuessVelocity = totalMomentum / totalMass;
	}
	else
	{
		// use velocity of collision objects as an initial guess:
		initialGuessVelocity = collidedMomentum / collidedMass;
	}

	for( int idx=0; idx<masses.size(); ++idx )
	{
		velocities.segment<3>( 3*idx ) = initialGuessVelocity;
	}

	// looks like this works pretty well for when things are in freefall actually, although
	// it still kind of sucks for resting contact. Maybe the solver's termination criterion
	// is wrong? Perhaps I can make a little class to make that configurable?
	// It'd probably help to include angular velocity too innit.

	// My understanding of why it works is that the low frequency modes are the ones the
	// solver finds most difficult to resolve, and this gets the lowest frequency mode
	// right first time

	ImplicitUpdateMatrix implicitMatrix( m_d, *this, constitutiveModel, collisionObjects.objects, fields, timeStep );
	
	// what's the forward update gonna look like?
	// v^(n) = collide ( v^(n+1) - M^-1 f^(n+1) dt - vc ) + vc
	// vr^(n) = collide ( vr^(n+1) - M^-1 f^(n+1) dt )
	// vr^(n) = collide ( vr^(n+1) - M^-1 F( x^n+1 ) dt )
	// vr^(n) = collide ( vr^(n+1) - M^-1 F( x^n + v^(n+1) dt ) dt )
	// vr^(n) = collide ( vr^(n+1) - M^-1 * F( x^n ) dt - M^-1 * DF * v^(n+1) * dt * dt )
	
	// We gotta solve this for vr^(n+1), then add of vc to get v^(n+1).
	// collide( v ) isn't linear though. What we'll do, (and I should really work out if this
	// is actually justified shouldn't I...) is just brainlessly take shit over the other side:

	// collide ( vr^(n) + M^-1 * F( x^n ) dt ) = P * ( vr^(n+1) - M^-1 * DF * v^(n+1) * dt * dt )
	// collide ( vr^(n) + M^-1 * F( x^n ) dt ) = P * ( vr^(n+1) - M^-1 * DF * ( vr^(n+1) + vc ) * dt * dt )
	// M * collide ( vr^(n) + M^-1 * F( x^n ) dt ) = P * ( M * vr^(n+1) - DF * ( vr^(n+1) + vc ) * dt * dt )
	// ... = P * ( M * vr^(n+1) - DF * vr^(n+1) * dt * dt - DF * vc * dt * dt )
	// ... = P * ( M * vr^(n+1) - DF * vr^(n+1) * dt * dt ) - P * DF * vc * dt * dt

	// P * ( M - DF * dt * dt ) * vr^(n+1) = M * collide ( vr^(n) + M^-1 * F( x^n ) dt ) + P * DF * vc * dt * dt
	// implicitMatrix * vr^(n+1) = explicitMomenta + P * DF * vc * dt * dt
	
	// todo: I guess we need to get the semi implicit stuff working too

	// work out the P * DF * vc * dt * dt term and add it onto explicitMomenta:
	VectorXf df( velocities.size() );
	calculateForceDifferentials( df, vc, constitutiveModel, fields );
	implicitMatrix.subspaceProject( df );

	// so subtract the extra term onto the explicit momenta:
	for( int idx=0; idx<masses.size(); ++idx )
	{
		explicitMomenta.segment<3>( 3 * idx ) +=
			timeStep * timeStep * df.segment<3>( 3 * idx );
	}
	
	
	// solve the linear system for the velocities relative to the collision objects:
	implicitSolver(
		implicitMatrix,
		explicitMomenta,
		velocities,
		d );
	
	// work out velocities relative to the grid:
	velocities += vc;
	
}

void Grid::updateDeformationGradients( float timeStep )
{
	
	const std::vector<Eigen::Vector3f>& particleX = dynamic_cast< VectorData* >( m_d["p"] )->m_data;
	std::vector<Eigen::Matrix3f>& particleF = dynamic_cast< MatrixData* >( m_d["F"] )->m_data;

	ShapeFunctionIterator& shIt = shapeFunctionIterator();
	Vector3f weightGrad;
	Vector3i particleCell;
	Matrix3f delV;
	
	Sim::ConstIndexIterator it = m_particleInds.begin();
	Sim::ConstIndexIterator end = m_particleInds.end();
	for( ; it != end; ++it )
	{
		int p = *it;
		delV.setZero();
		shIt.initialize( particleX[p], true );
		do
		{
			shIt.gridPos( particleCell );
			shIt.dw( weightGrad );
			int idx = coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
			delV += velocities.segment<3>( 3 * idx ) * weightGrad.transpose();
		} while( shIt.next() );
		
		Matrix3f newParticleF = ( Matrix3f::Identity() + timeStep * delV ) * particleF[p];
		particleF[p] = newParticleF;
	}
}

void Grid::updateParticleVelocities()
{
	const std::vector<Eigen::Vector3f>& particleX = dynamic_cast< VectorData* >( m_d["p"] )->m_data;
	std::vector<Eigen::Vector3f>& particleV = dynamic_cast< VectorData* >( m_d["v"] )->m_data;
	
	Grid::ShapeFunctionIterator& shIt = shapeFunctionIterator();
	Vector3i particleCell;
	
	const float alpha = 0.95f;

	// blend FLIP and PIC, as pure FLIP allows spurious particle motion
	// inside the cells:
	Sim::ConstIndexIterator it = m_particleInds.begin();
	Sim::ConstIndexIterator end = m_particleInds.end();
	for( ; it != end; ++it )
	{
		int p = *it;

		Vector3f vFlip = particleV[p];
		Vector3f vPic = Vector3f::Zero();
		shIt.initialize( particleX[p] );
		do
		{
			// soo... should I be interpolating momentum here instead? Need to experiment...
			shIt.gridPos( particleCell );
			int idx = coordsToIndex( particleCell[0], particleCell[1], particleCell[2] );
			float w = shIt.w();
			vFlip += w * ( velocities.segment<3>( 3 * idx ) - prevVelocities.segment<3>( 3 * idx ) );
			vPic += w * ( velocities.segment<3>( 3 * idx ) + m_frameVelocity );
		} while( shIt.next() );
		particleV[p] = alpha * vFlip + ( 1.0f - alpha ) * vPic;
	}
}

inline void Grid::minMax( float x, float& min, float& max )
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

inline int Grid::fixDim( float& min, float& max ) const
{
	int cellMin = int( floor( min / m_gridSize ) ) - m_shapeFunction.supportRadius() - 1;
	int cellMax = int( ceil( max / m_gridSize ) ) + m_shapeFunction.supportRadius() + 1;
	min = cellMin * m_gridSize;
	max = cellMax * m_gridSize;
	return cellMax - cellMin;
}

int Grid::coordsToIndex( int i, int j, int k ) const
{
	return i + m_n[0] * ( j + m_n[1] * k );
}

Grid::ShapeFunctionIterator::ShapeFunctionIterator( const Grid& g ) :
	m_grid( g ),
	m_diameter( 2 * g.m_shapeFunction.supportRadius() )
{
	for( int dim=0; dim < 3; ++dim )
	{
		m_w[dim].resize( m_diameter, 1 );
		m_dw[dim].resize( m_diameter, 0 );
	}
}

void Grid::ShapeFunctionIterator::initialize( const Vector3f& p, bool computeDerivatives )
{
	m_gradients = computeDerivatives;
	int r = (int)m_w[0].size() / 2;
	m_pos.setZero();
	m_base.setZero();

	for( int dim=0; dim < m_grid.m_dimension; ++dim )
	{
		float fracDimPos = ( p[dim] - m_grid.m_min[dim] ) / m_grid.m_gridSize;
		int dimPos = (int)floor( fracDimPos );
		fracDimPos -= dimPos;
		int j;
		
		m_base[ dim ] = dimPos + 1 - r;
		
		j = 1-r;
		for( int i = 0; i < m_diameter; ++i, ++j )
		{
			m_w[dim][i] = m_grid.m_shapeFunction.w( j - fracDimPos );
		}

		if( computeDerivatives )
		{
			j = 1-r;
			for( int i = 0; i < m_diameter; ++i, ++j )
			{
				m_dw[dim][i] = m_grid.m_shapeFunction.dw( j - fracDimPos ) / m_grid.m_gridSize;
			}
		}
	}
}

bool Grid::ShapeFunctionIterator::next()
{

	++m_pos[0];
	if( m_pos[0] >= m_diameter )
	{
		if( m_grid.m_dimension == 1 )
		{
			return false;
		}

		m_pos[0] = 0;
		++m_pos[1];
		
		if( m_pos[1] >= m_diameter )
		{
			
			if( m_grid.m_dimension == 2 )
			{
				return false;
			}

			m_pos[1] = 0;
			++m_pos[2];
		}
		if( m_pos[2] >= m_diameter )
		{
			return false;
		}
	}

	return true;
}

void Grid::ShapeFunctionIterator::gridPos( Eigen::Vector3i& pos ) const
{
	pos[0] = m_pos[0] + m_base[0];
	pos[1] = m_pos[1] + m_base[1];
	pos[2] = m_pos[2] + m_base[2];
}

void Grid::ShapeFunctionIterator::dw( Eigen::Vector3f& g ) const
{
	if( !m_gradients )
	{
		throw std::runtime_error( "Grid::ShapeFunctionIterator::dw(): derivatives not computed!" );
	}
	g[0] = -m_dw[0][ m_pos[0] ] *  m_w[1][ m_pos[1] ] *  m_w[2][ m_pos[2] ];
	g[1] = - m_w[0][ m_pos[0] ] * m_dw[1][ m_pos[1] ] *  m_w[2][ m_pos[2] ];
	g[2] = - m_w[0][ m_pos[0] ] *  m_w[1][ m_pos[1] ] * m_dw[2][ m_pos[2] ];
}

float Grid::ShapeFunctionIterator::w() const
{
	return m_w[0][ m_pos[0] ] * m_w[1][ m_pos[1] ] * m_w[2][ m_pos[2] ];
}

Grid::ShapeFunctionIterator& Grid::shapeFunctionIterator() const
{
	std::auto_ptr< ShapeFunctionIterator >& pIt = m_shapeFunctionIterators.local();
	if( pIt.get() == 0 )
	{
		pIt.reset( new ShapeFunctionIterator( *this ) );
	}
	return *pIt.get();
}
