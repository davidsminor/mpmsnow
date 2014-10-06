#ifndef MPMSIM_GRID_H
#define MPMSIM_GRID_H

#include "tbb/enumerable_thread_specific.h"

#include "ShapeFunction.h"
#include "MaterialPointData.h"
#include "CollisionObject.h"
#include "ConstitutiveModel.h"

#include <Eigen/Dense>

namespace MpmSimTest
{
class TestGrid;
class TestShapeFunction;
}

namespace MpmSim
{

class Grid
{

public:

	Grid(
		MaterialPointData& d,
		const Sim::IndexList& particleInds,
		float gridSize,
		const ShapeFunction& shapeFunction,
		Eigen::Vector3f frameVelocity = Eigen::Vector3f::Zero(),
		int dimension = 3
	);
	
	// work out particle volumes:
	void computeParticleVolumes() const;

	// evolve grid velocities
	void updateGridVelocities(
		float timeStep, 
		const ConstitutiveModel& constitutiveModel,
		const Sim::CollisionObjectSet& collisionObjects,
		const std::vector<const ForceField*>& fields,
		TerminationCriterion& termination,
		LinearSolver::Debug* d = 0
	);
	
	// update particle deformation gradients based on grid velocities
	void updateDeformationGradients( float timeStep );
	
	// transfer grid velocities to particles:
	void updateParticleVelocities();
	
private:
	
	// classes used by updateGridVelocities() in the implicit velocity solve:
	class ImplicitUpdateMatrix;
	class DiagonalPreconditioner;

	// map grid coordinates to cell index:
	int coordsToIndex( int i, int j, int k ) const;
	
	// the ShapeFunctionIterator is a class for iterating over all grid nodes
	// which have a specified point in their support, and evaluating the corresponding
	// shape functions and their derivatives at that point.
	class ShapeFunctionIterator
	{
	public:
		
		ShapeFunctionIterator( const Grid& g );
		
		// initialize for a specified point:
		void initialize( const Eigen::Vector3f& p, bool computeDerivatives = false );

		// move to next point - returns false when it's done:
		bool next();

		// what grid point are we on?
		void gridPos( Eigen::Vector3i& pos ) const;
		
		// compute shape function weight
		float w() const;
		
		// compute shape function weight derivative
		void dw( Eigen::Vector3f& g ) const;
		
	private:

		const Grid& m_grid;

		int m_diameter;
		
		std::vector<float> m_w[3];
		std::vector<float> m_dw[3];
		Eigen::Vector3i m_pos;
		Eigen::Vector3i m_base;
		bool m_gradients;
		
	};
	
	// thread specific storage for shape function iterators
	mutable tbb::enumerable_thread_specific< std::auto_ptr< ShapeFunctionIterator > > m_shapeFunctionIterators;

	// grabs a shape function iterator for this thread
	ShapeFunctionIterator& shapeFunctionIterator() const;
	
	// class for splatting quantities like mass onto the grid in paralell, using the shape function iterator:
	class GridSplatter;

	// splat a quantity onto the grid in paralell using the supplied splatter object
	template< class Splatter >
	void splat( Splatter& s ) const;
	

	// Paralell splatting requires that the particles be partitioned into chunks that can be processed in
	// paralell with the guarantee that the regions we write onto don't overlap. We do this by dividing
	// space up into voxels twice the size of the shape function's support radius, then taking all the
	// voxels that can be indexed (2i, 2j, 2k) for integers i,j and k as the first partition, all 
	// (2i+1,2j,2k) voxels as the second partition, etc for all the 8 possible offsets. That way we can
	// process all the particles in the first partition's voxels in paralell (seeing as no shape function
	// can possibly touch more than one of them) then we move onto the second partition, etc.
	
	// sort the specified index range in such a way that contiguous particles are in the same
	// voxel as described above.
	static void voxelSort(
		Sim::IndexIterator begin,
		Sim::IndexIterator end,
		float voxelSize,
		const std::vector<Eigen::Vector3f>& particleX,
		int dim=0 );

	// This variable consists of eight lists of voxels corresponding to the eight partitions described above
	typedef std::vector< std::pair<Sim::ConstIndexIterator, Sim::ConstIndexIterator> > ParticlesInVoxelList;
	ParticlesInVoxelList m_processingPartitions[2][2][2];
	
	
	// compute m_processingPartitions:
	void computeProcessingPartitions();
		
	// splatters for transferring mass and velocity onto the grid:
	class MassSplatter;
	class VelocitySplatter;

	
	// grid physics:
	
	// calculate the forces on all the grid nodes by considering how moving each node
	// will deform the particles in its shape function's support, and alter their energy
	// content:
	class ForceSplatter;
	void calculateForces(
		Eigen::VectorXf& forces, 
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const ForceField*>& fields ) const;
	
	// calculate the change in the forces on the grid nodes if their positions are offset
	// by df. Used in the implicit solve:
	class ForceDifferentialSplatter;
	void calculateForceDifferentials(
		Eigen::VectorXf& df,
		const Eigen::VectorXf& dx,
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const ForceField*>& fields ) const;
	
	// work out momenta in next frame using explicit Euler: calculate forces in this frame,
	// multiply by the time step and add onto the existing momenta
	void calculateExplicitMomenta(
		Eigen::VectorXf& explicitMomenta,
		std::vector<char>& nodeCollided,
		float timeStep,
		const ConstitutiveModel& constitutiveModel,
		const Sim::CollisionObjectSet& collisionObjects,
		const std::vector<const ForceField*>& fields
	);
	
	// velocities of the collision objects, for the cells in which they're active:
	void collisionVelocities(
		Eigen::VectorXf& vc,
		const std::vector<const CollisionObject*>& collisionObjects,
		const std::vector<char>& nodeCollided
	) const;
	
	// calculate how much the force changes for each node/axis per unit displacement
	// for that node on that axis. I'm currently trying to use this to make the implicit
	// solve converge faster:
	class dFidXiSplatter;
	void dForceidXi(
		Eigen::VectorXf& dfdxi, 
		const ConstitutiveModel& constitutiveModel ) const;
	
	// particle info:
	MaterialPointData& m_d;	
	Sim::IndexList m_particleInds;
	
	// grid variables:
	Eigen::VectorXf m_masses;
	Eigen::VectorXf m_velocities;
	Eigen::VectorXf m_prevVelocities;
	std::vector<char> m_nodeCollided;
	
	// shape function we're using for this grid:
	const ShapeFunction& m_shapeFunction;
	
	// grid geometry:
	float m_gridSize;
	Eigen::Vector3f m_min;
	Eigen::Vector3f m_max;
	Eigen::Vector3i m_n;
	int m_dimension;
	
	// grid motion:
	Eigen::Vector3f m_frameVelocity;

	friend class MpmSimTest::TestGrid;
	friend class MpmSimTest::TestShapeFunction;
	
	
};

} // namespace MpmSim

#include "Grid.inl"

#endif // MPMSIM_GRID_H
