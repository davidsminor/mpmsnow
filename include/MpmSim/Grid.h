#ifndef MPMSIM_GRID_H
#define MPMSIM_GRID_H

#include "tbb/enumerable_thread_specific.h"

#include "ShapeFunction.h"
#include "Sim.h"
#include "CollisionObject.h"
#include "ConstitutiveModel.h"

#include <Eigen/Dense>

namespace MpmSim
{

class Grid
{

public:

	Grid(
		Sim::MaterialPointDataMap& d,
		const Sim::IndexList& particleInds,
		float gridSize,
		const ShapeFunction& shapeFunction,
		Eigen::Vector3f frameVelocity = Eigen::Vector3f::Zero(),
		int dimension = 3
	);
	
	// map grid coordinates to cell index:
	int coordsToIndex( int i, int j, int k ) const;
	
	// work out particle volumes:
	void computeParticleVolumes() const;

	// evolve grid velocities
	void updateGridVelocities(
		float timeStep, 
		const ConstitutiveModel& constitutiveModel,
		const Sim::CollisionObjectSet& collisionObjects,
		const std::vector<const ForceField*>& fields,
		const LinearSolver& implicitSolver,
		LinearSolver::Debug* d = 0
	);
	
	// update particle deformation gradients based on grid velocities
	void updateDeformationGradients( float timeStep );
	
	// transfer grid velocities to particles:
	void updateParticleVelocities();

	typedef std::vector< std::pair<Sim::ConstIndexIterator, Sim::ConstIndexIterator> > PartitionList;
	const PartitionList& partition( int i, int j, int k ) const;
	
	void calculateExplicitMomenta(
		Eigen::VectorXf& explicitMomenta,
		std::vector<char>& nodeCollided,
		float timeStep,
		const ConstitutiveModel& constitutiveModel,
		const Sim::CollisionObjectSet& collisionObjects,
		const std::vector<const ForceField*>& fields
	);

	void collisionVelocities(
		Eigen::VectorXf& vc,
		const std::vector<const CollisionObject*>& collisionObjects,
		const std::vector<char>& nodeCollided
	) const;
	
	void calculateForces(
		Eigen::VectorXf& forces, 
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const ForceField*>& fields ) const;

	void calculateForceDifferentials(
		Eigen::VectorXf& df,
		const Eigen::VectorXf& dx,
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const ForceField*>& fields ) const;
	
	// iterator for looping over all the points whose shape function
	// p is in the support of:
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
	
	// grabs a shape function iterator
	ShapeFunctionIterator& shapeFunctionIterator() const;
	
	Eigen::VectorXf masses;
	Eigen::VectorXf velocities;
	Eigen::VectorXf prevVelocities;
	
	// class for splatting quantities like mass onto the grid in paralell:
	class GridSplatter
	{
		public:
			GridSplatter(
				const Grid& g,
				Eigen::VectorXf& result
			);
			
			virtual ~GridSplatter() {}
			
			void setPartition( int i, int j, int k );

			void operator()(const tbb::blocked_range<int> &r) const;

			template< class T >
			T* particleVariable( const std::string& name );
			
		protected:

			virtual void splat(
				Sim::ConstIndexIterator begin,
				Sim::ConstIndexIterator end,
				Eigen::VectorXf& result ) const = 0;

			const Grid& m_g;

		private:

			const PartitionList* m_partition;
			mutable Eigen::VectorXf& m_result;
			const void* m_args;

	};
	
	const Eigen::Vector3f& minCoord() const;
	const Eigen::Vector3f& maxCoord() const;
	const Eigen::Vector3i& n() const;
	float gridSize() const;
	
	Sim::MaterialPointDataMap& m_d;	
	Eigen::Vector3f m_frameVelocity;

private:

	template <class Splatter>
	void splat( Splatter& s ) const;
	
	const Sim::IndexList& m_particleInds;
	float m_gridSize;
	const ShapeFunction& m_shapeFunction;
	int m_dimension;
	
	Eigen::Vector3f m_min;
	Eigen::Vector3f m_max;
	Eigen::Vector3i m_n;
	
	void computeProcessingPartitions();
	PartitionList m_processingPartitions[2][2][2];

	friend class ImplicitUpdateMatrix;
	std::vector<char> m_nodeCollided;
	
	static inline void minMax( float x, float& min, float& max );
	inline int fixDim( float& min, float& max ) const;
	
	mutable tbb::enumerable_thread_specific< std::auto_ptr< ShapeFunctionIterator > > m_shapeFunctionIterators;
};

} // namespace MpmSim

#include "MpmSim/Grid.inl"

#endif // MPMSIM_GRID_H
