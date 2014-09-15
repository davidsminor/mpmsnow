#ifndef MPMSIM_SIM_H
#define MPMSIM_SIM_H

#include "Eigen/Dense"

#include "MpmSim/MaterialPointData.h"
#include "MpmSim/ShapeFunction.h"
#include "MpmSim/TerminationCriterion.h"
#include "MpmSim/LinearSolver.h"

#include <vector>
#include <map>
#include <memory>

namespace MpmSim
{

class ConstitutiveModel;
class CollisionObject;
class ForceField;

class Sim
{

public:

	class CollisionObjectSet
	{
	public:
		~CollisionObjectSet();
		void add( CollisionObject* o );
		std::vector<const CollisionObject*> objects;
		int collide(
			Eigen::Vector3f& v,
			const Eigen::Vector3f& x,
			const Eigen::Vector3f& frameVelocity,
			bool addCollisionVelocity = false
		) const;
	};

	class ForceFieldSet
	{
	public:
		~ForceFieldSet();
		void add( ForceField* f );
		std::vector<const ForceField*> fields;
	};

	// construct a sim from initial conditions:
	Sim(
		const std::vector<Eigen::Vector3f>& x,
		const std::vector<float>& masses,
		float gridSize,
		const ShapeFunction& shapeFunction,
		const ConstitutiveModel& model,
		const CollisionObjectSet& collisionObjects,
		const ForceFieldSet& forceFields,
		int dimension=3
	);
	
	~Sim();
	
	// complete a full simulation time step:
	void advance( float timeStep, TerminationCriterion& terminationCriterion, LinearSolver::Debug* d = 0 );

	typedef std::vector<int> IndexList;
	typedef IndexList::iterator IndexIterator;
	typedef IndexList::const_iterator ConstIndexIterator;

	typedef std::vector<IndexList> BodyList;
	typedef BodyList::iterator BodyIterator;
	typedef BodyList::const_iterator ConstBodyIterator;

	typedef std::map< std::string, MaterialPointDataBase* > MaterialPointDataMap;
	
	// SIM PARTITIONS:

	// number of contiguous bodies:
	size_t numBodies() const;
	
	// particle indices in body n:
	const IndexList& body( size_t n ) const;
	
	// particle indices for the ballistic particles:
	const IndexList& ballisticParticles() const;
	
	// PARTICLE VARIABLES:

	// number of variables per particle:
	size_t numParticleVariables() const;
	
	// material point data for all the particles
	MaterialPointDataMap particleData;

	// convenience methods for accessing per particle data by name:
	template< class T >
	T* particleVariable( const std::string& name );
	
	template< class T >
	const T* particleVariable( const std::string& name ) const;

	// class for doing neighbour queries on particles:
	class NeighbourQuery
	{
		public:
			NeighbourQuery( const std::vector<Eigen::Vector3f>& particleX, float r );
			void neighbours( const Eigen::Vector3f& p, std::vector<int>& neighbourInds ) const;
		private:

			void neighboursInCell( const Eigen::Vector3f& p, const Eigen::Vector3i& cell, std::vector<int>& neighbourInds ) const;

			size_t voxelOffset( const Eigen::Vector3i& cell ) const;
			
			Eigen::Vector3i m_cellMin;
			Eigen::Vector3i m_gridDims;

			const std::vector<Eigen::Vector3f>& m_particleX;
			float m_r;
			IndexList m_spatialSorting;
			std::vector<int> m_voxels;
	};
		
	// sort the specified index range into voxels:
	static void voxelSort(
		IndexIterator begin,
		IndexIterator end,
		float voxelSize,
		const std::vector<Eigen::Vector3f>& particleX,
		int dim=0 );


private:
	
	// partition the sim into contiguous bodies:
	void calculateBodies();
	
	// a list of particles with no neighbours, which are treated in isolation:
	IndexList m_ballisticParticles;

	// the rest of the particles are partitioned into contiguous bodies, whose
	// dynamics are calculated using the material point method. The indices
	// are spatially sorted so contiguous particles are in the same voxel:
	std::vector< IndexList > m_bodies;
	
	// computational grid cell size
	float m_gridSize;
	
	// shape function for interpolating particles back and forth between grids
	const ShapeFunction& m_shapeFunction;

	// constitutive model governing the material physics
	const ConstitutiveModel& m_constitutiveModel;
	
	// collision objects:
	const CollisionObjectSet& m_collisionObjects;

	// force fields:
	const ForceFieldSet& m_forceFields;

	// dimension:
	int m_dimension;
};

} // namespace MpmSim

#include "MpmSim/Sim.inl"

#endif //MPMSIM_PARTICLEDATA_H
