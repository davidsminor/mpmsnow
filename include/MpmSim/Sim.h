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

namespace MpmSimTest
{
class TestSimClass;
}

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
		ConstitutiveModel& model,
		const CollisionObjectSet& collisionObjects,
		const ForceFieldSet& forceFields,
		int dimension=3
	);
	
	// accessor for particle data:
	MaterialPointData& particleData();

	// complete a full simulation time step:
	void advance( float timeStep, TerminationCriterion& terminationCriterion, LinearSolver::Debug* d = 0 );
	
	typedef std::vector<int> IndexList;
	typedef IndexList::iterator IndexIterator;
	typedef IndexList::const_iterator ConstIndexIterator;
	
	typedef std::vector<IndexList> BodyList;
	typedef BodyList::iterator BodyIterator;
	typedef BodyList::const_iterator ConstBodyIterator;
	
private:
	
	// partition the sim into contiguous bodies:
	void calculateBodies();
	
	// material point data for all the particles
	MaterialPointData m_particleData;
	
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
	ConstitutiveModel& m_constitutiveModel;
	
	// collision objects:
	const CollisionObjectSet& m_collisionObjects;

	// force fields:
	const ForceFieldSet& m_forceFields;

	// dimension:
	int m_dimension;
	
	// testing:
	friend class MpmSimTest::TestSimClass;

};

} // namespace MpmSim

#endif //MPMSIM_SIM_H
