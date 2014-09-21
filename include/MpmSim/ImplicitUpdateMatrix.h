#ifndef MPMSIM_IMPLICITUPDATEMATRIX_H
#define MPMSIM_IMPLICITUPDATEMATRIX_H

#include "MpmSim/ProceduralMatrix.h"
#include "MpmSim/Sim.h"
#include "MpmSim/Grid.h"
#include "MpmSim/CollisionObject.h"

#include <Eigen/Dense>
namespace MpmSim
{

class ImplicitUpdateMatrix : public ProceduralMatrix
{

public:

	ImplicitUpdateMatrix(
		const MaterialPointData& d,
		const Grid& g,
		const ConstitutiveModel& constitutiveModel,
		const std::vector<const CollisionObject*>& collisionObjects,
		const std::vector<const ForceField*>& fields,
		float timeStep
	);
	
	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	virtual void multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	
	void subspaceProject( Eigen::VectorXf& x ) const;
	
private:
	
	const MaterialPointData& m_d;
	const Grid& m_g;
	const ConstitutiveModel& m_constitutiveModel;
	const std::vector< const CollisionObject* >& m_collisionObjects;
	const std::vector<const ForceField*>& m_fields;
	float m_timeStep;

	typedef std::vector< const CollisionObject* >::const_iterator CollisionIterator;
	typedef std::vector< const CollisionObject* >::const_reverse_iterator ReverseCollisionIterator;

};

}

#endif //MPMSIM_IMPLICITUPDATEMATRIX_H
