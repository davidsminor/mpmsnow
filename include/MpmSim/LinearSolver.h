#ifndef MPMSIM_LINEARSOLVER_H
#define MPMSIM_LINEARSOLVER_H

#include "ParticleData.h"
#include "CollisionObject.h"

#include <Eigen/Dense>
namespace MpmSim
{

class Grid;

class LinearSolver
{
public:
	virtual ~LinearSolver() {}
	
	virtual void operator()(
		const Grid* g,
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x ) const = 0;
	
};

} //namespace MpmSim

#endif
