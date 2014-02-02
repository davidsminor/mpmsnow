#ifndef MPMSIM_SOLVER_H
#define MPMSIM_SOLVER_H

#include "Grid.h"
#include "ParticleData.h"
#include "CollisionObject.h"

#include <Eigen/Dense>
namespace MpmSim
{

class Grid;

class Solver
{
public:
	virtual ~Solver() {}
	
	virtual void operator()(
		const Grid* g,
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x ) const = 0;
	
};

} //namespace MpmSim

#endif
