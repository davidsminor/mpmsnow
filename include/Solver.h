#ifndef MPMSNOW_SOLVER_H
#define MPMSNOW_SOLVER_H

#include "Grid.h"
#include "ParticleData.h"
#include "CollisionObject.h"

#include <Eigen/Dense>

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

#endif
