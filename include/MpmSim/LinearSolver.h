#ifndef MPMSIM_LINEARSOLVER_H
#define MPMSIM_LINEARSOLVER_H

#include "ProceduralMatrix.h"
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
		const ProceduralMatrix& mat,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x ) const = 0;
	
};

} //namespace MpmSim

#endif
