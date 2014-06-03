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
	
	class Debug
	{
		public:
			virtual void operator()( Eigen::VectorXf& x ) = 0;
	};
	
	virtual void operator()(
		const ProceduralMatrix& mat,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x,
		Debug* d=0 ) const = 0;
	
};

} //namespace MpmSim

#endif
