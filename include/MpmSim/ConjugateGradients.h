#ifndef MPMSIM_CONJUGATEGRADIENTS_H
#define MPMSIM_CONJUGATEGRADIENTS_H

#include "LinearSolver.h"

namespace MpmSim
{

class ConjugateGradients : public LinearSolver
{
public:
	
	ConjugateGradients( int iters, float tol_error );

	virtual void operator()(
		const Grid* g,
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x ) const;

private:

	int m_iters;
	float m_tolError;

};

} // namespace MpmSim

#endif
