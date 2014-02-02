#ifndef MPMSIM_CONJUGATERESIDUALS_H
#define MPMSIM_CONJUGATERESIDUALS_H

#include "LinearSolver.h"

namespace MpmSim
{

class ConjugateResiduals : public LinearSolver
{
public:
	
	ConjugateResiduals( int iters, float tol_error );

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
