#ifndef MPMSNOW_CONJUGATERESIDUALS_H
#define MPMSNOW_CONJUGATERESIDUALS_H

#include "Solver.h"

class ConjugateResiduals : public Solver
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

#endif
