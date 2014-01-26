#include "ConjugateResiduals.h"

#include <iostream>

using namespace Eigen;

ConjugateResiduals::ConjugateResiduals( int iters, float tol_error ) :
	m_iters( iters ),
	m_tolError( tol_error )
{
}

void ConjugateResiduals::operator()
(
		const Grid* g,
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x
) const
{
	// not implemented
}