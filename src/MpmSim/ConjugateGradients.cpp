#include "MpmSim/ConjugateGradients.h"
#include "MpmSim/Grid.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

ConjugateGradients::ConjugateGradients( int iters, float tol_error ) :
	m_iters( iters ),
	m_tolError( tol_error )
{
}

void ConjugateGradients::operator()
(
		const Grid* g,
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x
) const
{
	using std::sqrt;
	using std::abs;
	typedef float RealScalar;
	typedef float Scalar;

	RealScalar tol = m_tolError;
	int maxIters = m_iters;

	int n = (int)rhs.size();

	VectorXf residual;
	g->applyImplicitUpdateMatrix( d, collisionObjects, x, residual );
	residual = rhs - residual;

	RealScalar rhsNorm2 = rhs.squaredNorm();
	if(rhsNorm2 == 0) 
	{
		x.setZero();
		return;
	}
	RealScalar threshold = tol*tol*rhsNorm2;
	RealScalar residualNorm2 = residual.squaredNorm();
	if (residualNorm2 < threshold)
	{
		//iters = 0;
		//tol_error = sqrt(residualNorm2 / rhsNorm2);
		return;
	}

	VectorXf p = residual;		//initial search direction
	
	VectorXf tmp(n);
	RealScalar absNew = numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
	int i = 0;
	while(i < maxIters)
	{
		g->applyImplicitUpdateMatrix( d, collisionObjects, p, tmp ); // the bottleneck of the algorithm
		
		Scalar alpha = absNew / p.dot(tmp);   // the amount we travel on dir
		x += alpha * p;                       // update solution
		residual -= alpha * tmp;              // update residue

		residualNorm2 = residual.squaredNorm();
		if(residualNorm2 < threshold)
		  break;

		RealScalar absOld = absNew;
		absNew = numext::real(residual.squaredNorm());     // update the absolute value of r
		RealScalar beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
		p = residual + beta * p;                             // update search direction
		
		std::cerr << i << ":" << sqrt( residualNorm2 ) << "/" << sqrt( threshold ) << std::endl;

		i++;
	}
	//tol_error = sqrt(residualNorm2 / rhsNorm2);
}