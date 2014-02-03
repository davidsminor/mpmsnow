#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/Grid.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

ConjugateResiduals::ConjugateResiduals( int iters, float tol_error ) :
	m_iters( iters ),
	m_tolError( tol_error )
{
}

void ConjugateResiduals::operator()
(
		const ProceduralMatrix& mat,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x
) const
{
	typedef float RealScalar;

	RealScalar rhsNorm2 = rhs.squaredNorm();
	if(rhsNorm2 == 0) 
	{
		x.setZero();
		return;
	}
	
	RealScalar tol = m_tolError;
	int maxIters = m_iters;


	// r0 = b - A x0
	VectorXf r;
	mat.multVector( x, r );
	r = rhs - r;
	
	RealScalar threshold = tol*tol*rhsNorm2;
	RealScalar residualNorm2 = r.squaredNorm();
	if (residualNorm2 < threshold)
	{
		return;
	}

	// p0 = r0
	VectorXf p = r;
	
	VectorXf Ar;
	mat.multVector( r, Ar );

	VectorXf Ap = Ar;

	float rtAr = r.dot( Ar );
	float ApdotAp = Ap.squaredNorm();

	int i = 0;
	while(i < maxIters)
	{
		float alpha = rtAr / ApdotAp;
		x += alpha * p;
		r -= alpha * Ap;
		
		residualNorm2 = r.squaredNorm();
		if(residualNorm2 < threshold)
		{
			break;
		}
		
		mat.multVector( r, Ar );
		float rtArNew = r.dot( Ar );
		float beta = rtArNew/rtAr;
		rtAr = rtArNew;

		p = r + beta * p;
		Ap = Ar + beta * Ap;
		
		std::cerr << i << ":" << sqrt( residualNorm2 ) << "/" << sqrt( threshold ) << std::endl;

		++i;
	}

}