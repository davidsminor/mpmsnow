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
		const ProceduralMatrix& A,
		const Eigen::VectorXf& b,
		Eigen::VectorXf& x
) const
{
	// NB: there are papers with an extra bit that supposedly makes
	// this more numerically stable
	// (eg http://webmail.cs.yale.edu/publications/techreports/tr107.pdf).
	// Do we need this?
	x.setZero();

	float bNorm2 = b.squaredNorm();
	if(bNorm2 == 0) 
	{
		return;
	}
	float threshold = m_tolError*m_tolError*bNorm2;
	
	const int N = (int)b.size();
	
	VectorXf r = b;
	
	VectorXf p = r;
	
	VectorXf Ap( N );
	A.multVector( p, Ap );

	VectorXf Ar( N );
	A.multVector( r, Ar );

	float rAr = r.dot( Ar );

	int i = 0;
	while( i < m_iters )
	{
		// minimize along search direction:
		float alpha = rAr / Ap.dot( Ap );
		x += alpha * p;
		
		// update residual:
		r -= alpha * Ap;

		float rNorm2 = r.squaredNorm();
		std::cerr << i << ": " << sqrt( rNorm2 ) << " / " << sqrt( threshold ) << std::endl;
		if( rNorm2 < threshold )
		{
			return;
		}
		
		// find a new search direction that's A^2 orto to the previous ones:
		A.multVector( r, Ar );
		float rArOld = rAr;
		
		rAr = r.dot( Ar );
		float beta = rAr / rArOld;
		
		Ap = Ar + beta * Ap;
		p = r + beta * p;
		++i;
	}

}