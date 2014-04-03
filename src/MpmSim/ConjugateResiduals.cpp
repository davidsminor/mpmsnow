#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/Grid.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

ConjugateResiduals::ConjugateResiduals( int iters, float tol_error, bool log ) :
	m_iters( iters ),
	m_tolError( tol_error ),
	m_log( log )
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
	
	float bNorm2 = b.squaredNorm();
	if(bNorm2 == 0) 
	{
		return;
	}
	float threshold = m_tolError*m_tolError*bNorm2;
	
	const int N = (int)b.size();
	
	VectorXf r(N);
	A.multVector( x, r );
	r = b - r;
	
	VectorXf Ar( N );
	A.multVector( r, Ar );
	float rAr = r.dot( Ar );
	
	VectorXf p(N);
	VectorXf Ap(N);
	for( int i=0; i < b.size(); ++i )
	{
		p[i] = r[i];
		Ap[i] = Ar[i];
	}
	
	if( m_log )
	{
		residuals.push_back( r );
		searchDirections.push_back( p );
	}
	
	int i = 0;
	while( i < m_iters )
	{
		// minimize along search direction:
		float alpha = rAr / Ap.dot( Ap );
		x += alpha * p;
		
		// update residual:
		r -= alpha * Ap;

		float rNorm2 = r.squaredNorm();
		std::cerr << i << ": " << Ar.norm() << "," << sqrt( rNorm2 ) << " / " << sqrt( threshold ) << "  " << alpha << std::endl;
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
		
		if( m_log )
		{
			residuals.push_back( r );
			searchDirections.push_back( p );
		}
		
		++i;
	}

}

