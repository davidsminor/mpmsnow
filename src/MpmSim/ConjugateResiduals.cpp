#include "MpmSim/ConjugateResiduals.h"
#include "MpmSim/Grid.h"

#include <GL/gl.h>
#include <GL/glut.h>

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
		Eigen::VectorXf& x,
		Debug* d
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
	
	VectorXf r(N);
	VectorXf p(N);
	for( int i=0; i < b.size(); ++i )
	{
		r[i] = b[i];
		p[i] = b[i];
	}
	
	if( m_log )
	{
		residuals.push_back( r );
		searchDirections.push_back( p );
	}

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
		
		if( d )
		{
			(*d)( x );
		}

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

