#include "MpmSim/ConjugateResiduals.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

ConjugateResiduals::ConjugateResiduals(
	int iters,
	TerminationCriterion& terminationCriterion,
	const ProceduralMatrix* preconditioner,
	bool log
) :
	m_iters( iters ),
	m_terminationCriterion( terminationCriterion ),
	m_preconditioner( preconditioner ),
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
	// I copy pasted a lot of this from here:
	// https://github.com/cusplibrary/cusplibrary/commit/dc757ce17229243d3852994f3c6c9c789936620a
	
	const size_t N = b.size();
	const size_t recompute_r = 8;	     // interval to update r
	
	float bNorm2 = b.squaredNorm();
	if(bNorm2 == 0) 
	{
		x.setZero();
		return;
	}
	m_terminationCriterion.init( A, b );

	// allocate workspace
	Eigen::VectorXf y(N);
	Eigen::VectorXf z(N);
	Eigen::VectorXf r(N);
	Eigen::VectorXf p(N);
	Eigen::VectorXf Az(N);
	Eigen::VectorXf Ax(N);

	// y <- A*x
	A.multVector(x, Ax);

	// r <- b - A*x
	r = b - Ax;

	// z <- M^-1*r
	if( m_preconditioner )
	{
		m_preconditioner->multInverseVector( r, z );
	}
	else
	{
		z = r;
	}

	// p <- z
	p = z;

	// y <- A*p
	A.multVector( p, y );

	// Az <- A*z
	A.multVector( z, Az );

	// rz = <r^H, z>
	float rz = r.dot( Az );
	
	if( m_log )
	{
		residuals.push_back( r );
		searchDirections.push_back( p );
	}

	for( int i=0; i < m_iters; ++i )
	{
		// alpha <- <r,z>/<y,p>
		float alpha =  rz / y.dot(y);

		// x <- x + alpha * p
		x = x + alpha * p;
		A.subspaceProject( x );
		
		// debug output:
		if( d )
		{
			(*d)( x );
		}
		
		if( (i % recompute_r) && (i > 0) )
		{
			// r <- r - alpha * y
			r = r - alpha * y;
		}
		else
		{
			// y <- A*x
			A.multVector(x, Ax);

			// r <- b - A*x
			r = b - Ax;
		}
		
		if( m_terminationCriterion( r ) )
		{
			return;
		}
		
		// z <- M*r
		if( m_preconditioner )
		{
			m_preconditioner->multInverseVector( r, z );
		}
		else
		{
			z = r;
		}
		
		// Az <- A*z
		A.multVector(z, Az);
		
		float rz_old = rz;

		// rz = <r^H, z>
		rz = r.dot( Az );

		// beta <- <r_{i+1},r_{i+1}>/<r,r>
		float beta = rz / rz_old;

		// p <- r + beta*p
		p = r + beta * p;

		// y <- Az + beta*y
		y = Az + beta * y;
		
		if( m_log )
		{
			residuals.push_back( r );
			searchDirections.push_back( p );
		}
	}

	


	/*
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
	
	VectorXf r( N );
	VectorXf p( N );

	if( m_preconditioner )
	{
		m_preconditioner->multInverseVector( b, r );
		p = r;
	}
	else
	{
		r = b;
		p = b;
	}
	
	if( m_log )
	{
		residuals.push_back( r );
		searchDirections.push_back( p );
	}

	VectorXf Ap( N );
	A.multVector( p, Ap );

	VectorXf s;
	
	if( m_preconditioner )
	{
		s.resize( N );
		m_preconditioner->multInverseVector( Ap, s );
	}

	VectorXf Ar( N );
	A.multVector( r, Ar );

	float rAr = r.dot( Ar );

	int i = 0;
	while( i < m_iters )
	{
		// minimize along search direction:
		float alpha;
		if( m_preconditioner )
		{
			alpha = rAr / Ap.dot( s );
		}
		else
		{
			alpha = rAr / Ap.dot( Ap );
		}

		x += alpha * p;
		
		if( d )
		{
			(*d)( x );
		}

		// update residual:
		if( m_preconditioner )
		{
			r -= alpha * s;
		}
		else
		{
			r -= alpha * Ap;
			A.subspaceProject( r );
		}

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
	*/
}

