#include "MpmSim/ConjugateResiduals.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

ConjugateResiduals::ConjugateResiduals(
	TerminationCriterion& terminationCriterion,
	const ProceduralMatrix* preconditioner,
	bool log
) :
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
	Eigen::VectorXf r(N);
	Eigen::VectorXf r_precond(N);
	Eigen::VectorXf Ar(N);
	Eigen::VectorXf p(N);
	Eigen::VectorXf Ap(N);
	Eigen::VectorXf precond_Ap(N);
	Eigen::VectorXf Ax(N);

	// Ax <- A*x
	A.multVector(x, Ax);

	// r <- b - A*x
	r = b - Ax;

	// r_precond <- M^-1*r
	if( m_preconditioner )
	{
		m_preconditioner->multInverseVector( r, r_precond );
	}
	else
	{
		r_precond = r;
	}

	// p <- r_precond
	p = r_precond;

	// Ap <- A*p
	A.multVector( p, Ap );

	// Ar <- A*r
	r = r_precond;
	A.multVector( r, Ar );

	// rz = <r, Ar>
	float rz = r.dot( Ar );
	
	if( m_log )
	{
		residuals.push_back( r );
		searchDirections.push_back( p );
	}

	for( int i=0; ; ++i )
	{
		if( m_preconditioner )
		{
			m_preconditioner->multInverseVector( Ap, precond_Ap );
		}
		else
		{
			precond_Ap = Ap;
		}

		// alpha <- <r,Ar>/<Ap,P^-1 Ap>
		float alpha =  rz / Ap.dot(precond_Ap);

		// x <- x + alpha * p
		x = x + alpha * p;
		A.subspaceProject( x );
		
		// debug output:
		if( d )
		{
			(*d)( x );
		}
		
		// r <- r - alpha * P^-1 Ap
		r = r - alpha * precond_Ap;
		
		if( m_terminationCriterion( r, i ) )
		{
			return;
		}
		
		// Ar <- A*r
		A.multVector(r, Ar);
		
		float rz_old = rz;

		// rz = <r^H, r>
		rz = r.dot( Ar );

		// beta <- <r_{i+1},r_{i+1}>/<r,r>
		float beta = rz / rz_old;

		// p <- r + beta*p
		p = r + beta * p;

		// Ap <- Ar + beta*Ap
		Ap = Ar + beta * Ap;
		
		if( m_log )
		{
			residuals.push_back( r );
			searchDirections.push_back( p );
		}
	}
	
}

