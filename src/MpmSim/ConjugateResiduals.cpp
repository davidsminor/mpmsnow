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
	// Looks like the preconditioning stuff was broken though, so I changed that.
	const size_t N = b.size();
	
	float bNorm2 = b.squaredNorm();
	std::cerr << "conjugate residuals... rhs squared norm: " << bNorm2 << std::endl;
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
	A.multVector( r_precond, Ar );

	// rz = <r_precond, Ar>
	float rz = r_precond.dot( Ar );
	
	if( m_log )
	{
		m_residuals.push_back( r );
		m_searchDirections.push_back( p );
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

		// alpha <- <r_precond,Ar>/<Ap,P^-1 Ap>
		float ApdPAp = Ap.dot(precond_Ap);
		if( ApdPAp == 0 )
		{
			std::cerr << "terminating solve due to potential divide by zero" << std::endl;
			return;
		}
		float alpha =  rz / ApdPAp;

		// x <- x + alpha * p
		x = x + alpha * p;
		A.subspaceProject( x );
		
		// debug output:
		if( d )
		{
			(*d)( x );
		}
		
		// r_precond <- r_precond - alpha * P^-1 Ap
		r_precond = r_precond - alpha * precond_Ap;
		if( m_preconditioner )
		{
			m_preconditioner->multVector( r_precond, r );
		}
		else
		{
			r = r_precond;
		}
		
		if( m_terminationCriterion( r, i ) )
		{
			return;
		}
		
		// Ar <- A*r_precond
		A.multVector(r_precond, Ar);
		
		float rz_old = rz;

		// rz = <r_precond^H, r_precond>
		rz = r_precond.dot( Ar );
		
		if( rz_old == 0 )
		{
			std::cerr << "terminating solve due to potential divide by zero" << std::endl;
			return;
		}
		
		// beta <- <r_{i+1},r_{i+1}>/<r_precond,r_precond>
		float beta = rz / rz_old;

		// p <- r_precond + beta*p
		p = r_precond + beta * p;

		// Ap <- Ar + beta*Ap
		Ap = Ar + beta * Ap;
		
		if( m_log )
		{
			m_residuals.push_back( r );
			m_searchDirections.push_back( p );
		}
	}
	
}

