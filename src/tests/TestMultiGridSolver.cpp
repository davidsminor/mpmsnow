#include "tests/TestMultiGridSolver.h"
#include <iostream>
#include <Eigen/Dense>

namespace MpmSimTest
{

Eigen::VectorXf laplacian( const Eigen::VectorXf& x )
{
	Eigen::VectorXf l(x.size());
	for( int i=0; i < x.size(); ++i )
	{
		l[i] = 2 * x[i];
		if( i > 0 )
		{
			l[i] -= x[i-1];
		}
		if( i < x.size() - 1 )
		{
			l[i] -= x[i+1];
		}
	}
	return l;
}

void smooth( Eigen::VectorXf& x, const Eigen::VectorXf& b )
{
	x[0] = (1.0f/3)*( x[0] + x[1] + b[0] );
	for( int i=1; i < x.size()-1; ++i )
	{
		x[i] = (1.0f/3)*( x[i-1] + x[i] + x[i+1] + b[i] );
	}
	x[x.size()-1] = (1.0f/3)*( x[x.size()-1] + x[x.size()-2] + b[0] );
}

void restrict( Eigen::VectorXf& r )
{
	if( r.size() == 3 )
	{
		Eigen::VectorXf rNew( 1 );
		rNew[0] = 0.5f * r[1] + 0.25f * ( r[0] + r[2] );
		r = rNew;
		return;
	}

	int newSize = ( r.size() - 1 ) / 2 + 1;
	Eigen::VectorXf rNew( newSize );
	for( size_t i = 0; i < newSize; ++i )
	{
		rNew[i] = 0.5 * r[2*i];
		if( i > 0 )
		{
			rNew[i] += 0.25 * r[2*i-1];
		}
		
		if( i < newSize - 1 )
		{
			rNew[i] += 0.25 * r[2*i+1];
		}
	}
	r = rNew;
}

void interpolate( Eigen::VectorXf& r )
{
	if( r.size() == 1 )
	{
		Eigen::VectorXf rNew( 3 );
		rNew[0] = 0.5f * r[0];
		rNew[1] = r[0];
		rNew[2] = 0.5f * r[0];
		r = rNew;
		return;
	}

	int newSize = r.size() * 2 - 1;
	Eigen::VectorXf rNew( newSize );
	for( size_t i = 0; i < r.size(); ++i )
	{
		rNew[2*i] = r[i];
		if( i < r.size() - 1 )
		{
			rNew[2*i+1] = 0.5f * ( r[i] + r[i+1] );
		}
	}
	r = rNew;
}

void vCycle( Eigen::VectorXf& x, const Eigen::VectorXf& b )
{
	if( x.size() == 1 )
	{
		// solve exactly:
		x[0] = b[0]/2;
	}
	else
	{
		smooth(x,b);

		// solve for correction 'd' on a course grid:
		Eigen::VectorXf r = laplacian( x ) - b;
		
		restrict( r );
		Eigen::VectorXf d = Eigen::VectorXf::Zero( r.size() );
		vCycle( d, r );
		interpolate( d );
		x = x - d;
		smooth(x,b);
	}
}

void testMultiGridSolver()
{
	Eigen::VectorXf x = Eigen::VectorXf::Zero( 65 );
	Eigen::VectorXf b( 65 );
	
	for( size_t i=0; i < 65; ++i )
	{
		float f = float(i)/64;
		b[i] = f * (1-f);
	}

	// solve poisson via jacobi:
	for( size_t n = 0; n < 1000; ++n )
	{
		// smooth:
		smooth( x, b );
	}

	std::cerr << ( laplacian( x ) - b ).norm() / x.size() << std::endl;
	std::cerr << std::endl;

	// solve poisson equation via multigrid:
	x = Eigen::VectorXf::Zero( 65 );
	for( size_t n = 0; n < 1000; ++n )
	{
		// smooth:
		vCycle( x, b );
	}

	std::cerr << ( laplacian( x ) - b ).norm() / x.size() << std::endl;
	std::cerr << std::endl;
}

}
