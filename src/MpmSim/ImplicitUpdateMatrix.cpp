#include "MpmSim/ImplicitUpdateMatrix.h"

using namespace Eigen;
using namespace MpmSim;

ImplicitUpdateMatrix::ImplicitUpdateMatrix( const ParticleData& d, const Grid& g, const std::vector<CollisionObject*>& collisionObjects ) :
	m_d( d ),
	m_g(g),
	m_collisionObjects( collisionObjects )
{
}

void ImplicitUpdateMatrix::multVector( const Eigen::VectorXf& vNPlusOne, Eigen::VectorXf& result ) const
{
	// This method computes the forward momenta in this frame in terms of the velocities
	// in the next frame:
	// m * v^(n+1) - m_timeStep * dF(v^(n+1) * m_timeStep)
	
	// work out force differentials when you perturb the grid positions by v * m_timeStep:
	VectorXf df( vNPlusOne.size() );

	// so: this method effectively applies a symmetric matrix which is quite diagonally
	// dominant, with the diagonals largely controlled by the masses. Unfortunately, if
	// the masses vary significantly form cell to cell (which they almost always do),
	// this makes the eigenvalues SUCK, in that they range from something like 1.e-9 to 20.
	// The conjugate gradient solver hates this, so instead we transform the problem so that
	// those masses move off the main diagonal, but the matrix remains symmetric. This means
	// we need to divide both the input and the output of this function by the square roots
	// of the masses:
	
	// divide input:
	VectorXf vTransformed = vNPlusOne;
	for( int i=0; i < m_g.m_gridMasses.size(); ++i )
	{
		if( m_g.m_gridMasses[i] != 0 )
		{
			vTransformed.segment<3>( 3 * i ) /= sqrt( m_g.m_gridMasses[i] );
		}
	}

	for( int i=0; i < m_g.m_n[0]; ++i )
	{
		for( int j=0; j < m_g.m_n[1]; ++j )
		{
			for( int k=0; k < m_g.m_n[2]; ++k )
			{
				int idx = m_g.coordsToIndex( Vector3i( i, j, k ) );
				Vector3f v = vTransformed.segment<3>( 3 * idx );
				
				if( m_g.m_nodeCollided[idx] && m_g.m_gridMasses[ idx ] != 0 )
				{
					Vector3f x( m_g.m_gridH * i + m_g.m_min[0], m_g.m_gridH * j + m_g.m_min[1], m_g.m_gridH * k + m_g.m_min[2] );
					for( size_t objIdx = 0; objIdx < m_collisionObjects.size(); ++objIdx )
					{
						// intersecting the object
						Vector3f n;
						m_collisionObjects[objIdx]->grad( x, n );
						n.normalize();
						float nDotP = n.dot( v );
						
						// project out momentum perpendicular to the object
						v -= nDotP * n;
						
					}
				}
				
				vTransformed.segment<3>( 3 * idx ) = v;
			}
		}
	}
	
	m_g.calculateForceDifferentials( m_d, m_g.m_timeStep * vTransformed, df );
	
	result.resize( vTransformed.size() );
	for( int i=0; i < m_g.m_n[0]; ++i )
	{
		for( int j=0; j < m_g.m_n[1]; ++j )
		{
			for( int k=0; k < m_g.m_n[2]; ++k )
			{
				int idx = m_g.coordsToIndex( Vector3i( i, j, k ) );
				Vector3f resultMomentum = m_g.m_gridMasses[ idx ] * vTransformed.segment<3>( 3 * idx ) - m_g.m_timeStep * df.segment<3>( 3 * idx );
				
				// ok. So when you do this, is the matrix even symmetric any more? Maybe this should be in calculateForceDifferentials as well?
				if( m_g.m_nodeCollided[idx] )
				{
					Vector3f x( m_g.m_gridH * i + m_g.m_min[0], m_g.m_gridH * j + m_g.m_min[1], m_g.m_gridH * k + m_g.m_min[2] );
					for( size_t objIdx = 0; objIdx < m_collisionObjects.size(); ++objIdx )
					{
						// intersecting the object
						Vector3f n;
						m_collisionObjects[objIdx]->grad( x, n );
						n.normalize();
						float nDotP = n.dot( resultMomentum );
						
						// project out momentum perpendicular to the object
						resultMomentum -= nDotP * n;
						
					}
				}
				
				result.segment<3>( 3 * idx ) = resultMomentum;
			}
		}
	}
	
	// divide output:
	for( int i=0; i < m_g.m_gridMasses.size(); ++i )
	{
		if( m_g.m_gridMasses[i] != 0 )
		{
			result.segment<3>( 3 * i ) /= sqrt( m_g.m_gridMasses[i] );
		}
	}
}

