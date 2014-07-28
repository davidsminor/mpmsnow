#include "MpmSim/ImplicitUpdateMatrix.h"

using namespace Eigen;
using namespace MpmSim;

ImplicitUpdateMatrix::ImplicitUpdateMatrix(
	const Sim::MaterialPointDataMap& d,
	const Grid& g,
	const ConstitutiveModel& constitutiveModel,
	const std::vector<const CollisionObject*>& collisionObjects,
	const std::vector<const ForceField*>& fields,
	float timeStep
) :
	m_d( d ),
	m_g(g),
	m_constitutiveModel( constitutiveModel ),
	m_collisionObjects( collisionObjects ),
	m_fields( fields ),
	m_timeStep( timeStep )
{
}

void ImplicitUpdateMatrix::subspaceProject( Eigen::VectorXf& toProject ) const
{
	for( int i=0; i < m_g.m_n[0]; ++i )
	{
		for( int j=0; j < m_g.m_n[1]; ++j )
		{
			for( int k=0; k < m_g.m_n[2]; ++k )
			{
				int idx = m_g.coordsToIndex( i, j, k );
				if( m_g.m_nodeCollided[idx] == -1 )
				{
					// no collision
					continue;
				}
				else if( m_g.m_nodeCollided[idx] == -2 )
				{
					// more than one collision: set to zero
					toProject.segment<3>( 3 * idx ).setZero();
				}
				else
				{
					const CollisionObject* obj = m_collisionObjects[ m_g.m_nodeCollided[idx] ];
					Vector3f v = toProject.segment<3>( 3 * idx );
					
					// find object normal:
					Vector3f x( m_g.m_gridSize * i + m_g.m_min[0], m_g.m_gridSize * j + m_g.m_min[1], m_g.m_gridSize * k + m_g.m_min[2] );
					Vector3f n;
					obj->grad( x, n );
					n.normalize();
					float nDotP = n.dot( v );
					
					// project out component perpendicular to the object
					v -= nDotP * n;
					toProject.segment<3>( 3 * idx ) = v;
				}
				
			}
		}
	}
}

void ImplicitUpdateMatrix::multVector( const Eigen::VectorXf& vNPlusOne, Eigen::VectorXf& result ) const
{
	// This method computes the forward momenta in this frame in terms of the velocities
	// in the next frame:
	// m * v^(n+1) - m_timeStep * dF(v^(n+1) * m_timeStep)
	
	// apply collisions to input:
	result = vNPlusOne;
	subspaceProject( result );

	// work out force differentials when you perturb the grid positions by vTransformed * m_timeStep:
	VectorXf df( result.size() );
	m_g.calculateForceDifferentials( df, result, m_constitutiveModel, m_fields );
	
	// convert result to a momentum, and subtract off df multiplied by the time step:
	for( int idx=0; idx < m_g.masses.size(); ++idx )
	{
		result.segment<3>( 3*idx ) = m_g.masses[ idx ] * result.segment<3>( 3 * idx ) - m_timeStep * m_timeStep * df.segment<3>( 3 * idx );
	}

	// apply collisions to output:
	subspaceProject( result );
}

