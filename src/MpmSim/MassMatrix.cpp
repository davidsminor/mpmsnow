#include "MpmSim/MassMatrix.h"

using namespace MpmSim;

MassMatrix::MassMatrix( const Grid& g ) : m_g( g )
{
}

void MassMatrix::subspaceProject( Eigen::VectorXf& toProject ) const
{
	// not implemented
}

void MassMatrix::multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
{
	for( int i=0; i < m_g.masses.size(); ++i )
	{
		if( m_g.masses[i] != 0 )
		{
			result.segment<3>( 3 * i ) = x.segment<3>( 3 * i ) / m_g.masses[i];
		}
		else
		{
			result.segment<3>( 3 * i ).setZero();
		}
	}
}

void MassMatrix::multVector( const Eigen::VectorXf& vNPlusOne, Eigen::VectorXf& result ) const
{
	// not implemented
}
