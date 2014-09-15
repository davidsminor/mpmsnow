#include "MpmSim/DiagonalPreconditioner.h"

using namespace MpmSim;

DiagonalPreconditioner::DiagonalPreconditioner(
	const Grid& g,
	const ConstitutiveModel& constitutiveModel,
	float timeStep )
{
	g.dForceidXi( m_implicitUpdateDiagonal, constitutiveModel );
	m_implicitUpdateDiagonal *= - timeStep * timeStep;
	for( int i=0; i < g.masses.size(); ++i )
	{
		m_implicitUpdateDiagonal[3 * i + 0] += g.masses[i];
		m_implicitUpdateDiagonal[3 * i + 1] += g.masses[i];
		m_implicitUpdateDiagonal[3 * i + 2] += g.masses[i];
		if( m_implicitUpdateDiagonal[3 * i + 0] == 0 )
		{
			m_implicitUpdateDiagonal[3 * i + 0] = 1;
		}
		if( m_implicitUpdateDiagonal[3 * i + 1] == 0 )
		{
			m_implicitUpdateDiagonal[3 * i + 1] = 1;
		}
		if( m_implicitUpdateDiagonal[3 * i + 2] == 0 )
		{
			m_implicitUpdateDiagonal[3 * i + 2] = 1;
		}
	}
}

void DiagonalPreconditioner::subspaceProject( Eigen::VectorXf& toProject ) const
{
	// not implemented
}

void DiagonalPreconditioner::multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
{
	result = x.cwiseQuotient( m_implicitUpdateDiagonal );
}

void DiagonalPreconditioner::multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
{
	result = x.cwiseProduct( m_implicitUpdateDiagonal );
}
