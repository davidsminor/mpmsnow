
#include "MpmSim/SquareMagnitudeTermination.h"

#include <iostream>

using namespace MpmSim;

SquareMagnitudeTermination::SquareMagnitudeTermination( int maxIters, float tolError )
	: m_maxIters( maxIters ), m_tolError( tolError )
{
}

void SquareMagnitudeTermination::init( const ProceduralMatrix& A, const Eigen::VectorXf& b )
{
	float bNorm2 = b.squaredNorm();
	m_threshold = m_tolError*m_tolError*bNorm2;
}

bool SquareMagnitudeTermination::operator()( Eigen::VectorXf& r, int iterationNum ) const
{
	if( iterationNum >= m_maxIters )
	{
		return true;
	}
	float rNorm2 = r.squaredNorm();
	std::cerr << iterationNum << ":" << sqrt( rNorm2 ) << " / " << sqrt( m_threshold ) << std::endl;
	return rNorm2 < m_threshold;
}

bool SquareMagnitudeTermination::cancelled() const
{
	return false;
}