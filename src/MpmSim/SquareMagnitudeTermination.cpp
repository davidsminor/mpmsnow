
#include "MpmSim/SquareMagnitudeTermination.h"

#include <iostream>

using namespace MpmSim;

SquareMagnitudeTermination::SquareMagnitudeTermination( float tolError ) : m_tolError( tolError )
{
}

void SquareMagnitudeTermination::init( const ProceduralMatrix& A, const Eigen::VectorXf& b )
{
	float bNorm2 = b.squaredNorm();
	m_threshold = m_tolError*m_tolError*bNorm2;
}

bool SquareMagnitudeTermination::operator()( Eigen::VectorXf& r ) const
{
	float rNorm2 = r.squaredNorm();
	std::cerr << sqrt( rNorm2 ) << " / " << sqrt( m_threshold ) << std::endl;
	return rNorm2 < m_threshold;
}