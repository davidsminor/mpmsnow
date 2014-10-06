
#include "houdiniPlugin/HoudiniSolveTermination.h"

#include <iostream>

using namespace MpmSim;

HoudiniSolveTermination::HoudiniSolveTermination( int maxIters, float tolError, UT_Interrupt* utInterrupt )
	: SquareMagnitudeTermination( maxIters, tolError ), m_utInterrupt( utInterrupt )
{
}

bool HoudiniSolveTermination::cancelled() const
{
	return m_utInterrupt->opInterrupt();
}

bool HoudiniSolveTermination::operator()( Eigen::VectorXf& r, int iterationNum ) const
{
	if( m_utInterrupt->opInterrupt() )
	{
		return true;
	}
	
	return SquareMagnitudeTermination::operator()( r, iterationNum );
}
