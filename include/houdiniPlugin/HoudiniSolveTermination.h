#ifndef MPMSIM_HOUDINISOLVETERMINATION_H
#define MPMSIM_HOUDINISOLVETERMINATION_H

#include "MpmSim/SquareMagnitudeTermination.h"
#include <UT/UT_Interrupt.h>

namespace MpmSim
{

class HoudiniSolveTermination : public SquareMagnitudeTermination
{

public:

	HoudiniSolveTermination( int maxIters, float tolError, UT_Interrupt* utInterrupt );
	virtual bool operator()( Eigen::VectorXf& r, int iterationNumber ) const;

private:
	
	UT_Interrupt* m_utInterrupt;

};

}; // namespace MpmSim

#endif // MPMSIM_HOUDINISOLVETERMINATION_H
