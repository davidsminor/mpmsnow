#ifndef MPMSIM_RELATIVEMAGNITUDETERMINATION_H
#define MPMSIM_RELATIVEMAGNITUDETERMINATION_H

#include "MpmSim/TerminationCriterion.h"

namespace MpmSim
{

class SquareMagnitudeTermination : public TerminationCriterion
{

public:

	SquareMagnitudeTermination( int maxIters, float tolError );
	virtual void init( const ProceduralMatrix& A, const Eigen::VectorXf& b );
	virtual bool operator()( Eigen::VectorXf& r, int iterationNumber ) const;
	virtual bool cancelled() const;

private:

	float m_maxIters;
	float m_tolError;
	float m_threshold;

};

}; // namespace MpmSim

#endif // MPMSIM_RELATIVEMAGNITUDETERMINATION_H
