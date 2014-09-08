#ifndef MPMSIM_RELATIVEMAGNITUDETERMINATION_H
#define MPMSIM_RELATIVEMAGNITUDETERMINATION_H

#include "MpmSim/TerminationCriterion.h"

namespace MpmSim
{

class SquareMagnitudeTermination : public TerminationCriterion
{

public:

	SquareMagnitudeTermination( float tolError );
	virtual void init( const ProceduralMatrix& A, const Eigen::VectorXf& b );
	virtual bool operator()( Eigen::VectorXf& r ) const;

private:

	float m_tolError;
	float m_threshold;

};

}; // namespace MpmSim

#endif // MPMSIM_RELATIVEMAGNITUDETERMINATION_H
