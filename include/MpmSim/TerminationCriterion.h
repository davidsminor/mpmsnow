#ifndef MPMSIM_TERMINATIONCRITERION_H
#define MPMSIM_TERMINATIONCRITERION_H

#include "MpmSim/ProceduralMatrix.h"
#include <Eigen/Dense>

namespace MpmSim
{

class TerminationCriterion
{
	public:
		virtual void init( const ProceduralMatrix& A, const Eigen::VectorXf& b ) = 0;
		virtual bool operator()( Eigen::VectorXf& r ) const = 0;
};

}; // namespace MpmSim

#endif // MPMSIM_TERMINATIONCRITERION_H
