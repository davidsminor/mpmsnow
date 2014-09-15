#ifndef MPMSIM_CONJUGATERESIDUALS_H
#define MPMSIM_CONJUGATERESIDUALS_H

#include "LinearSolver.h"
#include "TerminationCriterion.h"
#include <vector>

namespace MpmSim
{

class ConjugateResiduals : public LinearSolver
{
public:
	
	ConjugateResiduals( TerminationCriterion& terminationCriterion, const ProceduralMatrix* preconditioner = 0, bool log=false );

	virtual void operator()(
		const ProceduralMatrix& mat,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x,
		Debug* d=0 ) const;

	mutable std::vector<Eigen::VectorXf> residuals;
	mutable std::vector<Eigen::VectorXf> searchDirections;

private:
	
	TerminationCriterion& m_terminationCriterion;
	const ProceduralMatrix* m_preconditioner;
	bool m_log;

};

} // namespace MpmSim

#endif
