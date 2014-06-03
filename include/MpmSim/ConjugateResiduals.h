#ifndef MPMSIM_CONJUGATERESIDUALS_H
#define MPMSIM_CONJUGATERESIDUALS_H

#include "LinearSolver.h"

namespace MpmSim
{

class ConjugateResiduals : public LinearSolver
{
public:
	
	ConjugateResiduals( int iters, float tol_error, bool log=false );

	virtual void operator()(
		const ProceduralMatrix& mat,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x,
		Debug* d=0 ) const;

	mutable std::vector<Eigen::VectorXf> residuals;
	mutable std::vector<Eigen::VectorXf> searchDirections;

private:

	int m_iters;
	float m_tolError;
	bool m_log;

};

} // namespace MpmSim

#endif
