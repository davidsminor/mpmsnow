#ifndef MPMSIM_PRECONDITIONER_H
#define MPMSIM_PRECONDITIONER_H

#include <Eigen/Dense>

namespace MpmSim
{

class Preconditioner
{

public:

	virtual void apply( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;

};

} // namespace MpmSim

#endif // MPMSIM_PRECONDITIONER_H
