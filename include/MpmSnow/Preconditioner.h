#ifndef MPMSNOW_PRECONDITIONER_H
#define MPMSNOW_PRECONDITIONER_H

#include <Eigen/Dense>

namespace MpmSnow
{

class Preconditioner
{

public:

	virtual void apply( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;

};

} // namespace MpmSnow

#endif // MPMSNOW_PRECONDITIONER_H
