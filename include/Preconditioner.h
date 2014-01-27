#ifndef MPMSNOW_PRECONDITIONER_H
#define MPMSNOW_PRECONDITIONER_H

#include <Eigen/Dense>

class Preconditioner
{

public:

	virtual void apply( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;

};

#endif // MPMSNOW_PRECONDITIONER_H
