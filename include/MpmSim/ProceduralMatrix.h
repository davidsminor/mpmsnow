#ifndef MPMSIM_PROCEDURALMATRIX_H
#define MPMSIM_PROCEDURALMATRIX_H

#include <Eigen/Dense>
namespace MpmSim
{

class ProceduralMatrix
{

public:

	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;

};

}

#endif //MPMSIM_PROCEDURALMATRIX_H
