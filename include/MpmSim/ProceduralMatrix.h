#ifndef MPMSIM_PROCEDURALMATRIX_H
#define MPMSIM_PROCEDURALMATRIX_H

#include <Eigen/Dense>
namespace MpmSim
{

class ProceduralMatrix
{

public:
	
	virtual ~ProceduralMatrix() {}
	
	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;
	virtual void multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const = 0;
	virtual void subspaceProject( Eigen::VectorXf& x ) const = 0;

};

}

#endif //MPMSIM_PROCEDURALMATRIX_H
