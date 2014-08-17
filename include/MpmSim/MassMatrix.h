#ifndef MPMSIM_MASSMATRIX_H
#define MPMSIM_MASSMATRIX_H

#include "MpmSim/ProceduralMatrix.h"
#include "MpmSim/Grid.h"

#include <Eigen/Dense>
namespace MpmSim
{

class MassMatrix : public ProceduralMatrix
{

public:

	MassMatrix(
		const Grid& g
	);
	
	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	virtual void multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	
	void subspaceProject( Eigen::VectorXf& x ) const;
	
private:
	const Grid& m_g;
};

}

#endif //MPMSIM_MASSMATRIX_H
