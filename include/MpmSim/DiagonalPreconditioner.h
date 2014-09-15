#ifndef MPMSIM_DIAGONALPRECONDITIONER_H
#define MPMSIM_DIAGONALPRECONDITIONER_H

#include "MpmSim/ProceduralMatrix.h"
#include "MpmSim/Grid.h"

#include <Eigen/Dense>
namespace MpmSim
{

class DiagonalPreconditioner : public ProceduralMatrix
{

public:

	DiagonalPreconditioner(
		const Grid& g,
		const ConstitutiveModel& constitutiveModel,
		float timeStep
	);
	
	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	virtual void multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	
	void subspaceProject( Eigen::VectorXf& x ) const;
	
private:
	Eigen::VectorXf m_implicitUpdateDiagonal;
};

}

#endif //MPMSIM_DIAGONALPRECONDITIONER_H
