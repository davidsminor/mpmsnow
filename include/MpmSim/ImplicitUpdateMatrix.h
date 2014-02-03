#ifndef MPMSIM_IMPLICITUPDATEMATRIX_H
#define MPMSIM_IMPLICITUPDATEMATRIX_H

#include "MpmSim/ProceduralMatrix.h"
#include "MpmSim/ParticleData.h"
#include "MpmSim/Grid.h"
#include "MpmSim/CollisionObject.h"

#include <Eigen/Dense>
namespace MpmSim
{

class ImplicitUpdateMatrix : public ProceduralMatrix
{

public:

	ImplicitUpdateMatrix( const ParticleData& d, const Grid& g, const std::vector<CollisionObject*>& collisionObjects );
	
	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const;
	
private:
	
	const ParticleData& m_d;
	const Grid& m_g;
	const std::vector< CollisionObject* >& m_collisionObjects;

};

}

#endif //MPMSIM_PROCEDURALMATRIX_H
