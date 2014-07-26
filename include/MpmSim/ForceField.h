#ifndef MPMSIM_FORCEFIELD_H
#define MPMSIM_FORCEFIELD_H

#include "Eigen/Dense"

namespace MpmSim
{

class ForceField
{
public:

	virtual ~ForceField() {}

	virtual Eigen::Vector3f force( const Eigen::Vector3f& x, float m ) const = 0;
	virtual Eigen::Matrix3f dFdx( const Eigen::Vector3f& x, float m ) const = 0;
};

} //namespace MpmSim

#endif
