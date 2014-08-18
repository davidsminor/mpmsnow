#ifndef MPMSIM_GRAVITYFIELD_H
#define MPMSIM_GRAVITYFIELD_H

#include "MpmSim/ForceField.h"

namespace MpmSim
{

class GravityField : public ForceField
{
public:

	GravityField( const Eigen::Vector3f& gravity );
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	virtual Eigen::Vector3f force( const Eigen::Vector3f& x, float m ) const;
	virtual Eigen::Matrix3f dFdx( const Eigen::Vector3f& x, float m ) const;

private:

	Eigen::Vector3f m_gravity;

};

} //namespace MpmSim

#endif
