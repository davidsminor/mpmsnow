#include "MpmSim/GravityField.h"

using namespace MpmSim;

GravityField::GravityField( const Eigen::Vector3f& gravity ) : m_gravity( gravity )
{
}

Eigen::Vector3f GravityField::force( const Eigen::Vector3f& x, float m ) const
{
	return m_gravity * m;
}

Eigen::Matrix3f GravityField::dFdx( const Eigen::Vector3f& x, float m ) const
{
	return Eigen::Matrix3f::Zero();
}
