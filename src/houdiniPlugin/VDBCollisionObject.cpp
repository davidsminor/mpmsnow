#include "houdiniPlugin/VDBCollisionObject.h"

using namespace Eigen;
using namespace MpmSim;
using namespace MpmSimHoudini;

VDBCollisionObject::VDBCollisionObject( const GEO_PrimVDB* vdb ) : m_vdb( vdb )
{
}

float VDBCollisionObject::phi( const Eigen::Vector3f& x ) const
{
	return m_vdb->getValueF( UT_Vector3( x[0], x[1], x[2] ) );
}

void VDBCollisionObject::grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const
{
	UT_Vector3 grad = m_vdb->getGradient( UT_Vector3( x[0], x[1], x[2] ) );
	dPhi[0] = grad.x();
	dPhi[1] = grad.y();
	dPhi[2] = grad.z();
}

void VDBCollisionObject::velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const
{
	v.setZero();
}

void VDBCollisionObject::draw() const
{
}
