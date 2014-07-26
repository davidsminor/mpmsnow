#ifndef MPMSIM_COLLISIONPLANE_H
#define MPMSIM_COLLISIONPLANE_H

#include "CollisionObject.h"

namespace MpmSim
{

class CollisionPlane : public CollisionObject
{

public:
	
	// takes a 4 component coefficient vector. Plane equation is:
	// phi = c[0] * x[0] + c[1] * x[1] + c[2] * x[2] + c[3]
	CollisionPlane( const Eigen::Vector4f& c );

	virtual float phi( const Eigen::Vector3f& x ) const;

	virtual void grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const;
	
	virtual void velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const;
	
	virtual float coulombFriction() const;

	virtual bool sticky() const;

	virtual void draw() const;
	
	void setCoeffs( const Eigen::Vector4f& c );
	void setV( const Eigen::Vector3f& v );

private:

	Eigen::Vector4f m_coeffs;
	Eigen::Vector3f m_v;

};

} // namespace MpmSim

#endif // MPMSIM_COLLISIONOBJECT_H
