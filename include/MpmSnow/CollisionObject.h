#ifndef MPMSNOW_COLLISIONOBJECT_H
#define MPMSNOW_COLLISIONOBJECT_H

#include <Eigen/Dense>

namespace MpmSnow
{

class CollisionObject
{
public:
	
	virtual ~CollisionObject() {}
	
	// return implicit surface value:
	virtual float phi( const Eigen::Vector3f& x ) const = 0;

	// unnormalized gradient for calculating the normal:
	virtual void grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const = 0;
	
	// surface velocity:
	virtual void velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const = 0;
	
	virtual void draw() const = 0;
};

} // namespace MpmSnow

#endif // MPMSNOW_COLLISIONOBJECT_H
