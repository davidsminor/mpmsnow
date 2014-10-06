#ifndef MPMSIM_COLLISIONOBJECT_H
#define MPMSIM_COLLISIONOBJECT_H

#include <Eigen/Dense>
#include <vector>

namespace MpmSim
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
	
	// coulomb friction:
	virtual float coulombFriction() const = 0;
	
	// is this sticky?
	virtual bool sticky() const = 0;

	virtual void draw() const = 0;
	
	class CollisionObjectSet
	{
	public:
		~CollisionObjectSet();
		void add( CollisionObject* o );
		const CollisionObject* object( size_t i ) const;
		size_t numObjects() const;
		int collide(
			Eigen::Vector3f& v,
			const Eigen::Vector3f& x,
			const Eigen::Vector3f& frameVelocity,
			bool addCollisionVelocity = false
		) const;
	private:
		std::vector<const CollisionObject*> m_objects;
	};
	
};

} // namespace MpmSim

#endif // MPMSIM_COLLISIONOBJECT_H
