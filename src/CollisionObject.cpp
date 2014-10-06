#include "MpmSim/CollisionObject.h"

using namespace Eigen;
using namespace MpmSim;

CollisionObject::CollisionObjectSet::~CollisionObjectSet()
{
	for( size_t i=0; i < m_objects.size(); ++i )
	{
		delete m_objects[i];
	}
}

void CollisionObject::CollisionObjectSet::add( CollisionObject* o )
{
	m_objects.push_back( o );
}

int CollisionObject::CollisionObjectSet::collide(
	Vector3f& v,
	const Vector3f& x,
	const Vector3f& frameVelocity,
	bool addCollisionVelocity
) const
{
	int collisionObject(-1);
	bool intersectedAlready(false);
	for( size_t objIdx = 0; objIdx < m_objects.size(); ++objIdx )
	{
		float phi = m_objects[objIdx]->phi( x );
		if( phi <= 0 )
		{
			if( intersectedAlready )
			{
				// colliding with more than one object: set
				// velocity to zero and bail
				collisionObject = -2;
				v.setZero();
				break;
			}
			intersectedAlready = true;

			// intersecting the object
			Vector3f vObj;
			m_objects[objIdx]->velocity( x, vObj );
			
			// express object velocity relative to moving frame:
			vObj -= frameVelocity;

			// subtract off object velocity:
			v -= vObj;
			
			Vector3f n;
			m_objects[objIdx]->grad( x, n );
			n.normalize();
			
			float nDotV = n.dot( v );
			if( nDotV < 0 )
			{
				// trying to move into the object:
				collisionObject = (int)objIdx;

				if( m_objects[objIdx]->sticky() )
				{
					v.setZero();
				}
				else
				{
					
					// velocity perpendicular to the object
					Vector3f vPerp = nDotV * n;
					
					// remaining component is velocity paralell to the object:
					Vector3f vTangent = v - vPerp;
					float vtNorm = vTangent.norm();
					float coulombFriction = m_objects[objIdx]->coulombFriction();
					if( vtNorm >= -nDotV * coulombFriction )
					{
						v = vTangent * ( 1 + coulombFriction * nDotV / vTangent.norm() );
					}
					else
					{
						v.setZero();
					}
				}
			}
			
			if( addCollisionVelocity )
			{
				v += vObj;
			}
		}
	}
	
	return collisionObject;
}


const CollisionObject* CollisionObject::CollisionObjectSet::object( size_t i ) const
{
	return m_objects[i];
}

size_t CollisionObject::CollisionObjectSet::numObjects() const
{
	return m_objects.size();
}
