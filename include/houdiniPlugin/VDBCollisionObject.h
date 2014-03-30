#ifndef MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H
#define MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H

#pragma warning(push, 0) 
#pragma warning( disable : 4503 )
#include <GEO/GEO_PrimVDB.h>
#pragma warning(pop)

#include "MpmSim/CollisionObject.h"

namespace MpmSimHoudini
{

class VDBCollisionObject : public MpmSim::CollisionObject
{

public:
	
	VDBCollisionObject( const GEO_PrimVDB* vdb );

	virtual float phi( const Eigen::Vector3f& x ) const;

	virtual void grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const;
	
	virtual void velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const;

	virtual void draw() const;

private:
	
	const GEO_PrimVDB* m_vdb;

};

} // namespace MpmSim

#endif // MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H
