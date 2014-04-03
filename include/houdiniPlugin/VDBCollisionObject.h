#ifndef MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H
#define MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H

#include <GEO/GEO_PrimVDB.h>

#include "MpmSim/CollisionObject.h"

namespace MpmSimHoudini
{

class VDBCollisionObject : public MpmSim::CollisionObject
{

public:
	
	VDBCollisionObject( const GEO_PrimVDB* vdb, const GEO_PrimVDB* velocityVdb = 0 );

	virtual float phi( const Eigen::Vector3f& x ) const;

	virtual void grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const;
	
	virtual void velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const;

	virtual void draw() const;

private:
	
	const GEO_PrimVDB* m_vdb;
	const GEO_PrimVDB* m_velocityVdb;

};

} // namespace MpmSim

#endif // MPMSIMHOUDINI_VDBCOLLISIONOBJECT_H
