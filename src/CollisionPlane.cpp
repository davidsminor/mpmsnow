#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>

#include "CollisionPlane.h"

using namespace Eigen;

CollisionPlane::CollisionPlane( const Eigen::Vector4f& c ) : m_coeffs( c )
{
}

float CollisionPlane::phi( const Eigen::Vector3f& x ) const
{
	return m_coeffs[0] * x[0] + m_coeffs[1] * x[1] + m_coeffs[2] * x[2] + m_coeffs[3];
}

void CollisionPlane::grad( const Eigen::Vector3f& x, Eigen::Vector3f& dPhi ) const
{
	dPhi = m_coeffs.segment<3>(0);
}

void CollisionPlane::velocity( const Eigen::Vector3f& x, Eigen::Vector3f& v ) const
{
	v.setZero();
}

void CollisionPlane::draw() const
{
	Vector3f basisY = m_coeffs.segment<3>(0);
	Vector3f origin = - basisY * m_coeffs[3] / (basisY.dot(basisY));
	basisY.normalize();

	float mincomponent(1.e10);
	int minidx;
	for( int i=0; i < 3; ++i )
	{
		if( fabs( basisY[i] ) < mincomponent )
		{
			minidx = i;
			mincomponent = fabs( basisY[i] );
		}
	}
	
	Vector3f basisX = Vector3f::Zero();
	basisX[minidx] = 1.0f;
	Vector3f basisZ = basisX.cross( basisY );
	basisZ.normalize();
	basisX = basisY.cross( basisZ );
	
	glBegin(GL_LINES);
	for( int x=0; x <= 10; ++x )
	{
		Vector3f start = origin  + ( x-5 ) * basisX - 5 * basisZ;
		Vector3f end = origin  + ( x-5 ) * basisX + 5 * basisZ;
		glVertex3f( start[0], start[1], start[2] );
		glVertex3f( end[0], end[1], end[2] );
	}
	for( int x=0; x <= 10; ++x )
	{
		Vector3f start = origin  + ( x-5 ) * basisZ - 5 * basisX;
		Vector3f end = origin  + ( x-5 ) * basisZ + 5 * basisX;
		glVertex3f( start[0], start[1], start[2] );
		glVertex3f( end[0], end[1], end[2] );
	}
	glEnd();
}
