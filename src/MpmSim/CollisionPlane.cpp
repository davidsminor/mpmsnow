#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>

#include "MpmSim/CollisionPlane.h"

using namespace Eigen;
using namespace MpmSim;

CollisionPlane::CollisionPlane(
	const Eigen::Vector4f& c,
	float coulombFriction,
	bool sticky
) :
	m_coeffs( c ),
	m_coulombFriction( coulombFriction ),
	m_sticky( sticky )
{
	m_v.setZero();
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
	v = m_v;
}

float CollisionPlane::coulombFriction() const
{
	return m_coulombFriction;
}

bool CollisionPlane::sticky() const
{
	return m_sticky;
}

void CollisionPlane::setCoeffs( const Eigen::Vector4f& c )
{
	m_coeffs = c;
}

void CollisionPlane::setV( const Eigen::Vector3f& v )
{
	m_v = v;
}

void CollisionPlane::draw() const
{
	Vector3f basisY = m_coeffs.segment<3>(0);
	Vector3f origin = - basisY * m_coeffs[3] / (basisY.dot(basisY));
	basisY.normalize();

	float mincomponent(1.e10);
	int minidx(0);
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
	
	glEnable( GL_LIGHTING );
	
	GLfloat colorWhite[] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colorWhite);
	
	
	glBegin( GL_QUADS );

	glNormal3f( basisY[0], basisY[1], basisY[2] );

	Vector3f x0 = origin + 5 * (   basisX + basisZ );
	Vector3f x1 = origin + 5 * (   basisX - basisZ );
	Vector3f x2 = origin + 5 * ( - basisX - basisZ );
	Vector3f x3 = origin + 5 * ( - basisX + basisZ );
	
	glVertex3f( x0[0], x0[1], x0[2] );
	glVertex3f( x1[0], x1[1], x1[2] );
	glVertex3f( x2[0], x2[1], x2[2] );
	glVertex3f( x3[0], x3[1], x3[2] );

	glEnd();
	
	
	glDisable( GL_LIGHTING );
	glBegin(GL_LINES);
	glColor3f( 1,1,1 );
	for( int x=0; x <= 10; ++x )
	{
		Vector3f start = origin  + ( x-5.f ) * basisX - 5.f * basisZ;
		Vector3f end = origin  + ( x-5.f ) * basisX + 5.f * basisZ;
		glVertex3f( start[0], start[1], start[2] );
		glVertex3f( end[0], end[1], end[2] );
	}
	for( int x=0; x <= 10; ++x )
	{
		Vector3f start = origin  + ( x-5.f ) * basisZ - 5.f * basisX;
		Vector3f end = origin  + ( x-5.f ) * basisZ + 5.f * basisX;
		glVertex3f( start[0], start[1], start[2] );
		glVertex3f( end[0], end[1], end[2] );
	}
	glEnd();
}
