#include "MpmSim/CubicBsplineShapeFunction.h"

using namespace MpmSim;
using namespace Eigen;

int CubicBsplineShapeFunction::supportRadius() const
{
	return 2;
}

float CubicBsplineShapeFunction::w( float x ) const
{
	float ax = fabs(x);
	if( ax < 1 )
	{
		return ( 0.5f * ax - 1.0f ) * ax * ax + 2.0f/3;
	}
	else if( ax < 2 )
	{
		return ( ( -ax / 6 + 1.0f ) *ax - 2.0f ) * ax + 4.0f/3;
	}
	else
	{
		return 0;
	}
}
	
float CubicBsplineShapeFunction::dw( float x ) const
{
	if( x < 0 )
	{
		return -dw( -x );
	}
	
	if( x < 1 )
	{
		return x * ( 1.5f * x - 2 );
	}
	else if( x < 2 )
	{
		x -= 2;
		return -0.5f * x * x;
	}
	else
	{
		return 0;
	}
}
