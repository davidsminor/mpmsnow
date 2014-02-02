#ifndef MPMSIM_CUBICBSPLINESHAPEFUNCTION_H
#define MPMSIM_CUBICBSPLINESHAPEFUNCTION_H

#include "MpmSim/ShapeFunction.h"

namespace MpmSim
{

class CubicBsplineShapeFunction : public ShapeFunction
{
public:
	
	// support radius of shape function in cells:
	virtual int supportRadius() const;

	// weight value at normalized coordinate x:
	virtual float w( float x ) const;
	
	// dw/dx:
	virtual float dw( float x ) const;

};

} //namespace MpmSim

#endif

