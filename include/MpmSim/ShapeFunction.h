#ifndef MPMSIM_SHAPEFUNCTION_H
#define MPMSIM_SHAPEFUNCTION_H

#include <vector>

#include <Eigen/Dense>

namespace MpmSim
{

class ShapeFunction
{
public:
	
	virtual ~ShapeFunction();
	
	// support radius of shape function in cells:
	virtual int supportRadius() const = 0;

	// weight value at normalized coordinate x:
	virtual float w( float x ) const = 0;
	
	// dw/dx:
	virtual float dw( float x ) const = 0;

};

} //namespace MpmSim

#endif
