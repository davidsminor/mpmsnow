#ifndef MPMSIM_SHAPEFUNCTION_H
#define MPMSIM_SHAPEFUNCTION_H

#include <vector>

#include <Eigen/Dense>

namespace MpmSim
{

class ShapeFunction
{
public:
	
	class PointToGridIterator
	{
	public:
		
		PointToGridIterator( const ShapeFunction& shapeFunction, float gridH, const Eigen::Vector3f& gridOrigin );
		
		void initialize( const Eigen::Vector3f& p, bool computeDerivatives = false );
		bool next();
		void gridPos( Eigen::Vector3i& pos ) const;
		void dw( Eigen::Vector3f& g ) const;
		float w() const;

	private:

		const ShapeFunction& m_shapeFunction;
		
		float m_gridH;
		Eigen::Vector3f m_gridOrigin;
		int m_diameter;
		
		std::vector<float> m_w[3];
		std::vector<float> m_dw[3];
		Eigen::Vector3i m_pos;
		Eigen::Vector3i m_base;
		bool m_gradients;
		
	};
	
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
