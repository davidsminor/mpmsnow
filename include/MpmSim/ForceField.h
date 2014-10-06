#ifndef MPMSIM_FORCEFIELD_H
#define MPMSIM_FORCEFIELD_H

#include "Eigen/Dense"

#include <vector>

namespace MpmSim
{

class ForceField
{
public:

	virtual ~ForceField() {}

	virtual Eigen::Vector3f force( const Eigen::Vector3f& x, float m ) const = 0;
	virtual Eigen::Matrix3f dFdx( const Eigen::Vector3f& x, float m ) const = 0;
	
	// class that assumes ownership of a list of force fields and accumulates forces:
	class ForceFieldSet
	{
	public:

		~ForceFieldSet();
		
		// add a force field to the list:
		void add( ForceField* f );
		
		// accumulate force on a mass m at point x:
		void force( Eigen::Vector3f& f, const Eigen::Vector3f& x, float m ) const;
		
		// accumulate dFdx on a mass m at point x:
		void dFdx( Eigen::Matrix3f& df, const Eigen::Vector3f& x, float m ) const;
		
	private:
		std::vector<const ForceField*> m_fields;
	};

};

} //namespace MpmSim

#endif
