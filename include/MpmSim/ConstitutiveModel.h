#ifndef MPMSIM_ConstitutiveMODEL_H
#define MPMSIM_ConstitutiveMODEL_H

#include "Sim.h"

#include <Eigen/Dense>

namespace MpmSim
{

class ConstitutiveModel
{
public:

	virtual ~ConstitutiveModel() {}
	
	// set up particles: must be called before using the following methods
	virtual void setParticles( MaterialPointData& p ) = 0;
	
	// constitutive model specific update step, so you can track things like plastic deformation and
	// history dependant material properties:
	virtual void updateParticleData() = 0;
	
	// energy density for particle p:
	virtual float energyDensity( size_t p ) const = 0;
	
	// matrix of derivatives of energy density with respect to components of the deformation gradient:
	virtual Eigen::Matrix3f dEnergyDensitydF( size_t p ) const = 0;
	
	// computes the change in dEnergyDensitydF when you change F by dFp:
	virtual Eigen::Matrix3f dEdFDifferential( const Eigen::Matrix3f& dFp, size_t p ) const = 0;

};

} // namespace MpmSim

#endif
