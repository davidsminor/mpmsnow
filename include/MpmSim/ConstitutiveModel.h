#ifndef MPMSIM_ConstitutiveMODEL_H
#define MPMSIM_ConstitutiveMODEL_H

#include "ParticleData.h"

#include <Eigen/Dense>

namespace MpmSim
{

class ConstitutiveModel
{
public:
	virtual ~ConstitutiveModel() {}
	
	virtual void initParticles( ParticleData& p ) const = 0;
	
	// update deformation at particle p:
	virtual void updateDeformation( ParticleData& d ) const = 0;

	// energy density for particle p:
	virtual float energyDensity( const ParticleData& d, size_t p ) const = 0;
	
	// matrix of derivatives of energy density with respect to components of the deformation gradient at particle p:
	virtual void dEnergyDensitydF( Eigen::Matrix3f& result, const ParticleData& d, size_t p ) const = 0;
	
	// computes the change in dEnergyDensitydF when you change F by dFp:
	virtual void dEdFDifferential( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const = 0;

};

} // namespace MpmSim

#endif
