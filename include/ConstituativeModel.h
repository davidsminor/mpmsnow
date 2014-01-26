#ifndef MPMSNOW_CONSTITUATIVEMODEL_H
#define MPMSNOW_CONSTITUATIVEMODEL_H

#include "ParticleData.h"

#include <Eigen/Dense>

class ConstituativeModel
{
public:
	
	virtual void initParticles( ParticleData& p ) const = 0;

	// energy density for particle p:
	virtual float energyDensity( const ParticleData& d, size_t p ) const = 0;
	
	// derivative of energy density with respect to the deformation gradient at particle p:
	virtual void dEnergyDensitydF( Eigen::Matrix3f& result, const ParticleData& d, size_t p ) const = 0;
	
	// energy differentials , using second derivatives of energy function
	virtual void forceDifferentialDensity( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const = 0;

};

#endif
