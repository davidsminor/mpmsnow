#ifndef MPMSNOW_SNOWCONSTITUATIVEMODEL_H
#define MPMSNOW_SNOWCONSTITUATIVEMODEL_H

#include "ConstituativeModel.h"

namespace MpmSnow
{

class SnowConstituativeModel : public ConstituativeModel
{
public:
	
	SnowConstituativeModel(
			float youngsModulus,
			float poissonRatio,
			float hardening,
			float compressiveStrength,
			float tensileStrength );

	virtual void initParticles( ParticleData& p ) const;

	// update deformation at particle p:
	virtual void updateDeformation( ParticleData& d ) const;

	// energy density for particle p:
	virtual float energyDensity( const ParticleData& d, size_t p ) const;
	
	// derivative of energy density with respect to the deformation gradient at particle p:
	virtual void dEnergyDensitydF( Eigen::Matrix3f& result, const ParticleData& d, size_t p ) const;
	
	// energy differentials, using second derivatives of energy function
	virtual void forceDifferentialDensity( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const;

private:

	static float matrixDoubleDot( const Eigen::Matrix3f& a, const Eigen::Matrix3f& b );
	static Eigen::Matrix3f computeRdifferential( const Eigen::Matrix3f& dF, const Eigen::Matrix3f& R, const Eigen::Matrix3f& S );
	
	float m_youngsModulus;
	float m_poissonRatio;
	float m_hardening;
	float m_compressiveStrength;
	float m_tensileStrength;

	float m_mu;
	float m_lambda;

};

} //namespace MpmSnow

#endif
