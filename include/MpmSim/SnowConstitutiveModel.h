#ifndef MPMSIM_SNOWConstitutiveMODEL_H
#define MPMSIM_SNOWConstitutiveMODEL_H

#include "ConstitutiveModel.h"

namespace MpmSim
{

class SnowConstitutiveModel : public ConstitutiveModel
{
public:
	
	SnowConstitutiveModel(
			float youngsModulus,
			float poissonRatio,
			float hardening,
			float compressiveStrength,
			float tensileStrength );
	
	virtual void setParticles( MaterialPointData& p );
	
	// update deformation at particle p:
	virtual void updateParticleData();

	// energy density for particle p:
	virtual float energyDensity( size_t p ) const;
	
	// matrix of derivatives of energy density with respect to components of the deformation gradient:
	virtual Eigen::Matrix3f dEnergyDensitydF( size_t p ) const;
	
	// computes the change in dEnergyDensitydF when you change F by dFp:
	virtual Eigen::Matrix3f dEdFDifferential( const Eigen::Matrix3f& dFp, size_t p ) const;

private:

	static float matrixDoubleDot( const Eigen::Matrix3f& a, const Eigen::Matrix3f& b );
	static Eigen::Matrix3f computeRdifferential( const Eigen::Matrix3f& dF, const Eigen::Matrix3f& R, const Eigen::Matrix3f& Ginv );
	
	float m_youngsModulus;
	float m_poissonRatio;
	float m_hardening;
	float m_compressiveStrength;
	float m_tensileStrength;

	float m_mu;
	float m_lambda;
	
	MaterialPointData* m_p;

	std::vector<Eigen::Matrix3f>* m_particleF;
	std::vector<Eigen::Matrix3f>* m_particleR;
	std::vector<Eigen::Matrix3f>* m_particleFinvTrans;
	std::vector<Eigen::Matrix3f>* m_particleGinv;
	
	std::vector<float>* m_particleJ;
	std::vector<float>* m_particleMu;
	std::vector<float>* m_particleLambda;
	
};

} //namespace MpmSim

#endif
