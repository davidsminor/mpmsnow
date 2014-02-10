#include "MpmSim/SnowConstitutiveModel.h"

#include <iostream>

using namespace Eigen;
using namespace MpmSim;

SnowConstitutiveModel::SnowConstitutiveModel(
	float youngsModulus,
	float poissonRatio,
	float hardening,
	float compressiveStrength,
	float tensileStrength
) :
	m_youngsModulus( youngsModulus ),
	m_poissonRatio( poissonRatio ),
	m_hardening( hardening ),
	m_compressiveStrength( compressiveStrength ),
	m_tensileStrength( tensileStrength )
{
	// calculate Lame parameters:
	m_mu = ( youngsModulus / ( 2 * ( 1 + poissonRatio ) ) );
	m_lambda = ( youngsModulus * poissonRatio / ( ( 1 + poissonRatio ) * ( 1 - 2 * poissonRatio ) ) );
}

void SnowConstitutiveModel::updateDeformation( ParticleData& d ) const
{
	for( size_t p=0; p < d.particleX.size(); ++p )
	{

		JacobiSVD<Matrix3f> svd(d.particleF[p], ComputeFullU | ComputeFullV );

		Vector3f singularValues = svd.singularValues();

		// apply plastic yeild:
		Matrix3f diagonalMat = Matrix3f::Zero();
		Matrix3f diagonalMatInv = Matrix3f::Zero();
		for( int i=0; i < 3; ++i )
		{
			// stretching:
			if( singularValues[i] > 1 + m_tensileStrength )
			{
				singularValues[i] = 1 + m_tensileStrength;
			}

			// compression:
			if( singularValues[i] < 1 - m_compressiveStrength )
			{
				singularValues[i] = 1 - m_compressiveStrength;
			}
			diagonalMat(i,i) = singularValues[i];
			diagonalMatInv(i,i) = 1.0f / singularValues[i];
		}
		
		Matrix3f FNplusOne = d.particleF[p] * d.particleFplastic[p];
		
		d.particleFplastic[p] = svd.matrixV() * diagonalMatInv * svd.matrixU().transpose() * FNplusOne;
		d.particleF[p] = svd.matrixU() * diagonalMat * svd.matrixV().transpose();
		d.particleFinvTrans[p] = svd.matrixU() * diagonalMatInv * svd.matrixV().transpose();
		d.particleR[p] = svd.matrixU() * svd.matrixV().transpose();
		d.particleS[p] = svd.matrixV() * diagonalMat * svd.matrixV().transpose();
		d.particleJ[p] = diagonalMat(0,0) * diagonalMat(1,1) * diagonalMat(2,2);
		
		
		// apply hardening:
		float hardeningFactor = exp( m_hardening * ( 1 - d.particleFplastic[p].determinant() ) );
		if( hardeningFactor > 1 )
		{
			hardeningFactor = 1;
		}

		d.particleMu[p] = m_mu * hardeningFactor;
		d.particleLambda[p] = m_lambda * hardeningFactor;
		
		if( d.particleJ[p] <= 0 )
		{
			std::cerr << "warning: inverted deformation gradient!" << std::endl;
		}
	}
}

void SnowConstitutiveModel::initParticles( ParticleData& d ) const
{
	size_t nParticles = d.particleX.size();
	d.particleFinvTrans.resize( nParticles, Eigen::Matrix3f::Identity() );
	d.particleJ.resize( nParticles, 1.0f );
	d.particleFplastic.resize( nParticles, Eigen::Matrix3f::Identity() );
	d.particleR.resize( nParticles, Eigen::Matrix3f::Identity() );
	d.particleS.resize( nParticles, Eigen::Matrix3f::Identity() );
	d.particleMu.resize( nParticles, m_mu );
	d.particleLambda.resize( nParticles, m_lambda );
}

float SnowConstitutiveModel::energyDensity( const ParticleData& d, size_t p ) const
{
	Matrix3f rigidDeviation = d.particleF[p] - d.particleR[p];
	float JminusOne = d.particleJ[p] - 1;
	return ( d.particleMu[p] * matrixDoubleDot( rigidDeviation, rigidDeviation ) + 0.5f * d.particleLambda[p] * JminusOne * JminusOne );
}
	
void SnowConstitutiveModel::dEnergyDensitydF( Eigen::Matrix3f& result, const ParticleData& d, size_t p ) const
{
	result = 2 * d.particleMu[p] * ( d.particleF[p] - d.particleR[p] ) + d.particleLambda[p] * ( d.particleJ[p] - 1 ) * d.particleJ[p] * d.particleFinvTrans[p];
}

void SnowConstitutiveModel::forceDifferentialDensity( Eigen::Matrix3f& result, const Eigen::Matrix3f& dFp, const ParticleData& d, size_t p ) const
{
	
		// work out energy derivatives with respect to the deformation gradient at this particle:
		// Ap = d2Psi / dF dF : dF (see the tech report). We've got dF, so just plug that into the
		// formulae...
		
		// if you look in dEnergyDensitydF, you'll find it's computing this:
		
		// 2 * MU * ( d.particleF[p] - d.particleR[p] ) + LAMBDA * ( d.particleJ[p] - 1 ) * d.particleJ[p] * d.particleFinvTrans[p];
		
		// what we're doing here is just assuming dFp is small and working out the corresponding variation in
		// that expression...
		
		float J = d.particleJ[p];

		// work out a couple of basic differentials:
		float dJ = J * matrixDoubleDot( d.particleFinvTrans[p], dFp );
		Matrix3f dFInvTrans = - d.particleFinvTrans[p] * dFp.transpose() * d.particleFinvTrans[p];
		
		Matrix3f dR = computeRdifferential( dFp, d.particleR[p], d.particleS[p] );
		
		// start with differential of 2 * MU * ( F - R )...
		result = 2 * d.particleMu[p] * ( dFp - dR );
		
		// add on differential of LAMBDA * ( J - 1 ) * J * F^-t
		// = LAMBDA * ( d( J - 1 ) * J F^-T + ( J - 1 ) * d( J F^-t ) )
		// = LAMBDA * ( dJ * J F^-T + ( J - 1 ) * ( dJ F^-t + J * d( F^-t ) )
		result += d.particleLambda[p] * ( dJ * J * d.particleFinvTrans[p] + ( J - 1 ) * ( dJ * d.particleFinvTrans[p] + J * dFInvTrans ) );
		
}


float SnowConstitutiveModel::matrixDoubleDot( const Matrix3f& a, const Matrix3f& b )
{
	return
		a(0,0) * b(0,0) + a(0,1) * b(0,1) + a(0,2) * b(0,2) + 
		a(1,0) * b(1,0) + a(1,1) * b(1,1) + a(1,2) * b(1,2) + 
		a(2,0) * b(2,0) + a(2,1) * b(2,1) + a(2,2) * b(2,2);
}

Matrix3f SnowConstitutiveModel::computeRdifferential( const Matrix3f& dF, const Matrix3f& R, const Matrix3f& S )
{
	Matrix3f M = R.transpose() * dF - dF.transpose() * R;
	Vector3f w( M(0,1), M(0,2), M(1,2) );
	
	Matrix3f G;
	G(0,0) = S(0,0) + S(1,1);
	G(1,1) = S(0,0) + S(2,2);
	G(2,2) = S(1,1) + S(2,2);
	
	G(0,1) = G(1,0) = S(1,2);
	G(0,2) = G(2,0) = -S(0,2);
	G(1,2) = G(2,1) = S(0,1);
	
	w = G.inverse() * w;
	
	Matrix3f RtdR;
	RtdR(0,0) = RtdR(1,1) = RtdR(2,2) = 0;
	
	RtdR(0,1) = w[0];
	RtdR(1,0) = -w[0];
	
	RtdR(0,2) = w[1];
	RtdR(2,0) = -w[1];
	
	RtdR(1,2) = w[2];
	RtdR(2,1) = -w[2];
	
	return R * RtdR;
}
