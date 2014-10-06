#include "MpmSim/SnowConstitutiveModel.h"

#include <stdexcept>
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

void SnowConstitutiveModel::updateParticleData()
{
	std::vector<Eigen::Matrix3f>& particleF = m_p->variable<Matrix3f>( "F" );
	std::vector<Eigen::Matrix3f>& particleFplastic = m_p->variable<Matrix3f>( "Fp" );
	std::vector<Eigen::Matrix3f>& particleFinvTrans = m_p->variable<Matrix3f>( "FinvTrans" );
	std::vector<Eigen::Matrix3f>& particleR = m_p->variable<Matrix3f>( "R" );
	std::vector<Eigen::Matrix3f>& particleGinv = m_p->variable<Matrix3f>( "Ginv" );
	
	std::vector<float>& particleJ = m_p->variable<float>( "J" );
	std::vector<float>& particleMu = m_p->variable<float>( "mu" );
	std::vector<float>& particleLambda = m_p->variable<float>( "lambda" );
	
	for( size_t p=0; p < particleF.size(); ++p )
	{
		JacobiSVD<Matrix3f> svd(particleF[p], ComputeFullU | ComputeFullV );

		Vector3f singularValues = svd.singularValues();

		// apply plastic yeild:
		Matrix3f diagonalMat = Matrix3f::Zero();
		Matrix3f diagonalMatInv = Matrix3f::Zero();
		bool modifiedSVD = false;
		for( int i=0; i < 3; ++i )
		{
			// stretching:
			if( singularValues[i] > 1 + m_tensileStrength )
			{
				modifiedSVD = true;
				singularValues[i] = 1 + m_tensileStrength;
			}

			// compression:
			if( singularValues[i] < 1 - m_compressiveStrength )
			{
				modifiedSVD = true;
				singularValues[i] = 1 - m_compressiveStrength;
			}
			diagonalMat(i,i) = singularValues[i];
			diagonalMatInv(i,i) = 1.0f / singularValues[i];
		}
		
		if( modifiedSVD )
		{
			Matrix3f FNplusOne = particleF[p] * particleFplastic[p];
			particleFplastic[p] = svd.matrixV() * diagonalMatInv * svd.matrixU().transpose() * FNplusOne;
			particleF[p] = svd.matrixU() * diagonalMat * svd.matrixV().transpose();
		}

		particleFinvTrans[p] = svd.matrixU() * diagonalMatInv * svd.matrixV().transpose();
		particleR[p] = svd.matrixU() * svd.matrixV().transpose();
		
		Matrix3f S = svd.matrixV() * diagonalMat * svd.matrixV().transpose();
		Matrix3f G;
		G(0,0) = S(0,0) + S(1,1);
		G(1,1) = S(0,0) + S(2,2);
		G(2,2) = S(1,1) + S(2,2);
		
		G(0,1) = G(1,0) = S(1,2);
		G(0,2) = G(2,0) = -S(0,2);
		G(1,2) = G(2,1) = S(0,1);
		particleGinv[p] = G.inverse();

		particleJ[p] = diagonalMat(0,0) * diagonalMat(1,1) * diagonalMat(2,2);
		
		
		// apply hardening:
		float hardeningFactor = m_hardening * ( 1 - particleFplastic[p].determinant() );
		if( hardeningFactor > 2 )
		{
			// don't let it harden by more than a factor of about 7.4
			hardeningFactor = 2;
		}
		hardeningFactor = exp( hardeningFactor );
		#ifdef WIN32
		if( !_finite(hardeningFactor) )
		#else
		if( isinff(hardeningFactor) || isnanf(hardeningFactor) )
		#endif
		{
			std::cerr << "plastic deformation: " << std::endl << particleFplastic[p] << std::endl;
			std::cerr << "det: " << std::endl << particleFplastic[p].determinant() << std::endl;
			std::cerr << "my log is hardening: " << m_hardening * ( 1 - particleFplastic[p].determinant() ) << std::endl;
			throw std::runtime_error( "infinite hardness!!!" );
		}

		
		particleMu[p] = m_mu * hardeningFactor;
		particleLambda[p] = m_lambda * hardeningFactor;
		
		if( particleJ[p] <= 0 )
		{
			std::cerr << "warning: inverted deformation gradient!" << std::endl;
		}
	}
}

void SnowConstitutiveModel::setParticles( MaterialPointData& p )
{
	m_p = &p;
	
	
	// create variables:
	// plastic deformation:
	m_p->createVariable<Matrix3f>("Fp");

	// inverse transpose of the deformation gradient:
	m_p->createVariable<Matrix3f>("FinvTrans");
	
	// jacobian of the deformation gradient:
	m_p->createVariable<float>("J");

	// rotational component of the deformation gradient
	m_p->createVariable<Matrix3f>("R");

	// matrix used to calculate the differential of R when you change F by a small amount
	m_p->createVariable<Matrix3f>("Ginv");
	
	// lame parameters
	m_p->createVariable<float>("mu");
	m_p->createVariable<float>("lambda");
	
	// initialize variables:
	size_t nParticles = m_p->variable<float>("m").size();
	m_p->variable<Matrix3f>("Fp").resize( nParticles, Eigen::Matrix3f::Identity() );
	m_p->variable<Matrix3f>("FinvTrans").resize( nParticles, Eigen::Matrix3f::Identity() );
	m_p->variable<Matrix3f>("R").resize( nParticles, Eigen::Matrix3f::Identity() );
	m_p->variable<Matrix3f>("Ginv").resize( nParticles, Eigen::Matrix3f::Identity() );
	
	m_p->variable<float>("J").resize( nParticles, 1.0f );
	m_p->variable<float>("mu").resize( nParticles, m_mu );
	m_p->variable<float>("lambda").resize( nParticles, m_lambda );
	
	
	// keep some pointers around for convenience/performance reasons:
	m_particleF = &p.variable<Matrix3f>( "F" );
	m_particleFinvTrans = &p.variable<Matrix3f>( "FinvTrans" );
	m_particleR = &p.variable<Matrix3f>( "R" );
	m_particleGinv = &p.variable<Matrix3f>( "Ginv" );
	
	m_particleJ = &p.variable<float>( "J" );
	m_particleMu = &p.variable<float>( "mu" );
	m_particleLambda = &p.variable<float>( "lambda" );
}

float SnowConstitutiveModel::energyDensity( size_t p ) const
{
	const std::vector<Eigen::Matrix3f>& particleF = *m_particleF;
	const std::vector<Eigen::Matrix3f>& particleR = *m_particleR;
	
	const std::vector<float>& particleJ = *m_particleJ;
	const std::vector<float>& particleMu = *m_particleMu;
	const std::vector<float>& particleLambda = *m_particleLambda;
	
	Matrix3f rigidDeviation = particleF[p] - particleR[p];
	float JminusOne = particleJ[p] - 1;
	float doubledot = matrixDoubleDot( rigidDeviation, rigidDeviation );
	float ret = ( particleMu[p] * doubledot + 0.5f * particleLambda[p] * JminusOne * JminusOne );
	return ret;
}

Eigen::Matrix3f SnowConstitutiveModel::dEnergyDensitydF( size_t p ) const
{
	const std::vector<Eigen::Matrix3f>& particleF = *m_particleF;
	const std::vector<Eigen::Matrix3f>& particleFinvTrans = *m_particleFinvTrans;
	const std::vector<Eigen::Matrix3f>& particleR = *m_particleR;
	
	const std::vector<float>& particleJ = *m_particleJ;
	const std::vector<float>& particleMu = *m_particleMu;
	const std::vector<float>& particleLambda = *m_particleLambda;

	Matrix3f rigidDeviation = particleF[p] - particleR[p];
	Matrix3f ret = 2 * particleMu[p] * (rigidDeviation ) + particleLambda[p] * ( particleJ[p] - 1 ) * particleJ[p] * particleFinvTrans[p];
	
	float n = ret.norm();
	#ifdef WIN32
	if( !_finite(n) )
	#else
	if( isinff(n) || isnanf(n) )
	#endif
	{
		std::cerr << "particle def grad: " << std::endl << particleF[p] << std::endl;
		std::cerr << "inverse transpose: " << std::endl << particleFinvTrans[p] << std::endl;
		std::cerr << "rigidDeviation: " << std::endl << rigidDeviation << std::endl;
		std::cerr << "determinant: " << particleJ[p] << std::endl;
		std::cerr << "lame parameters: " << particleMu[p] << " " << particleLambda[p] << std::endl;
		throw std::runtime_error( "nans in dEnergyDensitydF matrix!" );
	}
	
	return ret;
}

Eigen::Matrix3f SnowConstitutiveModel::dEdFDifferential( const Eigen::Matrix3f& dFp, size_t p ) const
{
	const std::vector<Eigen::Matrix3f>& particleFinvTrans = *m_particleFinvTrans;
	const std::vector<Eigen::Matrix3f>& particleR = *m_particleR;
	const std::vector<Eigen::Matrix3f>& particleGinv = *m_particleGinv;
	
	const std::vector<float>& particleJ = *m_particleJ;
	const std::vector<float>& particleMu = *m_particleMu;
	const std::vector<float>& particleLambda = *m_particleLambda;
	

	// work out energy derivatives with respect to the deformation gradient at this particle:
	// Ap = d2Psi / dF dF : dF (see the tech report). We've got dF, so just plug that into the
	// formulae...
	
	// if you look in dEnergyDensitydF, you'll find it's computing this:
	
	// 2 * MU * ( particleF[p] - particleR[p] ) + LAMBDA * ( particleJ[p] - 1 ) * particleJ[p] * particleFinvTrans[p];
	
	// what we're doing here is just assuming dFp is small and working out the corresponding variation in
	// that expression...
	
	float J = particleJ[p];

	// work out a couple of basic differentials:
	float dJ = J * matrixDoubleDot( particleFinvTrans[p], dFp );
	Matrix3f dFInvTrans = - particleFinvTrans[p] * dFp.transpose() * particleFinvTrans[p];
	
	Matrix3f dR = computeRdifferential( dFp, particleR[p], particleGinv[p] );
	
	return
		// start with differential of 2 * MU * ( F - R )...
			2 * particleMu[p] * ( dFp - dR )
		// add on differential of LAMBDA * ( J - 1 ) * J * F^-t
		// = LAMBDA * ( d( J - 1 ) * J F^-T + ( J - 1 ) * d( J F^-t ) )
		// = LAMBDA * ( dJ * J F^-T + ( J - 1 ) * ( dJ F^-t + J * d( F^-t ) )
			+ particleLambda[p] * ( dJ * J * particleFinvTrans[p] + ( J - 1 ) * ( dJ * particleFinvTrans[p] + J * dFInvTrans ) );

}

float SnowConstitutiveModel::matrixDoubleDot( const Matrix3f& a, const Matrix3f& b )
{
	return
		a(0,0) * b(0,0) + a(0,1) * b(0,1) + a(0,2) * b(0,2) + 
		a(1,0) * b(1,0) + a(1,1) * b(1,1) + a(1,2) * b(1,2) + 
		a(2,0) * b(2,0) + a(2,1) * b(2,1) + a(2,2) * b(2,2);
}

Matrix3f SnowConstitutiveModel::computeRdifferential( const Matrix3f& dF, const Matrix3f& R, const Matrix3f& Ginv )
{
	// \todo: compute w without computing the whole of M
	Matrix3f M = R.transpose() * dF;
	Vector3f w( M(0,1) - M(1,0), M(0,2) - M(2,0), M(1,2) - M(2,1) );
	
	w = Ginv * w;
	
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
