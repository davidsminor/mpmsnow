#ifndef MPMSNOW_GRID_H
#define MPMSNOW_GRID_H

#include "ParticleData.h"
#include "CollisionObject.h"

#include <Eigen/Dense>

#define GRID_H 0.1f
#define TIME_STEP 0.01f
#define INITIALDENSITY 400
#define GRAVITY -9.8f
#define COULOMBFRICTION 0.5f

#define YOUNGSMODULUS 1.4e5f
#define POISSONRATIO 0.2f

#define MU ( YOUNGSMODULUS / ( 2 * ( 1 + POISSONRATIO ) ) )
#define LAMBDA ( YOUNGSMODULUS * POISSONRATIO / ( ( 1 + POISSONRATIO ) * ( 1 - 2 * POISSONRATIO ) ) )
#define BETA 1

class Grid
{

public:

	Grid( const ParticleData& d );

	void draw();
	void computeDensities( ParticleData& d );
	void updateGridVelocities( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects );
	void updateDeformationGradients( ParticleData& d );
	void updateParticleVelocities( ParticleData& d );
	
	// testing:
	void testForces( const ParticleData& d );
	void testForceDifferentials( const ParticleData& d );
	unsigned matrixTexture( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects );

private:

	static float matrixDoubleDot( const Eigen::Matrix3f& a, const Eigen::Matrix3f& b );
	static Eigen::Matrix3f computeRdifferential( const Eigen::Matrix3f& dF, const Eigen::Matrix3f& R, const Eigen::Matrix3f& S );

	void applyImplicitUpdateMatrix( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects, const Eigen::VectorXf& v, Eigen::VectorXf& result );

	// stabilised biconjugate gradient solver copy pasted out of Eigen
	bool bicgstab(
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x,
		int& iters,
		float& tol_error );

	// conjugate gradient method from Eigen:
	void conjugate_gradient(
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const Eigen::VectorXf& rhs,
		Eigen::VectorXf& x,
		int& iters,
		float& tol_error );

	void calculateForces( const ParticleData& d, Eigen::VectorXf& forces );
	void calculateForceDifferentials( const ParticleData& d, const Eigen::VectorXf& dx, Eigen::VectorXf& df );

	// testing
	float calculateEnergy( const ParticleData& d );
	
	// shape functions and derivatives:
	static inline float N( float x );
	static inline float DN( float x );

	inline int coordsToIndex( int x, int y, int z );

	void cellAndWeights( const Eigen::Vector3f& particleX, Eigen::Vector3i& particleCell, float *w[], float** dw = 0 );

	inline void minMax( float x, float& min, float& max );

	inline int fixDim( float& min, float& max );

private:

	float m_xmin;
	float m_ymin;
	float m_zmin;

	float m_xmax;
	float m_ymax;
	float m_zmax;
	
	int m_nx;
	int m_ny;
	int m_nz;

	Eigen::VectorXf m_gridMasses;
	Eigen::VectorXf m_gridVelocities;
	Eigen::VectorXf m_prevGridVelocities;
	std::vector<bool> m_nodeCollided;

};


#endif // MPMSNOW_GRID_H
