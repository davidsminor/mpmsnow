#ifndef MPMSIM_GRID_H
#define MPMSIM_GRID_H

#include "ParticleData.h"
#include "CollisionObject.h"
#include "LinearSolver.h"
#include "ConstituativeModel.h"

#include <Eigen/Dense>

namespace MpmSim
{

class Grid
{

public:

	Grid( const ParticleData& d, float gridH, float timeStep, const ConstituativeModel& model );

	void draw() const;
	void computeDensities( ParticleData& d ) const;

	void updateGridVelocities(
		const ParticleData& d,
		const std::vector<CollisionObject*>& collisionObjects,
		const LinearSolver& implicitSolver );

	void updateDeformationGradients( ParticleData& d );
	void updateParticleVelocities( ParticleData& d );
	
	// testing:
	void testForces( const ParticleData& d );
	void testForceDifferentials( const ParticleData& d );
	void outputDiagnostics( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects ) const;

	void applyImplicitUpdateMatrix( const ParticleData& d, const std::vector<CollisionObject*>& collisionObjects, const Eigen::VectorXf& v, Eigen::VectorXf& result ) const;

private:

	void calculateForces( const ParticleData& d, Eigen::VectorXf& forces ) const;
	void calculateForceDifferentials( const ParticleData& d, const Eigen::VectorXf& dx, Eigen::VectorXf& df ) const;

	// testing
	float calculateEnergy( const ParticleData& d ) const;
	
	// shape functions and derivatives:
	static inline float N( float x );
	static inline float DN( float x );

	static inline void minMax( float x, float& min, float& max );
	inline int fixDim( float& min, float& max ) const;

	inline int coordsToIndex( int x, int y, int z ) const;

	void cellAndWeights( const Eigen::Vector3f& particleX, Eigen::Vector3i& particleCell, float *w[], float** dw = 0 ) const;

private:
	
	float m_gridH;
	float m_timeStep;

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

	const ConstituativeModel& m_constituativeModel;

};

} // namespace MpmSim


#endif // MPMSIM_GRID_H
