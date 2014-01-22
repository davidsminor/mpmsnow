#ifndef MPMSNOW_PARTICLEDATA_H
#define MPMSNOW_PARTICLEDATA_H

#include <vector>
#include <Eigen/Dense>

struct ParticleData
{
	std::vector< Eigen::Vector3f > particleX;
	std::vector< Eigen::Vector3f > particleV;
	std::vector< float > particleM;

	std::vector< Eigen::Matrix3f > particleF;
	std::vector< Eigen::Matrix3f > particleFplastic;
	std::vector< Eigen::Matrix3f > particleR;
	std::vector< Eigen::Matrix3f > particleS;
	std::vector< Eigen::Matrix3f > particleFinvTrans;
	
	std::vector< float > particleMu;
	std::vector< float > particleLambda;
	
	std::vector< float > particleJ;

	std::vector< float > particleVolumes;
	std::vector< float > particleDensities;
};

#endif //MPMSNOW_PARTICLEDATA_H
