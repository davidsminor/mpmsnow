#ifndef MPMSIM_PARTICLEDATA_H
#define MPMSIM_PARTICLEDATA_H

#include <vector>
#include <Eigen/Dense>

namespace MpmSim
{

struct ParticleData
{
	
	ParticleData( const std::vector<Eigen::Vector3f>& x, const std::vector<float>& masses, float gridSize );

	std::vector< Eigen::Vector3f > particleX;
	std::vector< float > particleM;

	std::vector< Eigen::Vector3f > particleV;
	std::vector< Eigen::Matrix3f > particleF;
	std::vector< Eigen::Matrix3f > particleFplastic;
	std::vector< Eigen::Matrix3f > particleR;
	std::vector< Eigen::Matrix3f > particleGinv;
	std::vector< Eigen::Matrix3f > particleFinvTrans;
	
	std::vector< float > particleMu;
	std::vector< float > particleLambda;
	
	std::vector< float > particleJ;

	std::vector< float > particleVolumes;
	std::vector< float > particleDensities;
	
	// indices permuting the particles. All particles in a voxel of size 4 * gridSize
	// will be adjacent in this table
	typedef std::vector<int> IndexList;
	typedef IndexList::iterator IndexIterator;
	IndexList spatialIndex;

	// 8 voxel lists, one for every corner of a cube. When you're splatting particles onto
	// a grid, it's safe to splat voxels in a partition in paralell, as their shape functions
	// will not overlap with other voxels in the partition.
	typedef std::vector< std::pair<IndexIterator, IndexIterator> > PartitionList;
	PartitionList processingPartitions[2][2][2];

	float gridSize;

	void advance( float timeStep );
	void computeProcessingPartitions();
	void voxelSort( std::vector<int>::iterator begin, std::vector<int>::iterator end, int dim = 0 );
};

} // namespace MpmSim

#endif //MPMSIM_PARTICLEDATA_H
