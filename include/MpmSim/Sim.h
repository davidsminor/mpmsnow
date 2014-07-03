#ifndef MPMSIM_SIM_H
#define MPMSIM_SIM_H

#include "MpmSim/MaterialPointData.h"
#include <map>
#include <memory>

namespace MpmSim
{

class Sim
{

public:

	// construct a sim from initial conditions:
	Sim( const std::vector<Eigen::Vector3f>& x, const std::vector<float>& masses, float gridSize );
	
	~Sim();

	size_t numParticleVariables() const;
	
	typedef std::vector<int> IndexList;
	typedef IndexList::iterator IndexIterator;
	typedef std::vector< std::pair<IndexIterator, IndexIterator> > PartitionList;
	
	// partition the sim into contiguous bodies:
	void calculateBodies();
	
	// number of contiguous bodies:
	size_t numBodies() const;
	
	// particle indices in body n:
	const IndexList& body( size_t n ) const;
	
	// checkerboard partitions of body n for paralell processing:
	const PartitionList& bodyPartitionList( size_t n, unsigned i, unsigned j, unsigned k ) const;

	// particle indices for the ballistic particles:
	const IndexList& ballisticParticles() const;
	
	// access point data by name:
	template< class T >
	T* particleVariable( const std::string& name );
	
	template< class T >
	const T* particleVariable( const std::string& name ) const;
	
	// class for doing neighbour queries on particles:
	class NeighbourQuery
	{
		public:
			NeighbourQuery( const std::vector<Eigen::Vector3f>& particleX, float r );
			void neighbours( const Eigen::Vector3f& p, std::vector<int>& neighbourInds ) const;
		private:

			void neighboursInCell( const Eigen::Vector3f& p, const Eigen::Vector3i& cell, std::vector<int>& neighbourInds ) const;

			size_t voxelOffset( const Eigen::Vector3i& cell ) const;
			
			Eigen::Vector3i m_cellMin;
			Eigen::Vector3i m_gridDims;

			const std::vector<Eigen::Vector3f>& m_particleX;
			float m_r;
			IndexList m_spatialSorting;
			std::vector<int> m_voxels;
	};
	

private:
	
	// a list of particles with no neighbours, which are treated in isolation:
	IndexList m_ballisticParticles;
	
	// sort the specified index range into voxels:
	static void voxelSort(
		IndexIterator begin,
		IndexIterator end,
		float voxelSize,
		const std::vector<Eigen::Vector3f>& particleX,
		int dim=0 );

	// the rest of the particles are partitioned into contiguous bodies, whose
	// dynamics are calculated using the material point method:
	struct Body
	{
		// spatially sorted list identifying the particles in this
		// body:
		IndexList particleInds;
		
		// 8 voxel lists, one for every corner of a cube. When you're splatting particles onto
		// a grid, it's safe to splat voxels in a partition in paralell, as their shape functions
		// will not overlap with other voxels in the partition.
		PartitionList processingPartitions[2][2][2];
		void computeProcessingPartitions(
			const std::vector<Eigen::Vector3f>& particleX,
			float voxelSize
		);
	};

	std::vector< Body > m_bodies;
	
	// material point data for all the particles
	typedef std::map< std::string, MaterialPointDataBase* > MaterialPointDataMap;
	MaterialPointDataMap m_pointData;

	// computational grid cell size
	float m_gridSize;

};

} // namespace MpmSim

#include "MpmSim/Sim.inl"

#endif //MPMSIM_PARTICLEDATA_H
