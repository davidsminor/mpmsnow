#ifndef MPMSIM_MATERIALPOINTDATA_H
#define MPMSIM_MATERIALPOINTDATA_H

#include <vector>
#include <map>
#include <Eigen/Dense>

namespace MpmSim
{

class MaterialPointDataBase
{
public:
	virtual size_t dataSize() = 0;
};

template <typename T> 
class MaterialPointData : public MaterialPointDataBase
{
public:
	MaterialPointData() {}
	MaterialPointData( size_t n ) : m_data(n) {}
	MaterialPointData( size_t n, const T& value ) : m_data(n,value) {}
	virtual size_t dataSize() { return m_data.size(); };
	typedef std::vector<T> Data;
	Data m_data;
};

typedef MaterialPointData<float> ScalarData;
typedef MaterialPointData<Eigen::Vector3f> VectorData;
typedef MaterialPointData<Eigen::Matrix3f> MatrixData;

typedef std::map< std::string, MaterialPointDataBase* > MaterialPointDataMap;

}
#endif
