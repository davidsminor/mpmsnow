#ifndef MPMSIM_MATERIALPOINTDATA_H
#define MPMSIM_MATERIALPOINTDATA_H

#include <vector>
#include <Eigen/Dense>

namespace MpmSim
{

class MaterialPointDataBase
{
public:
	virtual std::string dataType() = 0;
};

template <typename T> 
class MaterialPointData : public MaterialPointDataBase
{
public:
	virtual std::string dataType() { return "idontknowwhatimdoing"; };
	typedef std::vector<T> Data;
	Data m_data;
};

typedef MaterialPointData<float> ScalarData;
typedef MaterialPointData<Eigen::Vector3f> VectorData;
typedef MaterialPointData<Eigen::Matrix3f> MatrixData;

}
#endif
