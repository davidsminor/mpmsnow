#ifndef MPMSIM_MATERIALPOINTDATA_H
#define MPMSIM_MATERIALPOINTDATA_H

#include <vector>
#include <map>
#include <Eigen/Dense>

namespace MpmSim
{

class MaterialPointVariableBase
{
public:
	virtual size_t dataSize() = 0;
};

template <typename T> 
class MaterialPointVariable : public MaterialPointVariableBase
{
public:
	MaterialPointVariable() {}
	MaterialPointVariable( size_t n ) : m_data(n) {}
	MaterialPointVariable( size_t n, const T& value ) : m_data(n,value) {}
	virtual size_t dataSize() { return m_data.size(); };
	typedef std::vector<T> Data;
	Data m_data;
};

typedef MaterialPointVariable<float> ScalarVariable;
typedef MaterialPointVariable<Eigen::Vector3f> VectorVariable;
typedef MaterialPointVariable<Eigen::Matrix3f> MatrixVariable;

typedef std::map< std::string, MaterialPointVariableBase* > MaterialPointDataMap;

}
#endif
