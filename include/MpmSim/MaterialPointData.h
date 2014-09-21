#ifndef MPMSIM_MATERIALPOINTDATA_H
#define MPMSIM_MATERIALPOINTDATA_H

#include <vector>
#include <map>
#include <Eigen/Dense>

namespace MpmSim
{

// class containing per particle data of different types
class MaterialPointData
{
public:
	
	MaterialPointData();
	~MaterialPointData();
	
	// create a variable of the specified name with the specified type:
	template<typename T>
	void createVariable( const std::string& name );
	
	// retrieve a variable of the specified name with the specified type:
	template<typename T>
	std::vector<T>& variable( const std::string& name );
	
	template<typename T>
	const std::vector<T>& variable( const std::string& name ) const;
	
	// number of variables:
	size_t numVariables() const;
	
	// returns a list of variable names:
	std::vector<std::string> variableNames() const;
	
private:

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
	
	typedef std::map< std::string, MaterialPointVariableBase* > VariableMap;
	VariableMap m_variables;
};

}

#include "MaterialPointData.inl"

#endif
