#include "MpmSim/MaterialPointData.h"

#include <Eigen/Dense>

using namespace MpmSim;
using namespace Eigen;

MaterialPointData::MaterialPointData()
{
	// create default variables:
	createVariable<Vector3f>("p");
	createVariable<Vector3f>("v");
	createVariable<Matrix3f>("F");
	createVariable<float>("m");
	createVariable<float>("volume");
}

MaterialPointData::~MaterialPointData()
{
	for( VariableMap::iterator it = m_variables.begin(); it != m_variables.end(); ++it )
	{
		delete it->second;
	}
}

size_t MaterialPointData::numVariables() const
{
	return m_variables.size();
}

std::vector<std::string> MaterialPointData::variableNames() const
{
	std::vector<std::string> names;
	for( VariableMap::const_iterator it = m_variables.begin(); it != m_variables.end(); ++it )
	{
		names.push_back( it->first );
	}
	return names;
}
