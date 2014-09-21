#include <stdexcept>
#include "MpmSim/ConstitutiveModel.h"

using namespace MpmSim;

std::vector<float>& ConstitutiveModel::scalarData( MaterialPointDataMap& p, const std::string& name )
{
	MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	ScalarData* data = dynamic_cast<ScalarData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	return data->m_data;
}

std::vector<Eigen::Vector3f>& ConstitutiveModel::vectorData( MaterialPointDataMap& p, const std::string& name )
{
	MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	VectorData* data = dynamic_cast<VectorData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	return data->m_data;
}

std::vector<Eigen::Matrix3f>& ConstitutiveModel::matrixData( MaterialPointDataMap& p, const std::string& name )
{
	MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	MatrixData* data = dynamic_cast<MatrixData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	return data->m_data;

}

