#include <stdexcept>
#include "MpmSim/ConstitutiveModel.h"

using namespace MpmSim;

std::vector<float>& ConstitutiveModel::scalarData( MaterialPointData& p, const std::string& name )
{
	MaterialPointData::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	ScalarVariable* data = dynamic_cast<ScalarVariable*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	return data->m_data;
}

std::vector<Eigen::Vector3f>& ConstitutiveModel::vectorData( MaterialPointData& p, const std::string& name )
{
	MaterialPointData::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	VectorVariable* data = dynamic_cast<VectorVariable*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	return data->m_data;
}

std::vector<Eigen::Matrix3f>& ConstitutiveModel::matrixData( MaterialPointData& p, const std::string& name )
{
	MaterialPointData::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	MatrixVariable* data = dynamic_cast<MatrixVariable*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	return data->m_data;

}

