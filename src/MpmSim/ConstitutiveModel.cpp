#include <stdexcept>
#include "MpmSim/ConstitutiveModel.h"

using namespace MpmSim;

const std::vector<float>& ConstitutiveModel::scalarData( const Sim::MaterialPointDataMap& p, const std::string& name )
{
	Sim::MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	const ScalarData* data = dynamic_cast<const ScalarData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find scalar data " + name );
	}
	return data->m_data;
}

const std::vector<Eigen::Vector3f>& ConstitutiveModel::vectorData( const Sim::MaterialPointDataMap& p, const std::string& name )
{
	Sim::MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	const VectorData* data = dynamic_cast<const VectorData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find vector data " + name );
	}
	return data->m_data;
}

const std::vector<Eigen::Matrix3f>& ConstitutiveModel::matrixData( const Sim::MaterialPointDataMap& p, const std::string& name )
{
	Sim::MaterialPointDataMap::const_iterator it = p.find( name );
	if( it == p.end() )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	const MatrixData* data = dynamic_cast<const MatrixData*>( it->second );
	if( !data )
	{
		throw std::runtime_error( "couldn't find matrix data " + name );
	}
	return data->m_data;

}

