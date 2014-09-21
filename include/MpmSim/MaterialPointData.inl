#ifndef MPMSIM_MATERIALPOINTDATA_INL
#define MPMSIM_MATERIALPOINTDATA_INL

#include <stdexcept>

namespace MpmSim
{

template<typename T>
void MaterialPointData::createVariable( const std::string& name )
{
	if( m_variables.find( name ) != m_variables.end() )
	{
		throw std::runtime_error( "Material point data already has a variable named '" + name + "'" );
	}
	m_variables[name] = new MaterialPointVariable<T>;
}

template<typename T>
std::vector<T>& MaterialPointData::variable( const std::string& name )
{
	VariableMap::iterator it = m_variables.find( name );
	if( it == m_variables.end() )
	{
		throw std::runtime_error( "Material point data has no variable named '" + name + "'" );
	}
	
	MaterialPointVariable<T>* variable = dynamic_cast< MaterialPointVariable<T>* >( it->second );
	if( !variable )
	{
		throw std::runtime_error( "Material point data has no variable named '" + name + "' of the specified type" );
	}
	return variable->m_data;
}
	
template<typename T>
const std::vector<T>& MaterialPointData::variable( const std::string& name ) const
{
	return const_cast<MaterialPointData*>( this )->variable<T>( name );
}

}
#endif
