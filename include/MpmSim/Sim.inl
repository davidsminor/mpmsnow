#ifndef MPMSIM_SIM_INL
#define MPMSIM_SIM_INL

namespace MpmSim
{

template< class T >
T* Sim::particleVariable( const std::string& name )
{
	MaterialPointDataMap::iterator it = particleData.find( name );
	if( it == particleData.end() )
	{
		return 0;
	}
	return dynamic_cast<T*>( it->second );
}

template< class T >
const T* Sim::particleVariable( const std::string& name ) const
{
	return const_cast<Sim*>(this)->particleVariable<T>( name );
}

} // namespace MpmSim

#endif //MPMSIM_SIM_INL
