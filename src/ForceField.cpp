#include "MpmSim/ForceField.h"

using namespace Eigen;
using namespace MpmSim;

ForceField::ForceFieldSet::~ForceFieldSet()
{
	for( size_t i=0; i < m_fields.size(); ++i )
	{
		delete m_fields[i];
	}
}

void ForceField::ForceFieldSet::add( ForceField* f )
{
	m_fields.push_back( f );
}

void ForceField::ForceFieldSet::force( Eigen::Vector3f& f, const Eigen::Vector3f& x, float m ) const
{
	for( size_t i=0; i < m_fields.size(); ++i )
	{
		f += m_fields[i]->force( x, m );
	}
}

void ForceField::ForceFieldSet::dFdx( Eigen::Matrix3f& df, const Eigen::Vector3f& x, float m ) const
{
	for( size_t i=0; i < m_fields.size(); ++i )
	{
		df += m_fields[i]->dFdx( x, m );
	}
}