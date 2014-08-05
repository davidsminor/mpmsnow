
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_VDBUtils.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPart.h>
#include <GEO/GEO_PrimVDB.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_Director.h>
#include <SOP/SOP_Guide.h>

#include "houdiniPlugin/SOP_MPMSim.h"
#include "houdiniPlugin/VDBCollisionObject.h"

#include "MpmSim/SnowConstitutiveModel.h"


void
newSopOperator(OP_OperatorTable *table)
{
     table->addOperator(new OP_Operator("hdk_mpmsnow",
					"MPM SIM",
					 SOP_MPMSim::myConstructor,
					 SOP_MPMSim::myTemplateList,
					 2,
					 2,
					 0));
}

static PRM_Name        names[] = {
    PRM_Name("gridSize",		"Grid Size"),
    PRM_Name("subSteps",		"Sub Steps"),
    PRM_Name("startFrame",		"Start Frame"),
    PRM_Name("tolerance",		"Tolerance"),
    PRM_Name("maxIterations",		"Max Iterations"),
    PRM_Name("youngsModulus",		"Young's Modulus"),
    PRM_Name("poissonRatio",		"Poisson Ratio"),
    PRM_Name("hardening",		"Hardening"),
    PRM_Name("compressiveStrength",	"Compressive Strength"),
    PRM_Name("tensileStrength",		"Tensile Strength"),
};

static PRM_Default      toleranceDefault(1.e-4);         // Default to 5 divisions
static PRM_Default      iterationsDefault(60);         // Default to 5 divisions

static PRM_Default      youngsModulusDefault( 1.4e5f );
static PRM_Default      poissonRatioDefault( 0.2f );
static PRM_Default      hardeningDefault( 10 );
static PRM_Default      compressiveStrengthDefault( 2.5e-2 );
static PRM_Default      tensileStrengthDefault( 7.5e-3 );

PRM_Template
SOP_MPMSim::myTemplateList[] = {
    PRM_Template(PRM_FLT_J,	1, &names[0], PRMpointOneDefaults),
    PRM_Template(PRM_INT,	1, &names[1], PRMfourDefaults),
    PRM_Template(PRM_INT,	1, &names[2], PRMoneDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[3], &toleranceDefault),
    PRM_Template(PRM_INT,	1, &names[4], &iterationsDefault),
    PRM_Template(PRM_FLT_J,	1, &names[5], &youngsModulusDefault),
    PRM_Template(PRM_FLT_J,	1, &names[6], &poissonRatioDefault),
    PRM_Template(PRM_FLT_J,	1, &names[7], &hardeningDefault),
    PRM_Template(PRM_FLT_J,	1, &names[8], &compressiveStrengthDefault),
    PRM_Template(PRM_FLT_J,	1, &names[9], &tensileStrengthDefault),
    PRM_Template(),
};


OP_Node *
SOP_MPMSim::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_MPMSim(net, name, op);
}

SOP_MPMSim::SOP_MPMSim(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op), m_prevCookTime(-1000000), myGroup(0)
{
}

SOP_MPMSim::~SOP_MPMSim() {}

OP_ERROR
SOP_MPMSim::cookInputGroups(OP_Context &context, int alone)
{
    // The SOP_Node::cookInputPointGroups() provides a good default
    // implementation for just handling a point selection.
    return cookInputPointGroups(context, myGroup, myDetailGroupPair, alone);
}


void SOP_MPMSim::findVDBs( const GU_Detail *detail, const GEO_PrimVDB *&pVdb, const GEO_PrimVDB *&vVdb )
{
	int nprims = detail->primitives().entries();
	std::cerr << "input prims: " <<detail->primitives().entries() << std::endl;
	for( int i=0; i < nprims; ++i )
	{
		const GEO_Primitive *prim = detail->primitives().entry(i);
		if( !prim )
		{
			std::cerr << "prim is null?" << std::endl;
			continue;
		}

		const GEO_PrimVDB* vdb = dynamic_cast<const GEO_PrimVDB *>(prim);
		std::cerr << vdb->getGridName() << std::endl;
		switch( vdb->getStorageType() )
		{
			case UT_VDB_INVALID: std::cerr << "UT_VDB_INVALID" << std::endl; break;
			case UT_VDB_FLOAT: std::cerr << "UT_VDB_FLOAT" << std::endl; break;
			case UT_VDB_DOUBLE: std::cerr << "UT_VDB_DOUBLE" << std::endl; break;
			case UT_VDB_INT32: std::cerr << "UT_VDB_INT32" << std::endl; break;
			case UT_VDB_INT64: std::cerr << "UT_VDB_INT64" << std::endl; break;
			case UT_VDB_BOOL: std::cerr << "UT_VDB_BOOL" << std::endl; break;
			case UT_VDB_VEC3F: std::cerr << "UT_VDB_VEC3F" << std::endl; break;
			case UT_VDB_VEC3D: std::cerr << "UT_VDB_VEC3D" << std::endl; break;
			case UT_VDB_VEC3I: std::cerr << "UT_VDB_VEC3I" << std::endl; break;
		}

		if( vdb->getGridName() == std::string("surface") && vdb->getStorageType() == UT_VDB_FLOAT )
		{
			std::cerr << "got surface vdb" << std::endl;
			pVdb = vdb;
		}
		else if( vdb->getGridName() == std::string("v") && vdb->getStorageType() == UT_VDB_VEC3F )
		{
			std::cerr << "got velocity vdb" << std::endl;
			vVdb = vdb;
		}
	}
}

OP_ERROR
SOP_MPMSim::initSim(OP_Context &context)
{
	// cook first input with supplied context
	if (lockInput(0, context) >= UT_ERROR_ABORT)
	{
		std::cerr << "couldn't lock input " << std::endl;
		return UT_ERROR_ABORT;
	}
	
	const GU_Detail *initialStateGdp = inputGeo(0);
	const GEO_PrimVDB *pVdb = 0;
	const GEO_PrimVDB *vVdb = 0;
	
	findVDBs( initialStateGdp, pVdb, vVdb );
	
	if( pVdb )
	{
		fpreal startTime = context.getTime();
		m_snowModel.reset( new MpmSim::SnowConstitutiveModel(
				youngsModulus( startTime ),
				poissonRatio( startTime ),
				hardening( startTime ),
				compressiveStrength( startTime ),
				tensileStrength( startTime )
		) );

		m_forceFields = MpmSim::Sim::ForceFieldSet();
		m_forceFields.fields.push_back( new MpmSim::GravityField( Eigen::Vector3f( 0,-9.8f,0 ) ) );
		
		float h = gridSize( context.getTime() ) / 2;
		UT_BoundingBox bbox;
		pVdb->getBBox( &bbox );
		
		int nx = (int)ceil( bbox.sizeX() / h );
		int ny = (int)ceil( bbox.sizeY() / h );
		int nz = (int)ceil( bbox.sizeZ() / h );
		
		const float initialDensity = 400;
		std::vector<Eigen::Vector3f> x;
		std::vector<float> m;
		for( int i=0; i < nx; ++i )
		{
			for( int j=0; j < ny; ++j )
			{
				for( int k=0; k < nz; ++k )
				{
					float vdbValue = pVdb->getValueF( UT_Vector3( bbox.xmin() + ( i + 0.5 ) * h, bbox.ymin() + ( j + 0.5 ) * h, bbox.zmin() + ( k + 0.5 ) * h ) );
					
					if( vdbValue < 0 )
					{
						Eigen::Vector3f pos( bbox.xmin() + ( i + 0.5 ) * h, bbox.ymin() + ( j + 0.5 ) * h, bbox.zmin() + ( k + 0.5 ) * h );
						x.push_back( pos );
						m.push_back( initialDensity * h * h * h );
					}
				}
			}
		}
		
		m_sim.reset(
			new MpmSim::Sim(
				x, m, gridSize( context.getTime() ), m_shapeFunction, *m_snowModel, m_collisionObjects, m_forceFields
			)
		);
		
		if( vVdb )
		{
			
			const std::vector<Eigen::Vector3f>& x = m_sim->particleVariable<MpmSim::VectorData>( "p" )->m_data;
			std::vector<Eigen::Vector3f>& v = m_sim->particleVariable<MpmSim::VectorData>( "v" )->m_data;
			for( size_t i=0; i < x.size(); ++i )
			{
				UT_Vector3D vdbValue = vVdb->getValueV3( UT_Vector3( x[i][0], x[i][1], x[i][2] ) );
				v[i][0] = vdbValue.x();
				v[i][1] = vdbValue.y();
				v[i][2] = vdbValue.z();
			}
		}
		unlockInput(0);
		return UT_ERROR_NONE;
	}
	else
	{
		std::cerr << "couldn't get vdb" << std::endl;
		opError(OP_BAD_OPINPUT_READ, "Please supply one VDB primitive to input 0");
		unlockInput(0);
		return UT_ERROR_ABORT;
	}
}

OP_ERROR
SOP_MPMSim::cookMySop(OP_Context &context)
{

	std::cerr << "cookMySop" << std::endl;
	OP_Node::flags().timeDep = 1;
	// this is when we evaluate the collision objects:
	fpreal t = context.getTime();
	
	// The context takes seconds, not frame, so we convert.
	fpreal startTime = OPgetDirector()->getChannelManager()->getTime( startFrame(t) );
	
	std::cerr << "time/start: "<< t << "/" << startTime << std::endl;
	
	if( startTime == t )
	{
		std::cerr << "wipe out sim" << std::endl;
		m_sim.reset(0);
	}
	
	if( !m_sim.get() )
	{
		std::cerr << "create sim" << std::endl;
		OP_Context particlesContext = context;
		particlesContext.setTime( startTime );
		if( initSim( particlesContext ) >= UT_ERROR_ABORT )
		{
			std::cerr << "create sim failed" << std::endl;
			unlockInputs();
			return error();
		}
		
		m_prevCookTime = t - 1;
	}
	
	std::cerr << "yeah! lock inputs" << std::endl;
	// Before we do anything, we must lock our inputs.  Before returning,
	//	we have to make sure that the inputs get unlocked.
	if (lockInputs(context) >= UT_ERROR_ABORT)
	{
		std::cerr << "Lock inputs failed" << std::endl;
		return error();
	}

	if( m_prevCookTime >= t )
	{
		std::cerr << "previous cook time >= cook time" << std::endl;
		unlockInputs();
		return error();
	}
	
	std::cerr << "clear gdp" << std::endl;
	gdp->clearAndDestroy();

	if( !m_sim.get() )
	{
		std::cerr << "no particles" << std::endl;
		opError(OP_BAD_OPINPUT_READ, "no particles!");
		unlockInputs();
		return error();
	}
	
	if( t != startTime )
	{
		UT_Interrupt *boss = UTgetInterrupt();
		if (boss->opStart("Advancing Sim..."))
        {
			int steps = subSteps(t);
			float dt = float(t - m_prevCookTime) / steps;
			float h = (float)gridSize(t);
			
			m_collisionObjects = MpmSim::Sim::CollisionObjectSet();
			
			const GU_Detail *collisionsGdp = inputGeo(1);
			
			const GEO_PrimVDB *pVdb = 0;
			const GEO_PrimVDB *vVdb = 0;
			findVDBs( collisionsGdp, pVdb, vVdb );
			
			if( !pVdb )
			{
				opError(OP_BAD_OPINPUT_READ, "Couldn't find VDB for collision SDF!");
				boss->opEnd();
				unlockInputs();
				gdp->notifyCache(GU_CACHE_ALL);
				return error();
			}
			
			if( vVdb )
			{
				std::cerr << "found collision vdb with velocity" << std::endl;
				m_collisionObjects.objects.push_back( new MpmSimHoudini::VDBCollisionObject( pVdb, vVdb ) );
			}
			else
			{
				m_collisionObjects.objects.push_back( new MpmSimHoudini::VDBCollisionObject( pVdb ) );
				std::cerr << "found static collision vdb" << std::endl;
			}
			
			std::cerr << m_collisionObjects.objects.size() << " vdb collisions!" << std::endl;
			
			for( int i=0; i < steps; ++i )
			{
				std::cerr << "update particles " << dt << " " << h << " (" << i+1 << " of " << steps << ")" << std::endl;
				try
				{
					MpmSim::ConjugateResiduals solver( maxIterations(t), tolerance(t) );
					m_sim->advance( dt, solver );
					if( boss->opInterrupt() )
					{
						break;
					}
				}
				catch( const std::exception& e )
				{
					std::cerr << e.what() << std::endl;
				}
			}
			boss->opEnd();
		}
	}

	const std::vector<Eigen::Vector3f>& x = m_sim->particleVariable<MpmSim::VectorData>( "p" )->m_data;
	
	std::cerr << "create " << x.size() << " particles" << std::endl;
	
	GU_PrimParticle::build( gdp, x.size() );
	for( size_t i=0; i < x.size(); ++i )
	{
		gdp->setPos3( i, UT_Vector3(x[i][0],x[i][1],x[i][2]) );
	}
	
	m_prevCookTime = t;
	unlockInputs();
	gdp->notifyCache(GU_CACHE_ALL);
	return error();
}
