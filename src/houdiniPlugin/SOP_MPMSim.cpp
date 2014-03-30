
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <UT/UT_Interrupt.h>
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

#include "MpmSim/Grid.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/ConjugateResiduals.h"



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
    PRM_Name("maxIterations",	"Max Iterations"),
};

static PRM_Default      toleranceDefault(1.e-4);         // Default to 5 divisions
static PRM_Default      iterationsDefault(60);         // Default to 5 divisions


PRM_Template
SOP_MPMSim::myTemplateList[] = {
    PRM_Template(PRM_FLT_J,	1, &names[0], PRMpointOneDefaults),
    PRM_Template(PRM_INT,	1, &names[1], PRMfourDefaults),
    PRM_Template(PRM_INT,	1, &names[2], PRMoneDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[3], &toleranceDefault),
    PRM_Template(PRM_INT,	1, &names[4], &iterationsDefault),
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
	m_snowModel.reset(
		new MpmSim::SnowConstitutiveModel(
			1.4e5f, // young's modulus
			0.2f, // poisson ratio
			10, // hardening
			2.5e-2f, // compressive strength
			7.5e-3f	// tensile strength
		)
	);
}

SOP_MPMSim::~SOP_MPMSim() {}

OP_ERROR
SOP_MPMSim::cookInputGroups(OP_Context &context, int alone)
{
    // The SOP_Node::cookInputPointGroups() provides a good default
    // implementation for just handling a point selection.
    return cookInputPointGroups(context, myGroup, myDetailGroupPair, alone);
}

OP_ERROR
SOP_MPMSim::createParticles(OP_Context &context)
{
	// cook first input with supplied context
	if (lockInput(0, context) >= UT_ERROR_ABORT)
	{
		std::cerr << "couldn't lock input " << std::endl;
		return UT_ERROR_ABORT;
	}
	
	const GU_Detail *initialStateGdp = inputGeo(0);
	const GEO_PrimVDB *vdb;
	std::cerr << "input prims: " <<initialStateGdp->primitives().entries() << std::endl;
	if( initialStateGdp->primitives().entries() == 1 )
	{
        const GEO_Primitive *prim = initialStateGdp->primitives().entry(0);
		if( !prim )
		{
			std::cerr << "prim is null?" << std::endl;
		}
		else
		{
			std::cerr << "yup - got prim 0" << std::endl;
		}
        vdb = dynamic_cast<const GEO_PrimVDB *>(prim);
	}
	
	if( vdb )
	{
		float h = gridSize( context.getTime() ) / 2;
		UT_BoundingBox bbox;
		vdb->getBBox( &bbox );
		
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
					float vdbValue = vdb->getValueF( UT_Vector3( bbox.xmin() + ( i + 0.5 ) * h, bbox.ymin() + ( j + 0.5 ) * h, bbox.zmin() + ( k + 0.5 ) * h ) );
					
					if( vdbValue < 0 )
					{
						Eigen::Vector3f pos( bbox.xmin() + ( i + 0.5 ) * h, bbox.ymin() + ( j + 0.5 ) * h, bbox.zmin() + ( k + 0.5 ) * h );
						x.push_back( pos );
						m.push_back( initialDensity * h * h * h );
					}
				}
			}
		}

		m_particleData.reset( new MpmSim::ParticleData( x, m, 2 * h ) );
		m_snowModel->initParticles( *( m_particleData.get() ) );
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
		std::cerr << "wipe out particles" << std::endl;
		m_particleData.reset(0);
	}
	
	if( !m_particleData.get() )
	{
		std::cerr << "construct particles" << std::endl;
		OP_Context particlesContext = context;
		particlesContext.setTime( startTime );
		if( createParticles( particlesContext ) >= UT_ERROR_ABORT )
		{
			std::cerr << "create particles failed" << std::endl;
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

	if( !m_particleData.get() )
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
			
			std::vector< MpmSimHoudini::VDBCollisionObject > vdbCollisions;
			std::vector< MpmSim::CollisionObject* > collisionObjects;
			
			const GU_Detail *collisionsGdp = inputGeo(1);
			int numPrims = collisionsGdp->primitives().entries();
			for( int i=0; i < numPrims; ++i )
			{
				const GEO_PrimVDB *vdb = dynamic_cast<const GEO_PrimVDB *>( collisionsGdp->primitives().entry(i) );
				if( vdb )
				{
					std::cerr << "add collsion object!" << std::endl;
					vdbCollisions.push_back( MpmSimHoudini::VDBCollisionObject( vdb ) );
				}
			}
			
			std::cerr << vdbCollisions.size() << " vdb collisions!" << std::endl;
			for( size_t i=0; i < vdbCollisions.size(); ++i )
			{
				collisionObjects.push_back( &vdbCollisions[i] );
			}
			
			for( int i=0; i < steps; ++i )
			{
				std::cerr << "update particles " << dt << " " << h << std::endl;

				std::cerr << "construct grid" << std::endl;
				MpmSim::CubicBsplineShapeFunction shapeFunction;
				MpmSim::Grid g(
					*(m_particleData.get()),
					dt,	// time step
					shapeFunction,
					*( m_snowModel.get() )
				);

				if( m_particleData->particleVolumes.empty() )
				{
					std::cerr << "compute particle volumes" << std::endl;
					g.computeDensities( *(m_particleData.get()) );
					m_particleData->particleVolumes.resize( m_particleData->particleX.size() );
					for( size_t i = 0; i < m_particleData->particleDensities.size(); ++i )
					{
						m_particleData->particleVolumes[i] = m_particleData->particleM[i] / m_particleData->particleDensities[i];
					}
				}
				
				// update grid velocities using internal stresses...
				g.updateGridVelocities( *(m_particleData.get()), collisionObjects, MpmSim::ConjugateResiduals( maxIterations(t), tolerance(t) ) );
				
				// transfer the grid velocities back onto the particles:
				g.updateParticleVelocities( *(m_particleData.get()) );
				
				// update particle deformation gradients:
				g.updateDeformationGradients( *(m_particleData.get()) );

				// update positions...
				m_particleData->advance( dt );
				if( boss->opInterrupt() )
				{
					break;
				}
			}
			boss->opEnd();
		}
	}

	std::vector<Eigen::Vector3f>& x = m_particleData->particleX;
	
	std::cerr << "create " << x.size() << " particles" << std::endl;
	
	GU_PrimParticle* particles = GU_PrimParticle::build( gdp, x.size() );
	for( size_t i=0; i < x.size(); ++i )
	{
		gdp->setPos3( i, UT_Vector3(x[i][0],x[i][1],x[i][2]) );
	}
	
	m_prevCookTime = t;
    unlockInputs();
	gdp->notifyCache(GU_CACHE_ALL);
    return error();
}
