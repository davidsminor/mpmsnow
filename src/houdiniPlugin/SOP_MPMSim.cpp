
#define M_PI 3.14159265359

#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include "houdiniPlugin/SOP_MPMSim.h"


void
newSopOperator(OP_OperatorTable *table)
{
     table->addOperator(new OP_Operator("hdk_mpmsnow",
					"MPM SIM",
					 SOP_MPMSim::myConstructor,
					 SOP_MPMSim::myTemplateList,
					 1,
					 1,
					 0));
}

static PRM_Name        names[] = {
    PRM_Name("amp",	"Amplitude"),
    PRM_Name("phase",	"Phase"),
    PRM_Name("period",	"Period"),
};

PRM_Template
SOP_MPMSim::myTemplateList[] = {
    PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu),
    PRM_Template(PRM_FLT_J,	1, &names[0], PRMoneDefaults, 0,
				   &PRMscaleRange),
    PRM_Template(PRM_FLT_J,	1, &names[1], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[2], PRMoneDefaults),
    PRM_Template(),
};


OP_Node *
SOP_MPMSim::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_MPMSim(net, name, op);
}

SOP_MPMSim::SOP_MPMSim(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op), myGroup(0)
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

OP_ERROR
SOP_MPMSim::cookMySop(OP_Context &context)
{
    // Before we do anything, we must lock our inputs.  Before returning,
    //	we have to make sure that the inputs get unlocked.
    if (lockInputs(context) >= UT_ERROR_ABORT)
	return error();

    // Duplicate our incoming geometry with the hint that we only
    // altered points.  Thus if we our input was unchanged we can
    // easily roll back our changes by copying point values.
    duplicatePointSource(0, context);

    fpreal t = context.getTime();

    // We evaluate our parameters outside the loop for speed.  If we
    // wanted local variable support, we'd have to do more setup
    // (see SOP_Flatten) and also move these inside the loop.
    float phase = PHASE(t);
    float amp = AMP(t);
    float period = PERIOD(t);

    // Here we determine which groups we have to work on.  We only
    // handle point groups.
    if (error() < UT_ERROR_ABORT && cookInputGroups(context) < UT_ERROR_ABORT)
    {
        GA_Offset ptoff;
	GA_FOR_ALL_GROUP_PTOFF(gdp, myGroup, ptoff)
	{
	    UT_Vector3 p = gdp->getPos3(ptoff);

	    p.y() += SYSsin( (p.x() / period + phase) * M_PI * 2 ) * amp;

	    gdp->setPos3(ptoff, p);
	}
    }

    unlockInputs();
    
    return error();
}
