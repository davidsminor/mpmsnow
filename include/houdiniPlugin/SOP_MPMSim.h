#ifndef SOP_MPMSIM_H
#define SOP_MPMSIM_H

#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimVDB.h>
#include <GU/GU_Detail.h>

#include "MpmSim/ParticleData.h"
#include "MpmSim/SnowConstitutiveModel.h"

class SOP_MPMSim : public SOP_Node
{
public:
	     SOP_MPMSim(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_MPMSim();

    static PRM_Template		 myTemplateList[];
    static OP_Node		*myConstructor(OP_Network*, const char *,
							    OP_Operator *);

    /// This method is created so that it can be called by handles.  It only
    /// cooks the input group of this SOP.  The geometry in this group is
    /// the only geometry manipulated by this SOP.
    virtual OP_ERROR		 cookInputGroups(OP_Context &context, 
						int alone = 0);

protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);

private:
	fpreal	gridSize(fpreal t)		{ return evalFloat("gridSize", 0, t); }
	int	startFrame(fpreal t)		{ return evalInt("startFrame", 0, t); }
	int	subSteps(fpreal t)		{ return evalInt("subSteps", 0, t); }
	float tolerance(fpreal t)		{ return evalFloat("tolerance", 0, t); }
	int maxIterations(fpreal t)		{ return evalInt("maxIterations", 0, t); }
	
	float youngsModulus(fpreal t)		{ return evalFloat("youngsModulus", 0, t); }
	float poissonRatio(fpreal t)		{ return evalFloat("poissonRatio", 0, t); }
	float hardening(fpreal t)		{ return evalFloat("hardening", 0, t); }
	float compressiveStrength(fpreal t)	{ return evalFloat("compressiveStrength", 0, t); }
	float tensileStrength(fpreal t)		{ return evalFloat("tensileStrength", 0, t); }
	
	void findVDBs( 	const GU_Detail *detail, const GEO_PrimVDB *&pVdb, const GEO_PrimVDB *&vVdb );
	
	fpreal m_prevCookTime;

	OP_ERROR createParticles(OP_Context &context);
	std::auto_ptr< MpmSim::ParticleData > m_particleData;
	std::auto_ptr< MpmSim::SnowConstitutiveModel > m_snowModel;

	/// This variable is used together with the call to the "checkInputChanged"
	/// routine to notify the handles (if any) if the input has changed.
	GU_DetailGroupPair	 myDetailGroupPair;

	/// This is the group of geometry to be manipulated by this SOP and cooked
	/// by the method "cookInputGroups".
	const GA_PointGroup	*myGroup;
};

#endif
