#ifndef SOP_MPMSIM_H
#define SOP_MPMSIM_H

#include <SOP/SOP_Node.h>

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
    void	getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
    fpreal	AMP(fpreal t)		{ return evalFloat("amp", 0, t); }
    fpreal	PHASE(fpreal t)		{ return evalFloat("phase", 0, t); }
    fpreal	PERIOD(fpreal t)	{ return evalFloat("period", 0, t); }

    /// This variable is used together with the call to the "checkInputChanged"
    /// routine to notify the handles (if any) if the input has changed.
    GU_DetailGroupPair	 myDetailGroupPair;

    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup	*myGroup;
};

#endif
