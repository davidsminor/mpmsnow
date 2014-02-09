#include "tests/TestParticleData.h"
#include "tests/TestShapeFunction.h"
#include "tests/TestConjugateResiduals.h"
#include "tests/TestSnowConstitutiveModel.h"

#include "MpmSim/CubicBsplineShapeFunction.h"

using namespace MpmSim;
using namespace MpmSimTest;

int main(int argc, char** argv)
{
	testShapeFunction( CubicBsplineShapeFunction() );
	testConjugateResiduals();
	testSnowConstitutiveModel();
	testParticleData();

	return 0;
}
