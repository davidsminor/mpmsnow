#include "tests/TestShapeFunction.h"
#include "tests/TestConjugateResiduals.h"
#include "tests/TestSnowConstitutiveModel.h"
#include "tests/TestGrid.h"
#include "tests/TestSimClass.h"

#include "MpmSim/CubicBsplineShapeFunction.h"

using namespace MpmSim;
using namespace MpmSimTest;

int main(int argc, char** argv)
{
	TestShapeFunction::test( CubicBsplineShapeFunction() );
	TestConjugateResiduals::test();
	TestSnowConstitutiveModel::test();
	TestGrid::test();
	TestSimClass::test();
	return 0;
}
