#include "tests/TestShapeFunction.h"
#include "tests/TestConjugateResiduals.h"

#include "MpmSim/CubicBsplineShapeFunction.h"

using namespace MpmSim;
using namespace MpmSimTest;

int main(int argc, char** argv)
{
	testShapeFunction( CubicBsplineShapeFunction() );
	testConjugateResiduals();

	return 0;
}
