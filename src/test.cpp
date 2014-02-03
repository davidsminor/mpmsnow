#include "tests/TestShapeFunction.h"
#include "MpmSim/CubicBsplineShapeFunction.h"

using namespace MpmSim;
using namespace MpmSimTest;

int main(int argc, char** argv)
{
	// test shape functions:
	testShapeFunction( CubicBsplineShapeFunction() );
	
	return 0;
}
