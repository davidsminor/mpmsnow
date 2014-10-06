#ifndef MPMSIMTEST_TESTGRID_H
#define MPMSIMTEST_TESTGRID_H

namespace MpmSimTest
{

class TestGrid
{
public:
	static void test();
private:
	static void testProcessingPartitions();
	static void testSplatting();
	static void testDeformationGradients();
	static void testForces();
	static void testImplicitUpdate();
	static void testMovingGrid();
	static void testDfiDxi();
};

}

#endif
