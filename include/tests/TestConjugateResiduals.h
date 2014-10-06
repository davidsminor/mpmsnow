#ifndef MPMSIMTEST_TESTCONJUGATERESIDUALS_H
#define MPMSIMTEST_TESTCONJUGATERESIDUALS_H

namespace MpmSimTest
{

class TestConjugateResiduals
{
public:
	static void test();
private:
	static void testSolve();
	static void testPreconditioner();
};

}
#endif
