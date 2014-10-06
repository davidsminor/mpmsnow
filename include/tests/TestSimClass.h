#ifndef MPMSIMTEST_TESTSIMCLASS_H
#define MPMSIMTEST_TESTSIMCLASS_H

namespace MpmSimTest
{
class TestSimClass
{
public:
	static void test();
private:
	static void testInitialization();
	static void testTimestepAdvance();
};
}

#endif // MPMSIMTEST_TESTSIMCLASS_H
