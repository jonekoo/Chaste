#ifndef _TESTMAGNETICFIELDCALCULATOR_HPP_
#define _TESTMAGNETICFIELDCALCULATOR_HPP_

#include "MagneticFieldCalculator.hpp"

class TestMagneticFieldCalculator : public CxxTest::TestSuite
{

public:

    void TestThatWillFail()
    {
        TS_ASSERT(1==0);
    }
};


#endif /* _TESTMAGNETICFIELDCALCULATOR_HPP_ */