// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace toolbox {
  namespace blockstructured {
    namespace tests {
      class InterpolationTest;
    }
  }
}


class toolbox::blockstructured::tests::InterpolationTest: public tarch::tests::TestCase {
  private:

    /**
     * Tests this averaging
     *
     * The test results from a real test with the Euler equations and four
     * unknowns per coordinate axis. The input is as follows:
     *
<pre>

parameter tarch::la::DynamicMatrix::vectorToString(fineGridFaceEulerQUpdate.value,4*2*5): (entries=40,
{
 0.100047, 0.000933957, 0,0,0.00493399,
 0.100045, 0.000870917, 0,0,0.0045799,
 0.100039, 0.000801622, 0,0,0.00422816,
 0.100045, 0.00081321,  0,0,0.0042177,
 0.100047, 0.000933957, 0,0,0.00493399,
 0.100045, 0.000870917, 2.71051e-20,0,0.0045799,
 0.100039, 0.000801622, 0,0,0.00422816,
 0.100045,0.00081321,0,0,0.0042177
})

</pre>
     *
     */
    void testRestrictHaloLayer_AoS_averaging();

    void testRestrictCellForBreakingDam();

    /**
     * Trivial interpolation test for constant data
     */
    void testInterpolateCellDataAssociatedToVolumesIntoOverlappingCell_linear();

  public:

    /**
     * Cosntructor.
     */
    InterpolationTest();

    /**
     * Destructor, empty.
     */
    virtual ~InterpolationTest() = default;

    /**
     * This routine is triggered by the TestCaseCollection
     */
    virtual void run() override;
};

