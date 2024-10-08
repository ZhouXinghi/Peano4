// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    namespace tests {
      class MatrixTest;
    }
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::tests::MatrixTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testConstruction ();

  /**
   * Tests methods from MatrixAssign.h.
   */
  void testAssignment ();

  /**
   * Tests methods from MatrixOperations.h.
   */
  void testMatrixOperations ();

  /**
   * Tests methods from MatrixMatrixOperations.h.
   */
  void testMatrixMatrixOperations ();

  /**
   * Tests operations with TransposedMatrix.h.
   */
  void testTransposedMatrix ();

public:

  /**
   * Cosntructor.
   */
  MatrixTest ();

  /**
   * Destructor, empty.
   */
  virtual ~MatrixTest () {};

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

