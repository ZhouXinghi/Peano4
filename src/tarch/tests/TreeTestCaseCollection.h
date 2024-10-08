// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

#include <list>
#include <map>

namespace tarch {
  namespace tests {
    class TreeTestCaseCollection;
  }
}


/**
 *
 * @author Tobias Weinzierl
 * @author Wolfgang Eckhardt
 * @version $Revision: 1.7 $
 */
class tarch::tests::TreeTestCaseCollection: public tarch::tests::TestCase {
  private:
    /**
     * Sequence of test cases that are executes on a run() call. The class is
     * not responsible for destroying them.
     */
    std::list<TestCase*>                            _testCases;
    std::map<std::string, TreeTestCaseCollection*>  _subTests;

    /**
     * Tells the collection wheather to log a progression bar.
     */
    bool _writeToLog;

    /**
     * Log interface the class writes to.
     */
    static tarch::logging::Log _log;

    /**
     * Tells whether the testcases contained should be deleted uppon destruction of this object.
     */
    bool _deleteTestCases;

    static bool isNameWithoutHierarchy(const std::string& testCaseName);
    static std::string getFirstIdentifierInHierarchy(const std::string& testCaseName);
    static std::string getRemainingPathWithoutIdentifier(const std::string& testCaseName);
  public:

    /**
     * Creates a test case collection.
     *
     * @param testCaseCollectionName Name of the test case collection. May not
     *   contain ::. If you pass in the empty string "", it is a mere
     *   container. It still runs all sub test cases, but it does not dump any
     *   progress info on this particular package: Usually you get
     *
     *         a.b.c ... ok
     *         a.b ... ok
     *         a ... pass (no local tests)
     *
     *   If a is "", then no message is written.
     */
    TreeTestCaseCollection(
      const std::string& testCaseCollectionName = "",
      bool deleteTestCases = true,
      bool writeToLog = true
    );

    /**
     * Destructor
     */
    virtual ~TreeTestCaseCollection();

    /**
     * Runs all test cases assigned.
     */
    virtual void run();

    /**
     * Runs all test cases assigned.
     */
    virtual void run( const std::string& prefix );

    /**
     * Adds a new test case.
     *
     * Although you pass a pointer, you are still responsible for the instance
     * management.
     *
     * You have to specify in the constructor of a TestCaseCollection whether
     * it should delete testcases contained on construction.
     *
     * @Note: default is set to false (i.e. no destruction)
     */
    void addTestCase( const std::string& fullQualifiedPath, TestCase* testCase );

    /**
     * Same as above, but fully qualified test case name is extracted from passed object.
     */
    void addTestCase( TestCase* testCase );
};


