#include "tarch/tests/TestCaseCollection.h"


tarch::logging::Log tarch::tests::TestCaseCollection::_log("tarch::tests::TestCaseCollection");


tarch::tests::TestCaseCollection::TestCaseCollection():
  TestCase( "" ) {
}


tarch::tests::TestCaseCollection::TestCaseCollection(const std::string& testCaseCollectionName, bool deleteTestCases, bool writeToLog):
  TestCase::TestCase( testCaseCollectionName ),
  _writeToLog(writeToLog),
  _deleteTestCases(deleteTestCases) {
}


tarch::tests::TestCaseCollection::~TestCaseCollection() {
  if (_deleteTestCases) {
    for (std::list<tarch::tests::TestCase*>::iterator p = _testCases.begin(); p!=_testCases.end(); p++ ) {
      tarch::tests::TestCase* currentTestCase = *p;
      delete currentTestCase;
    }
  }
}


void tarch::tests::TestCaseCollection::run() {
  logTraceInWith1Argument( "run()", _testCaseName );
  std::string logInformation = "running test case collection \"" + _testCaseName + "\" ";
  for (std::list<tarch::tests::TestCase*>::iterator p = _testCases.begin(); p!=_testCases.end(); p++ ) {
    tarch::tests::TestCase* currentTestCase = *p;
    logInformation += ".";
    currentTestCase->run();
    int additionalErrors = currentTestCase->getNumberOfErrors();
    if (additionalErrors>0) {
      logInformation += "x";
      _errors += currentTestCase->getNumberOfErrors();
    }
    else {
      logInformation += ".";
    }
  }
  if (_errors==0) {
    logInformation += " ok";
  }
  else {
    logInformation += " failed";
  }
  if (_writeToLog) {
    logInfo("run()",logInformation );
  }
  logTraceOutWith2Arguments( "run()", _testCaseName, _errors );
}


void tarch::tests::TestCaseCollection::addTestCase( TestCase* testCase ) {
  _testCases.push_back(testCase);
}
