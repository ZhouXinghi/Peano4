#include "tarch/multicore/Core.h"
#include "tarch/multicore/Tasks.h"
#include "tarch/multicore/multicore.h"
#include "tarch/compiler/CompilerSpecificSettings.h"

#ifdef CompilerHasSysinfo
#include <sched.h>
#endif

#if defined(SharedCPP)
#include <sstream>

#include "config.h"

tarch::logging::Log tarch::multicore::Core::_log( "tarch::multicore::Core" );

tarch::multicore::Core::Core():
  _numberOfThreads(std::thread::hardware_concurrency()) {

  internal::configureInternalTaskQueues(_numberOfThreads);
}

tarch::multicore::Core::~Core() {}

tarch::multicore::Core& tarch::multicore::Core::getInstance() {
  static Core instance;
  return instance;
}

void tarch::multicore::Core::configure( int numberOfThreads ) {
  if (numberOfThreads != UseDefaultNumberOfThreads) {
    _numberOfThreads = numberOfThreads;
  }
  else {
    _numberOfThreads = std::thread::hardware_concurrency();
  }

  if (_numberOfThreads>getNumberOfUnmaskedThreads()) {
    logWarning( "configure(int)", "number of configured threads (" << numberOfThreads << ") is bigger than available unmasked threads (" << getNumberOfUnmaskedThreads() << ")" );
    logWarning( "configure(int)", "unmasked threads: " << printUnmaskedThreads() );
  }

  // we don't change the number of thread queues, as it seems we cannot
  // control the number of threads in C++
}

void tarch::multicore::Core::shutdown() {}

bool tarch::multicore::Core::isInitialised() const {
  return true;
}

int tarch::multicore::Core::getNumberOfThreads() const {
  return _numberOfThreads;
}

int tarch::multicore::Core::getCoreNumber() const {
  #ifdef CompilerHasSysinfo
    return sched_getcpu();
  #else
    //  https://stackoverflow.com/questions/33745364/sched-getcpu-equivalent-for-os-x
    return 1;
  #endif
}

int tarch::multicore::Core::getThreadNumber() const {
  return getCoreNumber();
}


void tarch::multicore::Core::yield() {
  std::this_thread::yield();
}
#endif
