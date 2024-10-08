#include "tarch/multicore/Core.h"
#include "tarch/multicore/Tasks.h"
#include "tarch/multicore/multicore.h"
#include "tarch/compiler/CompilerSpecificSettings.h"

#ifdef CompilerHasSysinfo
#include <sched.h>
#endif

#if defined(SharedOMP)

#include <omp.h>

#include "config.h"

tarch::logging::Log  tarch::multicore::Core::_log( "tarch::multicore::Core" );


tarch::multicore::Core::Core():
  _numberOfThreads( omp_get_max_threads() ) {

  internal::configureInternalTaskQueues(_numberOfThreads);
}


tarch::multicore::Core::~Core() {}


tarch::multicore::Core& tarch::multicore::Core::getInstance() {
  static Core instance;
  return instance;
}


void tarch::multicore::Core::configure( int numberOfThreads) {
  if ( omp_get_num_procs() != omp_get_max_threads() ) {
    logWarning( "configure(int)", "omp_get_num_procs reports " << omp_get_num_procs() << " while omp_get_max_threads reports " << omp_get_max_threads() << ". Take smaller value to avoid overbooking" );
  }
  int maxThreads = std::min(omp_get_num_procs(), omp_get_max_threads());

  if (numberOfThreads!=UseDefaultNumberOfThreads) {
    if ( maxThreads!=numberOfThreads ) {
      logWarning( "configure(int)", "number of threads configured (" << numberOfThreads << ") does not match system thread level of " << maxThreads << ". OpenMP may ignore manual thread count reset");
    }

    omp_set_num_threads(numberOfThreads);
    _numberOfThreads = omp_get_max_threads();
    logInfo( "configure(int)", "manually reset number of threads used to " << numberOfThreads );
  } else {
    omp_set_num_threads(maxThreads);
    _numberOfThreads = omp_get_max_threads();
  }

  if (_numberOfThreads>getNumberOfUnmaskedThreads()) {
    logWarning( "configure(int)", "number of configured threads (" << _numberOfThreads << ") is bigger than available unmasked threads (" << getNumberOfUnmaskedThreads() << ")" );
    logWarning( "configure(int)", "unmasked threads: " << printUnmaskedThreads() );
  }

  internal::configureInternalTaskQueues(_numberOfThreads);
}


void tarch::multicore::Core::shutdown() {}


bool tarch::multicore::Core::isInitialised() const {
  return true;
}


int tarch::multicore::Core::getNumberOfThreads() const {
  return _numberOfThreads;
}


int tarch::multicore::Core::getThreadNumber() const {
  return omp_get_thread_num();
}


int tarch::multicore::Core::getCoreNumber() const {
  #ifdef CompilerHasSysinfo
    return sched_getcpu();
  #else
    //  https://stackoverflow.com/questions/33745364/sched-getcpu-equivalent-for-os-x
    return 0;
  #endif
}


void tarch::multicore::Core::yield() {
  #pragma omp taskyield
}
#endif
