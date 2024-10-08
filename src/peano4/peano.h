// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "peano4/utils/Globals.h"
#include "tarch/la/Vector.h"


/**
 * @dir src/peano4
 *
 * Peano's core layer. It sits on top of the technical architecture.
 */
namespace peano4 {
  /**
   * You can invoke this operation manually, but it will also implicitly be
   * triggered by the init routines.
   */
  void writeCopyrightMessage();

  /**
   * Fill Lookup Tables
   *
   * Fill all the lookup tables used within the application. As lookup
   * tables are used by many operations, I suggest to call this operation
   * as soon as possible.
   *
   * There shall no error occur in this operation. Thus, it does not return
   * any code.
   */
  void fillLookupTables();

  /**
   * Init Parallel Environment
   *
   * Inits the parallel environment. If the parallel mode is not set, the
   * operation detoriates to nop. The function returns 0 if everything
   * is o.k., it returns -2 otherwise. Please call this operation before
   * you call any other operation that could result in an error. I suggest
   * to call it right after fillLookupTables().
   *
   * Please note that Peano 4 considers both shared memory and distributed
   * memory to be a parallel environment.
   *
   * init might change the variables passed. If you want to parse the
   * command line arguments, use the values returned. If you use the
   * arguments directly without calling initParallelEnvironment() they
   * might contain MPI values not important to the program.
   *
   * <h2> Usage/implementation details </h2>
   *
   * You may not use the trace macros before this operation has invoked the init
   * operation. Otherwise, the getRank() assertion fails, as the node has not
   * been configured correctly.
   *
   * Invoke with an address operator before that.
   * <pre>
  peano4::initParallelEnvironment(&argc,&argv);
     </pre>
   * This has to be done as one of the very first things, i.e. before you init
   * the logging, or run tests, or ...
   */
  int initParallelEnvironment(int* argc, char*** argv);

  /**
   * Shutdown all the parallel environment, i.e. free all MPI datatypes and
   * close down MPI. This also turns off the shared memory environment.
   * Before this happens, you have to shutdown the node such that everybody
   * knows that we are going down. So you have to call Node::shutdown()
   * before you trigger this operation. This is your responsibility.
   *
   * The routine first adds a barrier. This barrier is necessary. If the very last
   * activity of all ranks is for example to plot stuff, they typically use
   * global semaphores as well. To make these semaphores work, we still require
   * that all nodes call receiveDanglingMessages(). It is only after everyone
   * has done their dump, that we can shut down the shared memory system. This is the
   * reason the barrier has to come after the node's shutdown and then has to
   * be a Peano 4 barrier which still invokes receiveDanglingMessages() on all
   * services.
   *
   * Once all shared memory tasks have terminated, we free the MPI datatypes.
   *
   * Eventually, we shut down the MPI rank.
   *
   * Once this routine has terminated, do not add any barrier() anymore!
   *
   * @see peano4::parallel::Node::shutdown()
   */
  void shutdownParallelEnvironment();

  /**
   * Fire up all the singletons.
   *
   *
   * Singletons that I don't touch are:
   *
   * - tarch::mpi::Rank, as the rank is not a service and is handled
   *   separately by the initParallelEnvironment().
   *
   * Note that there is a bug when passing in initializer lists
   * to this function. If the number of elements in "offset" is
   * greater than Dimensions, then "width" will borrow some elements that
   * were intended for "offset". Best practice is to initialize
   * vectors that are to be passed in before the call to initSingletons.
   * That way the compiler will throw an error if the sizes are not correct.
   */
  void initSingletons(
    const tarch::la::Vector<Dimensions, double>& offset,
    const tarch::la::Vector<Dimensions, double>& width,
    const std::bitset<Dimensions>&               periodicBC = 0
  );

  /**
   * The very first thing I have to do is to shut down Node. This
   * shutdown will tell all the other ranks to go down as well.
   */
  void shutdownSingletons();
} // namespace peano4
