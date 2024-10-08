// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


#include "config.h"
#include "tarch/compiler/CompilerSpecificSettings.h"


#ifdef Parallel
#include <mpi.h>
#endif




/**

 @namespace tarch::mpi

 <h2> Using DaStGen2's send and receive </h2>

 With DaStGen2, we have two variants of sends and receives: We can use the
 blocking variant and a non-blocking one. Furthermore, we can obviously use
 the MPI datatype that the generator yields as well. With the blocking variant,
 the usage is straightforward:

 <pre>
IntegerMessage::receive(message, rank, BarrierTag, tarch::mpi::Rank::getInstance().getCommunicator());
 </pre>

 With the option to insert a wait functor, you can do more complicated stuff.
 The minimalist version is:

 <pre>
tarch::mpi::IntegerMessage::receive(
  message, rank, BarrierTag,
  []() {},
  [&]() {
    tarch::services::ServiceRepository::getInstance().receiveDanglingMessages();
  },
  tarch::mpi::Rank::getInstance().getCommunicator()
);
 </pre>

 It ensures that MPI progress is made as it actively polls the MPI subsystem for
 further incoming messages. I do prefer the more sophisticated version which also
 has a timeout mechanism:

      tarch::mpi::IntegerMessage::receive(
        message, rank, BarrierTag,
        [&]() {
          tarch::mpi::Rank::getInstance().setDeadlockWarningTimeStamp();
          tarch::mpi::Rank::getInstance().setDeadlockTimeOutTimeStamp();
        },
        [&]() {
          tarch::mpi::Rank::getInstance().writeTimeOutWarning( "tarch::mpi::DoubleMessage", "sendAndPollDanglingMessages()",destination, tag );
          tarch::mpi::Rank::getInstance().triggerDeadlockTimeOut( "tarch::mpi::DoubleMessage", "sendAndPollDanglingMessages()", destination, tag );
          tarch::services::ServiceRepository::getInstance().receiveDanglingMessages();
        },
        tarch::mpi::Rank::getInstance().getCommunicator()
      );
 </pre>

 As inserting the code from above is cumbersome, I decided to create my own
 Peano4 aspect and to inject it into DaStGen. So just invoke the sendAndPoll
 or receiveAndPoll operation and you basically get the above aspect.

 */



