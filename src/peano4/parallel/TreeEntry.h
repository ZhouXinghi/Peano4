// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include <string>

#ifdef Parallel
  #include <mpi.h>
  #include <functional>
#endif

#include "tarch/la/Vector.h"
#include "tarch/mpi/Rank.h"
#include "tarch/services/ServiceRepository.h"
#include "peano4/grid/LoadStoreComputeFlag.h"
#include "peano4/utils/Globals.h"
#include "peano4/grid/TraversalObserver.h"


  
namespace peano4{
namespace parallel{
  struct TreeEntry;
}
}

struct peano4::parallel::TreeEntry {
  public:
    TreeEntry() {}
    TreeEntry(int  __id, int  __master);

    int   getId() const;
    void   setId(int value);
    int   getMaster() const;
    void   setMaster(int value);
    TreeEntry(const TreeEntry& copy) = default;


    #ifdef Parallel
    /**
     * Hands out MPI datatype if we work without the LLVM MPI extension.
     * If we work with this additional feature, this is the routine where
     * the lazy initialisation is done and the datatype is also cached.
     */
    /*[[clang::map_mpi_datatype]]*/
    static MPI_Datatype  getForkDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static MPI_Datatype  getJoinDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static MPI_Datatype  getBoundaryExchangeDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static MPI_Datatype  getMultiscaleDataExchangeDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static MPI_Datatype  getGlobalCommunciationDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static void  freeForkDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static void  freeJoinDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static void  freeBoundaryExchangeDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static void  freeMultiscaleDataExchangeDatatype();

    /*[[clang::map_mpi_datatype]]*/
    static void  freeGlobalCommunciationDatatype();

    /**
     * @return The rank of the sender of an object. It only make ssense to call
     *         this routine after you've invoked receive with MPI_ANY_SOURCE.
     */
    int getSenderRank() const;

    /**
     * Wrapper around getDatatype() to trigger lazy evaluation if we
     * use the lazy initialisation.
     */
    static void initDatatype();

    /**
     * Free the underlying MPI datatype.
     */
    static void shutdownDatatype();

    /**
     * In DaStGen (the first version), I had a non-static version of the send
     * as well as the receive. However, this did not work with newer C++11 
     * versions, as a member function using this as pointer usually doesn't 
     * see the vtable while the init sees the object from outside, i.e. 
     * including a vtable. So this routine now is basically an alias for a
     * blocking MPI_Send. 
     */
    static void send(const peano4::parallel::TreeEntry& buffer, int destination, int tag, MPI_Comm communicator );
    static void receive(peano4::parallel::TreeEntry& buffer, int source, int tag, MPI_Comm communicator );

    /**
     * Alternative to the other send() where I trigger a non-blocking send an 
     * then invoke the functor until the corresponding MPI_Test tells me that 
     * the message went through. In systems with heavy MPI usage, this can 
     * help to avoid deadlocks.
     */
    static void send(const peano4::parallel::TreeEntry& buffer, int destination, int tag, std::function<void()> startCommunicationFunctor, std::function<void()> waitFunctor, MPI_Comm communicator );
    static void receive(   peano4::parallel::TreeEntry& buffer, int source,      int tag, std::function<void()> startCommunicationFunctor, std::function<void()> waitFunctor, MPI_Comm communicator );
    #endif


    enum ObjectConstruction {
      NoData
    };

    TreeEntry( ObjectConstruction ):
    TreeEntry() {}

#ifdef Parallel
    static void sendAndPollDanglingMessages(const peano4::parallel::TreeEntry& message, int destination, int tag, MPI_Comm communicator=tarch::mpi::Rank::getInstance().getCommunicator());
    static void receiveAndPollDanglingMessages(peano4::parallel::TreeEntry& message, int source, int tag, MPI_Comm communicator=tarch::mpi::Rank::getInstance().getCommunicator() );
#endif


    std::string toString() const;

  private:
    /*[[clang::pack_range(0,std::numeric_limits<int>::max())]]*/  int   _id;
    /*[[clang::pack_range(-1,std::numeric_limits<int>::max())]]*/  int   _master;

    #ifdef Parallel
    private:
      int _senderDestinationRank;

      #if !defined(__MPI_ATTRIBUTES_LANGUAGE_EXTENSION__)
      /**
       * Whenever we use LLVM's MPI extension (DaStGe), we rely on lazy
       * initialisation of the datatype. However, Peano calls init explicitly
       * in most cases. Without the LLVM extension which caches the MPI 
       * datatype once constructed, this field stores the type.
       */
      static MPI_Datatype  Datatype;
      #endif
    #endif
};
