// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include <string>

#ifdef Parallel
#include <functional>
#include <mpi.h>
#endif

#include "peano4/grid/LoadStoreComputeFlag.h"
#include "peano4/grid/TraversalObserver.h"
#include "peano4/utils/Globals.h"
#include "tarch/la/Vector.h"
#include "tarch/mpi/Rank.h"
#include "tarch/services/ServiceRepository.h"

namespace peano4 {
  namespace grid {
    struct GridVertex;
  } // namespace grid
} // namespace peano4

struct peano4::grid::GridVertex {
public:
  enum class State : int {
    HangingVertex       = 0,
    New                 = 1,
    Unrefined           = 2,
    Refined             = 3,
    RefinementTriggered = 4,
    Refining            = 5,
    EraseTriggered      = 6,
    Erasing             = 7,
    Delete              = 8
  };

  GridVertex() {}
  GridVertex(
    State                                 __state,
    tarch::la::Vector<TwoPowerD, int>     __adjacentRanks,
    tarch::la::Vector<TwoPowerD, int>     __backupOfAdjacentRanks,
    bool                                  __hasBeenAntecessorOfRefinedVertexInPreviousTreeSweep,
    bool                                  __isAntecessorOfRefinedVertexInCurrentTreeSweep,
    bool                                  __hasBeenParentOfSubtreeVertexInPreviousTreeSweep,
    bool                                  __isParentOfSubtreeVertexInCurrentTreeSweep,
    int                                   __numberOfAdjacentRefinedLocalCells,
    tarch::la::Vector<Dimensions, double> __x,
    int                                   __level
  );

  peano4::grid::GridVertex::State   getState() const;
  void                              setState(State value);
  tarch::la::Vector<TwoPowerD, int> getAdjacentRanks() const;
  void                              setAdjacentRanks(const tarch::la::Vector<TwoPowerD, int>& value);
  int                               getAdjacentRanks(int index) const;
  void                              setAdjacentRanks(int index, int value);
  tarch::la::Vector<TwoPowerD, int> getBackupOfAdjacentRanks() const;
  void                              setBackupOfAdjacentRanks(const tarch::la::Vector<TwoPowerD, int>& value);
  int                               getBackupOfAdjacentRanks(int index) const;
  void                              setBackupOfAdjacentRanks(int index, int value);
  bool                              getHasBeenAntecessorOfRefinedVertexInPreviousTreeSweep() const;
  void                              setHasBeenAntecessorOfRefinedVertexInPreviousTreeSweep(bool value);
  bool                              getIsAntecessorOfRefinedVertexInCurrentTreeSweep() const;
  void                              setIsAntecessorOfRefinedVertexInCurrentTreeSweep(bool value);
  bool                              getHasBeenParentOfSubtreeVertexInPreviousTreeSweep() const;
  void                              setHasBeenParentOfSubtreeVertexInPreviousTreeSweep(bool value);
  bool                              getIsParentOfSubtreeVertexInCurrentTreeSweep() const;
  void                              setIsParentOfSubtreeVertexInCurrentTreeSweep(bool value);
  int                               getNumberOfAdjacentRefinedLocalCells() const;
  void                              setNumberOfAdjacentRefinedLocalCells(int value);
#if PeanoDebug > 0
  tarch::la::Vector<Dimensions, double> getX() const;
  void                                  setX(const tarch::la::Vector<Dimensions, double>& value);
  double                                getX(int index) const;
  void                                  setX(int index, double value);
#endif
  int  getLevel() const;
  void setLevel(int value);
  GridVertex(const GridVertex& copy);
  GridVertex& operator=(const GridVertex& other);

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
    int
    getSenderRank() const;

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
  static void send(const peano4::grid::GridVertex& buffer, int destination, int tag, MPI_Comm communicator);
  static void receive(peano4::grid::GridVertex& buffer, int source, int tag, MPI_Comm communicator);

  /**
   * Alternative to the other send() where I trigger a non-blocking send an
   * then invoke the functor until the corresponding MPI_Test tells me that
   * the message went through. In systems with heavy MPI usage, this can
   * help to avoid deadlocks.
   */
  static void send(
    const peano4::grid::GridVertex& buffer,
    int                             destination,
    int                             tag,
    std::function<void()>           startCommunicationFunctor,
    std::function<void()>           waitFunctor,
    MPI_Comm                        communicator
  );
  static void receive(
    peano4::grid::GridVertex& buffer,
    int                       source,
    int                       tag,
    std::function<void()>     startCommunicationFunctor,
    std::function<void()>     waitFunctor,
    MPI_Comm                  communicator
  );
#endif

  enum ObjectConstruction { NoData };

  GridVertex(ObjectConstruction):
  GridVertex() {}

#ifdef Parallel
  static void sendAndPollDanglingMessages(
    const peano4::grid::GridVertex& message,
    int                             destination,
    int                             tag,
    MPI_Comm                        communicator = tarch::mpi::Rank::getInstance().getCommunicator()
  );
  static void receiveAndPollDanglingMessages(
    peano4::grid::GridVertex& message,
    int                       source,
    int                       tag,
    MPI_Comm                  communicator = tarch::mpi::Rank::getInstance().getCommunicator()
  );
#endif

  std::string toString() const;

private:
  /*[[clang::pack]]*/ State _state;
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::pack_range(-1, std::numeric_limits<int>::max())]]*/ int _adjacentRanks[TwoPowerD];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  tarch::la::Vector<TwoPowerD, int> _adjacentRanks;
#endif
  tarch::la::Vector<TwoPowerD, int>       _backupOfAdjacentRanks;
  /*[[clang::pack]]*/ bool                    _hasBeenAntecessorOfRefinedVertexInPreviousTreeSweep;
  /*[[clang::pack]]*/ bool                    _isAntecessorOfRefinedVertexInCurrentTreeSweep;
  /*[[clang::pack]]*/ bool                    _hasBeenParentOfSubtreeVertexInPreviousTreeSweep;
  /*[[clang::pack]]*/ bool                    _isParentOfSubtreeVertexInCurrentTreeSweep;
  /*[[clang::pack_range(0, TwoPowerD)]]*/ int _numberOfAdjacentRefinedLocalCells;
#if PeanoDebug > 0
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::truncate_mantissa(23)]]*/ double _x[Dimensions];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  tarch::la::Vector<Dimensions, double> _x;
#endif
#endif

    int _level;

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
  static MPI_Datatype Datatype;
#endif
#endif
};
