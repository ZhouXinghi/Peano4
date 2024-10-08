
//
// Generated by DaStGen2 (C) 2020 Tobias Weinzierl
//
// For DaStGen's copyright, visit www.peano-framework.org. These generated files
// however are not subject of copyright, i.e. feel free to add your copyright in
// here
//
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
    struct AutomatonState;
  } // namespace grid
} // namespace peano4

struct peano4::grid::AutomatonState {
public:
  AutomatonState() = default;
  AutomatonState(
    int                                        __level,
    tarch::la::Vector<Dimensions, double>      __x,
    tarch::la::Vector<Dimensions, double>      __h,
    bool                                       __inverted,
    std::bitset<Dimensions>                    __evenFlags,
    tarch::la::Vector<DimensionsTimesTwo, int> __accessNumber
  );

  AutomatonState(const AutomatonState& copy);
  AutomatonState& operator=(const AutomatonState& other);

  ~AutomatonState() = default;

  int                                        getLevel() const;
  void                                       setLevel(int value);
  tarch::la::Vector<Dimensions, double>      getX() const;
  void                                       setX(const tarch::la::Vector<Dimensions, double>& value);
  double                                     getX(int index) const;
  void                                       setX(int index, double value);
  tarch::la::Vector<Dimensions, double>      getH() const;
  void                                       setH(const tarch::la::Vector<Dimensions, double>& value);
  double                                     getH(int index) const;
  void                                       setH(int index, double value);
  bool                                       getInverted() const;
  void                                       setInverted(bool value);
  std::bitset<Dimensions>                    getEvenFlags() const;
  void                                       setEvenFlags(const std::bitset<Dimensions>& value);
  bool                                       getEvenFlags(int index) const;
  void                                       setEvenFlags(int index, bool value);
  void                                       flipEvenFlags(int index);
  tarch::la::Vector<DimensionsTimesTwo, int> getAccessNumber() const;
  void                                       setAccessNumber(const tarch::la::Vector<DimensionsTimesTwo, int>& value);
  int                                        getAccessNumber(int index) const;
  void                                       setAccessNumber(int index, int value);

#ifdef Parallel
    /**
     * Hands out MPI datatype if we work without the LLVM MPI extension.
     * If we work with this additional feature, this is the routine where
     * the lazy initialisation is done and the datatype is also cached.
     */
    //[[clang::map_mpi_datatype]]
    static MPI_Datatype  getForkDatatype();

    //[[clang::map_mpi_datatype]]
    static MPI_Datatype  getJoinDatatype();

    //[[clang::map_mpi_datatype]]
    static MPI_Datatype  getBoundaryExchangeDatatype();

    //[[clang::map_mpi_datatype]]
    static MPI_Datatype  getMultiscaleDataExchangeDatatype();

    //[[clang::map_mpi_datatype]]
    static MPI_Datatype  getGlobalCommunciationDatatype();

    //[[clang::map_mpi_datatype]]
    static void  freeForkDatatype();

    //[[clang::map_mpi_datatype]]
    static void  freeJoinDatatype();

    //[[clang::map_mpi_datatype]]
    static void  freeBoundaryExchangeDatatype();

    //[[clang::map_mpi_datatype]]
    static void  freeMultiscaleDataExchangeDatatype();

    //[[clang::map_mpi_datatype]]
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
  static void send(const peano4::grid::AutomatonState& buffer, int destination, int tag, MPI_Comm communicator);
  static void receive(peano4::grid::AutomatonState& buffer, int source, int tag, MPI_Comm communicator);

  /**
   * Alternative to the other send() where I trigger a non-blocking send an
   * then invoke the functor until the corresponding MPI_Test tells me that
   * the message went through. In systems with heavy MPI usage, this can
   * help to avoid deadlocks.
   */
  static void send(
    const peano4::grid::AutomatonState& buffer,
    int                                 destination,
    int                                 tag,
    std::function<void()>               startCommunicationFunctor,
    std::function<void()>               waitFunctor,
    MPI_Comm                            communicator
  );
  static void receive(
    peano4::grid::AutomatonState& buffer,
    int                           source,
    int                           tag,
    std::function<void()>         startCommunicationFunctor,
    std::function<void()>         waitFunctor,
    MPI_Comm                      communicator
  );
#endif

  enum ObjectConstruction { NoData };

  AutomatonState(ObjectConstruction):
  AutomatonState() {}

#ifdef Parallel
  static void sendAndPollDanglingMessages(
    const peano4::grid::AutomatonState& message,
    int                                 destination,
    int                                 tag,
    MPI_Comm                            communicator = tarch::mpi::Rank::getInstance().getCommunicator()
  );
  static void receiveAndPollDanglingMessages(
    peano4::grid::AutomatonState& message,
    int                           source,
    int                           tag,
    MPI_Comm                      communicator = tarch::mpi::Rank::getInstance().getCommunicator()
  );
#endif


  std::string toString() const;

private:
  /*[[clang::pack_range(0, 63)]]*/ int _level;
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::truncate_mantissa(23)]]*/ double _x[Dimensions];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  tarch::la::Vector<Dimensions, double> _x;
#endif
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::truncate_mantissa(23)]]*/ double _h[Dimensions];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  tarch::la::Vector<Dimensions, double> _h;
#endif
  /*[[clang::pack]]*/ bool _inverted;
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::pack]]*/ bool _evenFlags[Dimensions];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  std::bitset<Dimensions> _evenFlags;
#endif
#if defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  /*[[clang::pack_range(-TwoTimesD, TwoTimesD)]]*/ int _accessNumber[DimensionsTimesTwo];
#endif
#if !defined(__PACKED_ATTRIBUTES_LANGUAGE_EXTENSION__)
  tarch::la::Vector<DimensionsTimesTwo, int> _accessNumber;
#endif

#ifdef Parallel
    private: int _senderDestinationRank;

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
