// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


#include "tarch/logging/Log.h"
#include "peano4/grid/GridVertex.h"
#include "tarch/multicore/BooleanSemaphore.h"


#include "TreeEntry.h"


#include <set>
#include <queue>
#include <map>


namespace peano4 {
  namespace parallel {
    class Node;

    namespace tests {
      class NodeTest;
    }
  }
}


/**
 * Node is Peano's abstraction for the hybrid of MPI/threads. It represents
 * one spacetree instance basically.
 *
 * @author Tobias Weinzierl
 */
class peano4::parallel::Node {
  public:
    static constexpr int        Terminate = -2;
    static constexpr int        MaxSpacetreesPerRank = 1024;

    /**
     * Value for _currentProgramStep.
     */
    static constexpr int        UndefProgramStep = -1;
  private:
    friend class peano4::parallel::tests::NodeTest;

    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    static Node                _singleton;

    tarch::multicore::BooleanSemaphore  _semaphore;

    /**
     * We need one stack for outgoing data, one for incoming
     * data, and we distinguish horizontal and vertical data
     * flow.
     */
    static constexpr int        StacksPerCommunicationPartner = 4;

    #if !defined(Parallel)
    typedef int MPI_Comm;
    #endif
    MPI_Comm                    _dataExchangeCommunicators[MaxSpacetreesPerRank];

    int                         _currentProgramStep;

    int                         _rankOrchestrationTag;

    std::map< int, TreeEntry >  _treeEntries;

    /**
     * The standard constructor assigns the attributes default values and
     * checks whether the program is compiled using the -DParallel option.
     * If this is not the case, a warning is logged.
     */
    Node();

    /**
     * The copy constructor is private.
     */
    Node( const Node& node ) = delete;

    /**
     * Peano maps grid instance threads + mpi ranks onto global IDs through
     * this routine. To identify a global stack number, please use this
     * function and add it to your local stack number.
     */
    int getId(int rank, int localTreeId) const;

    /**
     * The operation is not thread-safe as we call it only internally,
     * i.e. you are not supposed to call this function directly.
     *
     * @see getId(int,int)
     */
    void registerId(int id, int masterId);

    /**
     * Analyse to which of the 2d faces it is adjacent
     *
     * We enumerate again in the way we always enumerate, i.e. first the face
     * with the x axis as normal, then the one with the y axis, ... And after
     * these d faces, we enumerate all the opposite faces on the cube.
     *
     * The routine is non-trivial: We store the @f$ 2^d @f$ adjacency flags
     * per vertex, and now have to work backwards which periodic faces they
     * do represent.
     */
    static std::bitset<2*Dimensions> getPeriodicBoundaryNumber(const tarch::la::Vector<TwoPowerD,int>& flags);
  public:
    /**
     * @see getId(int,int)
     */
    int getLocalTreeId(int treeId) const;
    int getGlobalTreeId(int treeId) const;

    /**
     * I originally wanted to embed these guys into the singleton's
     * constructor. However, the singleton might not yet be up when
     * I run the (mpi) tests.
     *
     */
    static void initMPIDatatypes();
    static void shutdownMPIDatatypes();

    /**
     * This operation returns the singleton instance. Before using this
     * instance, one has to call the init() operation on the instance returned.
     *
     * @return The singleton instance
     */
    static Node& getInstance();

    /**
     * The standard destructor calls MPI_Finalize().
     */
    virtual ~Node();

    /**
     * Is this the global master?
     */
    static bool isGlobalMaster(int treeId);

    /**
     * @return Number of registered trees on this rank
     */
    int getNumberOfRegisteredTrees() const;

    /**
     * This operation is not const as it does some internal bookkeeping. It
     * internally invokes  registerId() on the result.
     *
     * @param forTreeId If we run a split, this is the rank that wants to
     *          split. In general, this always is the parent/master of the
     *          new rank. So if you move a rank, you have to ensure that
     *          the master is preserved.
     *
     * @return -1 If there are no ids left anymore. In principle, we do not
     *            limit the number of threads at all. However, we have to
     *            work with a limited number of tags in MPI, so we put a
     *            cap on the number of trees per node.
     * @see getId(int,int)
     */
    int reserveId(int rank, int forTreeId);

    /**
     * You hand in a tree number and the node tells you on which rank such a
     * tree is hosted.
     * @see getId(int,int)
     */
    int getRank(int treeId) const;

    /**
     * Only the SpacetreeSet should call this operation.
     * @see getId(int,int)
     */
    void deregisterId(int id);

    /**
     * Hand in a spacetree id and get back the number that we should use to
     * send something to this tree.
     */
    static int getOutputStackNumberForHorizontalDataExchange(int id);

    /**
     * Get the input stack where a tree writes all of its vertical
     * data to/from when it exchanges information with id. Such vertical
     * data is multiscale information or split/join data.
     */
    static int getInputStackNumberForVerticalDataExchange(int id);
    static int getOutputStackNumberForVerticalDataExchange(int id);

    /**
     * Counterpart of getOutputStackNumberOfBoundaryExchange(int)
     */
    static int getInputStackNumberForHorizontalDataExchange(int id);

    /**
     * A periodic boundary stack is basically a stack (an integer), but I do
     * augment it by a bitset which identifies which symmetry axis belong to
     * this stack number. So if we return (23,010), then data goes to stack 23
     * and it is a stack that does realise periodic boundary conditions along
     * the y axis.
     *
     * The second argument (as sketched above) is of type std::bitset<Dimensions>.
     * Unfortunately, bitsets do not work properly within the pair type as they
     * seem not to define the comparison operator properly. So I use an int
     * here and wrap it into a bitset whenever necessary.
     *
     * @see toString()
     */
    typedef std::pair<int, int > PeriodicBoundaryStackIdentifier;

    /**
     * You hand in a the flags of a vertex and you get the boundary stack identifiers
     * including their directions.
     *
     * ## Algorithm
     *
     * We first analyse which faces with a periodic boundary flag are adjacent to the
     * vertex. For this, we employ getPeriodicBoundaryNumber(). After that, we build
     * the power set over the result: If a face is adjacent to the upper face and the
     * left face, it also has to communicate along the diagonal, i.e. the logical
     * combination of the two.
     *
     * We now have to translate these powerset entries into stack indices. For this,
     * we run over the @f$ 3^d @f$ neighbour cells or a cell, construct the index for
     * this direction, and then look whether this index is contained within the power
     * set.
     *
     * ## Implementation
     *
     * For the powerset, I'd like to use
     *
     *      std::set< std::bitset<2*Dimensions> > powerSet;
     *
     * but C++ defines no comparison operator on bitsets and therefore you can't use
     * it in a set.
     */
    static std::set< PeriodicBoundaryStackIdentifier >  getOutputStacksForPeriodicBoundaryExchange(const tarch::la::Vector<TwoPowerD,int>& flags);

    static std::string toString( const std::set< PeriodicBoundaryStackIdentifier >& data );

    /**
     * Identify output stack for periodic boundary data written by face
     *
     * To make the mapping of face and vertex stacks consistent, we take the
     * face number and determine what the adjacent, i.e. face-connected
     * neighbour volume out of the @f$ 3^d @f$ adjacent, i.e. vertex-connected
     * volumes is.
     *
     * Take care that the result is not a mere linearisation of the index. You
     * have to explicitly exclude the centre volume.
     *
     *
     * @see getPeriodicBoundaryExchangeInputStackNumberForOutputStack() for a
     *   description how to compute an input stack index from the output
     *   index.
     *
     * @see mapPeriodicBoundaryExchangeOutputStackOntoInputStack() for a mapping
     *   of periodic output stacks onto input stacks
     *
     * @return The output stack for a periodic boundary data exchange for one face
     *
     */
    static int  getOutputStackForPeriodicBoundaryExchange(int faceNumber);

    /**
     * Hand in an output stack number of a face and you get back the input
     * stack number.
     *
     * If you want to use it for merging data, you have to be careful: Think about
     * the right face of your domain. It is assigned a BC output stack. This
     * routine tells you onto which input stack this output stack is mapped. In
     * return, if you read in the periodic BC data for the right face, you will
     * merge a data with its own! You however have to merge it with the input BC
     * stack from the left face.
     *
     * This routine is not made to orchestrate the stack exchange! It is meant
     * that you hand in number to which a vertex or face would write to. And then
     * you get the number back which it should read from. This is not how the
     * stacks are permuted after each step! If you want to achieve this, you have
     * to call mapPeriodicBoundaryExchangeOutputStackOntoInputStack().
     *
     * @see peano4::parallel::SpacetreeSet::exchangeAllPeriodicBoundaryDataStacks() for
     *   a routine which uses this one in return.
     *
     * @see peano4::grid::PeanoCurve:MaxNumberOfCoreStacksPerSpacetreeInstance for a
     *   discussion of the ordering of the stacks.
     */
    static int  getPeriodicBoundaryExchangeInputStackNumberForOutputStack(int outputStackNumber);
    static PeriodicBoundaryStackIdentifier  getPeriodicBoundaryExchangeInputStackNumberForOutputStack(PeriodicBoundaryStackIdentifier outputStackIdentifier);

    /**
     * @see tests::NodeTest::testGetOutputStacksForPeriodicBoundaryExchange()
     */
    static int mapPeriodicBoundaryExchangeOutputStackOntoInputStack(int outputStack);

    /**
     * See getOutputStackNumberOfBoundaryExchange().
     */
    static bool isHorizontalDataExchangeOutputStackNumber(int number);
    static bool isHorizontalDataExchangeInputStackNumber(int number);

    static bool isVerticalDataExchangeOutputStackNumber(int number);
    static bool isVerticalDataExchangeInputStackNumber(int number);

    static bool isPeriodicBoundaryExchangeOutputStackNumber(int number);

    /**
     * There are communication stacks and storage stacks. This is
     * a storage stack.
     */
    static bool isStorageStackNumber(int number);

    /**
     * Gives you back the id of a communication partner, i.e. a spacetree
     * number, which corresponds to stack number. This works for all kinds
     * of stacks.
     */
    static int getTreeNumberTiedToExchangeStackNumber(int number);

    /**
     * You should call this operation only on the ranks >0 to find out whether
     * you should do more iteration/sweeps. The spacetree set internally
     * hijacks this operation also on rank 0 to trigger the send out of startup
     * messages. Therefore, the operation is not const.
     */
    bool continueToRun();

    /**
     * The user tells the set which program step to use, i.e. which adapter
     * to use for the next grid run or which step to run next. Should only
     * be called on global master.
     */
    void setNextProgramStep( int number );

    int getCurrentProgramStep() const;

    void init();

    /**
     * The shutdown is invoked by peano4::shutdownSingletons()
     *
     * If we run it on the global master, we set the program state to terminate
     * and distribute this message. For this, we hijack the routine continueToRun().
     * This is really a hack. We in particular have to be careful and may not
     * call the continue operation on the other ranks. They typically are stuck
     * in a continue while loop in their main routine and just wait for the
     * global master to send out this terminate command.
     *
     * It is important that this is the very first call when you want your
     * application to go down. Invoked on the global master, it triggers the
     * shutdown on all the other ranks. You will get deadlocks if you try to free
     * any common memory (RDMA) or free datatypes if you don't call the node's
     * shutdown first.
     *
     * Other ranks might struggle to receive the shutdown command as they might
     * still ask for further subtrees, e.g., or issue load balancing. Therefore,
     * I insert a Peano 4 barrier right after we've issued the shutdown. It makes
     * the rank that goes down handle all incoming data until we can be sure that
     * everybody is in shutdown mode.
     */
    void shutdown();

    enum class ExchangeMode {
      HorizontalData,
      /**
       * Vertical data is data running from the master to the worker and the other
       * way round. Consequently, some tree management data are vertical data, too.
       * See for example peano4::parallel::SpacetreeSet::addSpacetree().
       */
      VerticalData
    };

    static std::string toString( ExchangeMode mode );

    typedef std::pair<int, MPI_Comm>  GridDataExchangeMetaInformation;

    /**
     * I use two tags per spacetree per rank: one for boundary data
     * (horizontal) and one for up-down and synchronous (forks) data
     * exchanges (vertical and fork/join).
     *
     * The nice thing about Peano's new data management is that it is
     * stack-only. Furthermore, all stack data exchange is triggered via the
     * spacetree set, i.e. the order of the stacks is known externally, too.
     * So we effectively do not need that many tags even though we need
     * different tags per tree pair.
     */
    GridDataExchangeMetaInformation getGridDataExchangeMetaInformation( int sendingTreeId, int receivingTreeId, ExchangeMode exchange ) const;

    /**
     * I use this for debugging. When I have a tag, I can ask the node whether
     * the node is aware what this tag actually means.
     */
    static std::string getSemanticsForTag( int tag );
};

