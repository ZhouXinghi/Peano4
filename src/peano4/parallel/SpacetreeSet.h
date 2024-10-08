// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "tarch/logging/Log.h"

#include "tarch/services/Service.h"

#include "tarch/multicore/Tasks.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include "peano4/maps/maps.h"
#include "peano4/grid/grid.h"

#include "TreeManagementMessage.h"

#include <list>
#include <map>
#include <set>

namespace peano4 {
  namespace grid {
    class Spacetree;
    class TraversalObserver;
    struct GridStatistics;
    struct AutomatonState;
  }

  namespace parallel {
    class SpacetreeSet;
  }
}

/**
 * The spacetree set has to be a singleton, as it is reponsible to accept
 * requests for new trees from remote ranks. So there is only one set.
 *
 * @author Tobias Weinzierl
 */
class peano4::parallel::SpacetreeSet: public tarch::services::Service {
  private:
    friend class peano4::grid::Spacetree;

    static SpacetreeSet  _singleton;

    std::vector<peano4::parallel::TreeManagementMessage>   _unansweredMessages;

    /**
     * Each task triggers the traversal of one specific spacetree. After
     * that, we might directly trigger the data exchanges. Yet, this is not a
     * good idea as other tasks might linger in the background not have sent
     * the data out yet. So we don't to anything here.
     */
    class TraverseTask: public tarch::multicore::Task {
      private:
        peano4::grid::Spacetree&          _spacetree;
        SpacetreeSet&                     _spacetreeSet;
        peano4::grid::TraversalObserver&  _observer;
        const bool                        _invertTreeTraversalDirectionBeforeWeStart;

      public:
        TraverseTask( peano4::grid::Spacetree&  tree, SpacetreeSet& set, peano4::grid::TraversalObserver&  observer, bool invertTreeTraversalDirectionBeforeWeStart );

  	    /**
         * I create the copy of the observer, run the traversal on my local
         * tree _spacetree and finally destroy the local observer copy.
         */
        virtual bool run() override;
    };

  public:
    /**
     * Realise domain boundary exchange (over multiple scales)
     *
     * This routine realises the boundary data exchange, i.e. it
     * send outs in iteration n and receives (see counterpart routine), but
     * the received data is used only in traversal n+1. The whole exchange
     * thus asynchronous (logically).
     *
     * Typically, this routine is called by user observers to for their
     * user-defined data. We however also use it within the grid to trigger
     * the actual vertex data exchange. As these routines all are invoked
     * via tasks, this routine is typically invoked by multiple threads
     * concurrently.
     *
     *
     * ## Deep copy semantics for shared memory
     *
     * Whenever we copy a container over to another thread, we use the
     * clone(). When we study the memory ownership, it is important to
     * memorise that the boundary exchange containers are copies of "real"
     * containers of the code. That is, the mesh looks at a face or vertex
     * and stores its content persistently within the tree. Once it has
     * done this, it analyses if this face or vertex is adjacent to a domain
     * boundary. If so, it copies the data of interest into the boundary
     * container. After the traversal, this copy is either sent out to
     * another rank, or we clone() it into the container of another tree
     * on the same rank (shared memory data exchange). After that, the
     * original boundary container is cleared.
     *
     *
     * If you work with containers over pointers, the clone()
     * should create a deep copy, whereas otherwise any push and pull only
     * moves pointers around. That is, a container holding boundary data
     * basically holds pointers to data that is also held by containers
     * administering interior data. After the traversal, it is the clone()
     * which ensures that we actually create copies, i.e. work with
     * real halo data which is replicated if data belongs to another tree -
     * even if that tree resides on the same rank. Therefore, we can clear()
     * the output stack. All data there had been pointers to data which are
     * still held by the local tree.
     *
     *
     * It is the responsibility of the actual merge to ensure that data on the
     * input containers are probably merged into the local data. If these data
     * are pointers, it is the responsibility of the merge to ensure that they
     * are properly deleted. The input stack will be empty after the traversal,
     * as the grid management will pull all data and pipe it into merges, but
     * if they hold pointers, you have to take special care.
     *
     * @see python4/peano4/toolbox/particles/ParticleSet.template.h and the
     *   merge() implementation therein for an example of a merge of a set of
     *   pointers.
     *
     *
     * ## Copy semantics for distributed memory
     *
     * If you launch a distributed memory exchange, it is the responsibility
     * of the container to create a real copy on the incoming rank. If your
     * container managed pointers, this is obviously a deep copy. It is also
     * the container's responsibility to ensure that all data are freed when
     * we call clear() after the send is complete.
     *
     * If your container holds pointers, it holds pointers which are also held
     * by some containers that are responsible for the interior. So the clear()
     * is just fine, as it eliminates redundant pointers.
     *
     *
     * @param symmetricDataCardinality If two ranks A and B are adjacent to
     *   each other, and if A send n entries to B, then B received exactly
     *   n entries in turn.
     *
     * @see finishAllOutstandingSendsAndReceives();
     */
    template <class Container>
    static void exchangeAllHorizontalDataExchangeStacks( Container& stackContainer, int spacetreeId, bool symmetricDataCardinality );

    /**
     * Exchange periodic BC data
     *
     * Close to trivial routine, as it basically examines whether anything
     * has been written to a stack that is associated with periodic BCs.
     * If so, the corresponding stack is copied over onto the respective
     * input stack for the subsequent iteration. This routine should
     * degenerate to nop on all ranks besides rank 0.
     */
    template <class Container>
    static void exchangeAllPeriodicBoundaryDataStacks( Container& stackContainer, int spacetreeId );

    /**
     * <h2> Prepare new children after fork </h2>
     *
     * When we prepare new children, we have to map output stack on output stack, as both
     * trees are not yet inverted at the time we invoke this routine. That is, the master's
     * output stack is mapped on to the new child's outpu stack. The master is inverted
     * automatically later on (that's the last thing we do throughout the traversal). For
     * the new child, we have to invert manually however. This is done in
     * exchangeDataBetweenExistingAndNewTreesAndRerunNewTrees().
     *
     * @param childrenIds Depending on the context, this is either the children or
     *   the new children that are just about to be kicked off
     */
    template <class Container>
    static void exchangeAllVerticalDataExchangeStacks(
      Container& stackContainer,
      int spacetreeId, int parentId
    );

    template <class Container>
    static void deleteAllStacks(
      Container& stackContainer,
      int spacetreeId
    );

    /**
     * This routine finishes all the sends and receives that are still active,
     * i.e. it searched for pending MPI requests and waits for them to finish.
     * After this is ensured, the routine runs over all stacks and ensures that
     * all temporary data is released. The last step is important, as we otherwise
     * quickly run out of memory - we replicate all data structures whenever we fork
     * and as C++ vectors don't release memory, this memory would be lost without
     * a manual freeing.
     */
    template <class Container>
    static void finishAllOutstandingSendsAndReceives( Container& stackContainer, int spacetreeId );

    /**
     * Copies (streams) data from the master to the worker.
     *
     * These are complete copies of the data set as we know that both trees afterwards
     * will coarsen their mesh and thus free resources.
     *
     * Invoked by user code from streamDataFromSplittingTreeToNewTree() and by
     * spacetree set through SpacetreeSet::streamDataFromSplittingTreesToNewTrees().
     * SpacetreeSet runs over the set of trees with the label EmptyRun, i.e. those
     * that haven't done any iteration yet. It then invokes this routine (indirectly)
     * on the master. That is, even if you are on a rank where the master does not exist,
     * the code will temporarily create an observer for the master and then ask this
     * observer to trigger the data exchange.
     *
     * We have to trigger the routine multiple times, to ensure we catch both the
     * on-rank and the between-rank case. As a consequence, we may only stream if the
     * destination stack is still empty.
     *
     * All data that are streamed are clones of the original data. These clones
     * are created by peano4::parallel::SpacetreeSet::streamDataFromSplittingTreesToNewTrees()
     * before we trigger this routine.
     *
     * ## On-rank realisation
     *
     * If source and destination rank are the same, a tree splits up into two trees
     * both residing on the same MPI rank. We therefore simply copy the stream. As
     * this routine is invoked on the master, it is the master that creates the
     * stream on the worker and befills it.
     *
     * Please note that in this particular case, the routine is called twice: We call
     * it per rank for the splitting trees and then for the empty ones. The second
     * call however degenerates to nop.
     *
     * Where to the clones happen?
     *
     * ## Data exchange between different ranks
     *
     * If the master rank is the local guy, then we have to trigger a send. Otherwise,
     * we trigger a receive. The message exchange consists of two phases. An integer
     * message first is exchanged. It carries the number of messages. After that, I
     * send out the actual data.
     *
     * We don't have to finish any sends, i.e. wait for Isends or Irecvs. SpacetreeSet
     * will call finishAllOutstandingSendsAndReceives() later on.
     *
     * The routine is idempotent on a single rank, i.e. you can call it multiple times.
     * Only the first one will copy, all the others will become nop. It is not idempotent
     * in a parallel sense. It has to be idempotent, as I indeed have to call it twice
     * in a distributed memory environment: I have to call it on the receiver side and
     * on the sender side.
     *
     * <h2> Data ownership </h2>
     *
     * The routine solely accesses vertical/stream stacks, i.e. checks
     * getOutputStackNumberForVerticalDataExchange(). That is, we assume that all data
     * are ready to be streamed out. For user data, the spacetree's
     * peano4::grid::Spacetree::sendUserData() will deposit all outgoing user data in
     * the respective output buffers. That is, by the time we hit
     * streamDataFromSplittingTreeToNewTree(), all user data is already in a (temporary)
     * output stream. For the actual tree data, we have to deposit it there manually.
     * This happens in
     * peano4::parallel::SpacetreeSet::streamDataFromSplittingTreesToNewTrees().
     *
     *
     * If we send out data via MPI, I assume that the data resides within a bespoke
     * buffer. tryToFinishSendOrReceive() after all clears the buffer once it has gone
     * out. So it is kind of safe to deposit data in a temporary buffer once. The send
     * implicitly will clear it. After that, the garbage collection will eventually
     * really free the underlying memory.
     *
     * If the target tree and the destination tree reside on the same rank, this routine
     * does a plain copy-over. In this case, we have to delete the (temporary) stack that
     * is used to stream out stuff after the copy. Otherwise, it will be dead memory in
     * the code. Even worse, should the source tree later on stream once again to the
     * same tree - as the tree meanwhile has degenerated or fused back into the master
     * code - then we get a memory corruption as you can't clone into an existing stack.
     *
     * In peano4::parallel::SpacetreeSet::streamDataFromSplittingTreesToNewTrees() we
     * might copy data multiple times if we immediately erase it after the clone.
     * Therefore, I check per call (the function is called twice) whether the target
     * buffer is empty. If so, I copy over. In any case, teh source buffer is deleted
     * if it is not empty - even though I have not cloned it in the second step. This
     * is kind of a hack, but no better one comes to my mind. The whole discussion is
     * not relevant for inter-rank communication as we know that the buffer won't
     * be deleted before the data has not been sent out - which is only once for the
     * whole cycle.
     */
    template <class Container>
    static void streamDataFromSplittingTreeToNewTree( Container& stackContainer, int master, int worker );

    template <class Container>
    static void streamDataFromJoiningTreeToMasterTree( Container& stackContainer, int master, int worker );

  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    /**
     * Semaphore to protect container holding all the local trees.
     */
    static tarch::multicore::BooleanSemaphore                      _semaphore;

    enum class SpacetreeSetState {
      Waiting,
      TraverseTreesAndExchangeData
    };

    static std::string toString( SpacetreeSetState state );

    /**
     * I use this tag to identify messages send from one tree to another rank.
     * All answers go through an answer tag. To identify the right one, please
     * use getAnswerTag(). The tag should be const, but I set it in init() and
     * therefore I had to remove the const - even though its logically
     * not possible to change it.
     */
    int     _requestMessageTag;

    /**
     * Never use this tag directly. It is the first tag of a series fo answer
     * tags. To find the right one for a particular application context, use
     * getAnswerTag().
     */
    int     _answerMessageTag;

    /**
     * These are the local spacetrees.
     */
    std::list< peano4::grid::Spacetree >  _spacetrees;

    /**
     * The state identifies what the set is doing right now. The flag is for example
     * used by receiveDanglingMessages() as we may not create new trees or change
     * tree states while we are running through them or realise their data exchange.
     */
    SpacetreeSetState                     _state;

    /**
     * I create/clone one observer per local tree. This is an
     * on-the-fly process. At the end of the set traversal, we delete all of
     * clones. Originally, I did this delete right after the creation and the
     * subsequent traversal. However, we need the observer clone once more for
     * the data exchange. To avoid reclones, I store all clones in this map and
     * then delete them en bloc.
     *
     * @see deleteClonedObservers
     */
    std::map< int, peano4::grid::TraversalObserver* >    _clonedObserver;

    peano4::grid::Spacetree& getSpacetree(int id);

    const peano4::grid::Spacetree& getSpacetree(int id) const;

    /**
     * @return tag that one should use to answer one particular spacetree
     */
    int getAnswerTag( int targetSpacetreeId ) const;

    /**
     * <h2> Multithreading </h2>
     *
     * I originally intended to make this routine use tasks. However, the data
     * transfer routines do modify the underlying map/stack structures of Peano,
     * as new containers are created, data is moved around, and so forth. So the
     * data exchange per se is not thread-safe. As we do not/can not use threads,
     * we have to be careful with the order. I originally had a loop over the
     * trees that triggered per tree the data exchange and then waited for MPI.
     * Obviously, this does not work with rendezvous protocol. We need two loops.
     * The first one triggers all the MPI sends/receives and it also realises the
     * local data exchange. The second one waits for all MPI operations to
     * terminate.
     */
    void exchangeHorizontalDataBetweenTrees(peano4::grid::TraversalObserver&  observer);
    void exchangeVerticalDataBetweenTrees(peano4::grid::TraversalObserver&  observer);

    /**
     * I do this after a join/after I've removed an empty tree. Have to call it
     * explicitly, as a join does not delete/throw away the data. It simply hands
     * on data ownership.
     */
    void deleteAllStacks( peano4::grid::TraversalObserver&  observer, int spacetreeId );

    /**
     * Copy the data from a splitting tree onto its new workers
     *
     * When we split a tree, we realise this split in two grid sweeps where the second
     * sweep breaks up the traversal into three logical substeps. In the first sweep, the
     * splitting master tells everybody around that it will split. No split is done though.
     * After this first sweep, the grid data structure is replicated on the new worker.
     * Please note that we only replicate the grid, i.e. there's no user information
     * associated with it yet.
     * 
     * In the second sweep, the master still "owns" all data, i.e. it receives all
     * boundary data and merges it in. However, it does not send out boundary data anymore.
     * Instead, it takes all data that should go over to the new worker and dumps it in the
     * vertical data exchange stack. That's part one of the second sweep. Part two means
     * that the new worker runs through its mesh for the first time and incorporates all
     * the streamed user data. After that, we have to run through the grid one more time
     * on the worker (without any action) to get all the data into the right order.
     *
     * ## Vertex grid data
     *
     * The copying of the whole tree data is literally a copying of the output stack of
     * the master. It happens after the master has finished its splitting traversal.
     * peano4::grid::Spacetree::traverse() does invert the traversal direction in an
     * epilogue automatically. Therefore, by the time we hit this routine, we have to
     * copy over the input stack - it is already in the right order for the "counter"-
     * traversal.
     *
     * We always have to invoke the data exchange for both the master and the worker, i.e.
     * we call it twice. This way, we invoke the data exchange on the destination rank
     * (MPI receive) and on the source rank (MPI send). For a single-node split, the
     * second invocation degenerates to nop automatically. See streamDataFromSplittingTreeToNewTree()
     * which implements a simple emptyness check.
     *
     * ## User data
     *
     * It is only the vertex grid data that is copied over in one rush prior to the splitting
     * state. User data is streamed within the splitting state. See Spacetree::sendUserData().
     * The latter routine dumps data in the vertical data exchange stacks, and these stacks
     * subsequently are then transferred over by the present routine, too.
     *
     * ## Stack types
     *
     * All the routines copy over the stream stack data, i.e. the stacks for the vertical
     * data exchange. It is hence important to recognise that the normal stacks are not
     * transferred.
     *
     * @see streamDataFromSplittingTreeToNewTree()
     * @see peano4::grid::Spacetree::sendUserData()
     */
    void streamDataFromSplittingTreesToNewTrees(peano4::grid::TraversalObserver&  observer);

    /**
     * This operation should be called pretty close towards the end of a traversal.
     * I recommend to invoke it after cleanUpTrees(). It creates new trees and thus
     * might invoke remote method calls on other ranks. It is thus important that
     * all the data exchange locally is done - otherwise we might ask another rank
     * to create a tree for us, while the other rank still waits for boundary data
     * from our side.
     *
     * @see exchangeDataBetweenTrees() for details.
     */
    void createNewTrees();

    /**
     * @see _clonedObserver
     */
    void deleteClonedObservers();

    /**
     * Adds a new spacetree to the set. The responsibility goes over to the
     * set. The operation clones the original spacetree handed in into a new
     * spacetree with the id newTreeId.
     *
     * <h2> Local node </h2>
     *
     * If the new tree will be a local tree, we simply add a new tree object.
     * The Spacetree's constructor takes care of all the "cloning". Clone here
     * means a setup. We traverse this tree afterwards separately to stream all
     * data in.
     *
     * <h2> Distributed memory <h2>
     *
     * In a distributed memory environment, we have to break up the creation
     * into a set of sends forth and back through TreeManagementMessages. This
     * routine is the one that sends out the requests.
     *
     * <h2> Multithreading </h2>
     *
     * I may not globally lock this routine, as I otherwise would block the
     * createObserverCloneIfRequired(). So I literally "only" protect the
     * actual push back in the vector.
     *
     * <h2> Call sequence </h2>
     *
     * The operation is used by createNewTrees() .
     */
    void addSpacetree( int masterId, int newTreeId );

    /**
     * Whenever we join two partitions, we have to stream data from the worker
     * to the master. In principle, I could omit this, as I realise a ``only
     * degenerated trees may join'' policy. However, it can happen that a tree
     * joins while other trees split or have split or join as well. While we
     * run through the joining phase, these other ranks have sent their updated
     * info to the worker that is about to join. The master that will hold the
     * data in the future is not aware (yet) of any topology changes. Therefore,
     * it is important that the worker still merges all incoming information
     * into its local tree and then forwards it (streams it) to its master. The
     * master in turn has to traverse later (in a secondary tree traversal) and
     * take this updated topological information into account. This is a stream
     * operation. Therefore, we not only have to exchange stacks, we also have
     * to revert them to transform the stack into a stream.
     *
     * All of this painful stuff is done for the vertex data only. For the user
     * data, we can rely on our traditional vertical data exchange mechanisms.
     *
     * @see peano4::stacks::STDVectorStack::revert()
     */
    void streamLocalVertexInformationToMasterThroughVerticalStacks(
      int spacetreeId, int parentId,
      const std::set<int>& joiningIds
    );

    /**
     * <h2> Merge process </h2>
     *
     * We make the worker first of all decide whether to merge or not. A
     * worker has the control/knowledge whether to join or not. When a rank
     * wants to join, we first of all run the workers. Then, all data is the
     * join buffers. After that, we sweep through the masters to take up the
     * data. So we may not merge both the worker and the master into their
     * workers at the same time. If we did so, we'd not be able to run all
     * the merging trees in parallel.
     */
    void cleanUpTrees(peano4::grid::TraversalObserver&  observer);

    /**
     * I need this routine for technical reasons: Prior to the sweep of trees,
     * I have to identify all of those trees which wanna merge with their
     * workers. This is an analysis I have to do before I actually traverse
     * any worker. Because of this traversal, more trees might denote their
     * workers as joining, so if I query the joining field after the traversal,
     * I intermix newly joining and old joining ranks.
     */
    std::set<int>                 getLocalTreesMergingWithWorkers() const;

    /**
     * Quick lookup whether an observer clone for this tree id does already exist.
     * If not, we create one quickly.
     *
     * \section Multithreading
     *
     * This operation uses a lock on the semaphore to ensure that no two threads
     * insert an observer into the global table at the same time.
     */
    void createObserverCloneIfRequired(peano4::grid::TraversalObserver& observer, int treeId);

    SpacetreeSet();
    SpacetreeSet(const SpacetreeSet& ) = delete;
    SpacetreeSet& operator=(const SpacetreeSet& ) = delete;

  public:
    /**
     * As the setis a singleton and a service, it has to deregister itself. Otherwise, the overall
     * service landscape still might try to call receiveDanglingMessages().
     */
    ~SpacetreeSet();

    static SpacetreeSet& getInstance();


    /**
     * Run through the set of unanswered questions and, well, answer them.
     *
     * The following messages are received/handled:
     *
     * - Action::RequestNewRemoteTree (invoked by split())
     * - Action::CreateNewRemoteTree (invoked by addSpacetree(int,int))
     * - Action::RemoveChildTreeFromBooksAsChildBecameEmpty (invoked by cleanUpTrees())
     *
     * split() is something any rank can trigger at any time. Most default load
     * balancers call it throughout the end of the grid sweep. addSpacetree() is
     * called by the set through createNewTrees() just before all observers are
     * deleted and the traversal terminates. cleanUpTrees() just arises before
     * that one.
     *
     * No all of these messages can be answered at any point of the local grid
     * sweeps. That is, we may never add a tree to the local tree collection while
     * we are right in the middle of traversals on this rank, e.g.
     * We may not insert or
     * remove trees from the global data structure while we are still traversing
     * some local trees or run some data exchange. I originally tried to protect
     * the whole code by a _state==SpacetreeSetState::Waiting check, i.e. to make
     * the set react if and only if we are not doing anything anyway. That did
     * not work properly, as tree requests for example are to be handled
     * immediately. So what I do now is that I pipe all requests into a vector.
     * Then, I run through the vector and, depending on the set state, do answer
     * messages or leave them in the queue for the time being.
     *
     * This approach leads to errors whenever a message send-out is followed by
     * a second message that provides further details. The first message might be
     * buffered locally, and, as we can't answer the first one immediately, the
     * second message (of another datatype) will be interpreted as yet another
     * request. So that means that every single message exchange with the set
     * always has to be follow a send-acknowledge pattern.
     *
     * <h2> Action::RequestNewRemoteTree </h2>
     *
     * Can be answered any time, as it literally just books a new number but
     * nothing more happens at this point. The acknowledgement message carries
     * the new tree's number.
     *
     * <h2> Action::CreateNewRemoteTree </h2>
     *
     * This message triggers the insertation of a new tree into the local set of
     * spacetrees. We thus may handle it if and only if we are the end of a
     * traversal. The message logically consists of two parts: The master of a
     * new tree (which has previously used Action::RequestNewRemoteTree to get a
     * new tree's number) sends out the create message, waits for the
     * acknowledgement and then sends out the tree state, to the rank hosting the
     * new tree can actually create the data structure. This last step is followed
     * by an acknowledge message which carries no particular information. That
     * is, this message exchange belongs to the one-way information flow, but we
     * have one acknowledgement message to say "go ahead", and one message to
     * say "ok, we got it".
     *
     * <h2> Action::RemoveChildTreeFromBooksAsChildBecameEmpty </h2>
     *
     * This is a simple one though we have again to ensure that it is handled if
     * and only if we have finishes the local traversals and thus can safely
     * manipulate the local spacetree.
     *
     * <h2> Data consistency </h2>
     *
     * In this routine, we have to be very careful with the data consistency. While
     * we are running through the set of unanswered messages, new ones might drop
     * in. However, C++ is not very happy if a vector is added further elements
     * while we are traversing it (in theory, it might happen that it is expanded
     * and thus copied). So what I do is something different: I first run through
     * the vector and take those out that I know that I can handle. Then I handle
     * these guys and return - well-aware that meanwhile further messages might
     * have dropped in. Which is something I don't care, as I rely on the calling
     * code just to invoke answer again if it is important for the code's progress.
     *
     *
     * @todo replyToUnansweredMessages() sollte der Name sein
     */
    void answerQuestions();

    /**
     * We poll the tree management messages.
     *
     * Messages are not answered directly. Instead, we buffer them in
     * _unansweredMessages and hand over to replyToUnansweredMessages().
     *
     * See the description in mpi.h for a howto where this routine is
     * used and how to use it.
     */
    virtual void receiveDanglingMessages() override;

    /**
     * @see Spacetree::Spacetree()
     * @see shutdown()
     */
    void init(
      const tarch::la::Vector<Dimensions,double>&  offset,
      const tarch::la::Vector<Dimensions,double>&  width,
      const std::bitset<Dimensions>&               periodicBC = 0
    );

    virtual void shutdown() override;

    /**
     * Invoke traverse on all spacetrees in parallel.
     *
     * <h2> Sequence </h2>
     *
     * It is important that I don't answer to spacetree request messages while I am
     * running through the spacetree set. Some of these request messages add further
     * elements to the set. However, the set should remain invariant while I run
     * through it. So I introduced the _state. See also _unansweredMessages.
     *
     * So when I'm done with the local traversal, I have to ensure that all other sets
     * have got their tree modifications through. That is: If another rank wants to
     * insert a tree, that has to happen in the right traversal. As the other rank
     * will wait for an acknowledgement of the addition, there is no risk that this
     * rank has already proceeded wrongly into the next grid sweep. However, the target
     * rank has to answer the request in the right one and may not proceed to early.
     *
     * So I introduce a barrier. The barrier has to do two things: On the one hand, it
     * has to receive dangling messages. On the other hand, it should answer messages
     * that it has not answered before. In principle, this is indirectly done through
     * receiveDanglingMessages(). However, I played around with running the dangling
     * thing if and only if iprobe tells me that there are new messages. So to be on
     * the safe side, I rather invoke the answer routines manually here.
     *
     * There are two options where to place the barrier: We could add it to the end of
     * traverse() or right after we've received the startup message. There are two different
     * kinds of messages: messages that can be answered straightaway (give me  a new
     * rank) or others which can be answered only in-between iterations. If we add
     * a barrier right at the begin of a traversal, then we are fine, as the typical
     * split (asking for further ranks) arises in-between traversals.
     */
    void traverse(peano4::grid::TraversalObserver& observer);

    /**
     * Return statistics object for primary spacetree.
     */
    peano4::grid::GridStatistics  getGridStatistics() const;

    peano4::grid::GridStatistics  getGridStatistics(int treeId) const;

    /**
     * Split a local tree.
     *
     * If the target tree shall be stored on the local node, then you pass
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * peano4::parallel::Node::getInstance().getRank(treeId)
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * as last argument. Splits on the local node allow Peano to exploit more
     * cores on the node. The details of the total splitting procedure are
     * discussed in @ref peano_domain_decomposition "Peano's generic domain decomposition" discussion.
     *
     *
     * @param treeId  Id of tree that should split, i.e. give away cells to a
     *   new tree
     * @param cells   Number of cells to deploy to a new tree. This count refers
     *   to the status-quo, i.e. if dynamic adaptivity makes an unrefined node
     *   of the tree unfold in the next step, the new @f$ 3^d @f$ children do
     *   count only with one towards this quota.
     * @param targetRank  Rank which should host the new tree.
     */
    bool split(int treeId, const peano4::SplitInstruction& instruction, int targetRank);

    /**
     * Codes hold one spacetree set per rank. With this routine, you can find
     * out whether a local set contains a particular id.
     */
    bool isLocalSpacetree(int treeId) const;

    /**
     * @return Set of ids of local spacetrees
     */
    std::set<int> getLocalSpacetrees() const;
};

#include "SpacetreeSet.cpph"
