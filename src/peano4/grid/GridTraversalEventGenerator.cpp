#include "GridTraversalEventGenerator.h"
#include "grid.h"
#include "PeanoCurve.h"
#include "Spacetree.h"
#include "TraversalObserver.h"

#include "peano4/utils/Loop.h"

tarch::logging::Log  peano4::grid::GridTraversalEventGenerator::_log( "peano4::grid::GridTraversalEventGenerator" );

peano4::grid::GridTraversalEventGenerator::GridTraversalEventGenerator(int id):
  _id(id) {}

std::bitset<TwoPowerD> peano4::grid::GridTraversalEventGenerator::areVerticesLocal(
  [[maybe_unused]] GridVertex  vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining
) const {
  std::bitset<TwoPowerD> bitset;
  for (int i=0; i<TwoPowerD; i++) {
    bitset.set(i,isVertexAdjacentToLocalSpacetree(vertices[i], splitTriggered, splitting, joinTriggered, joining, true, true));
  }
  logDebug( "areVerticesLocal(...)", bitset );
  return bitset;
}

std::bitset<TwoTimesD> peano4::grid::GridTraversalEventGenerator::areFacesLocal(
  [[maybe_unused]] GridVertex                  vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&   splitTriggered,
  [[maybe_unused]] const std::set<int>&        splitting,
  [[maybe_unused]] const std::set< int >&      joinTriggered,
  [[maybe_unused]] const std::set< int >&      joining
) const {
  std::bitset<TwoTimesD> result;
  for (int faceNumber=0; faceNumber<2*Dimensions; faceNumber++) {
    bool isLocal = false;

    const int   normal = faceNumber % Dimensions;
    for (int i=0; i<TwoPowerD; i++) {
      std::bitset<Dimensions> studiedVertex = i;
      studiedVertex.set(normal,faceNumber>=Dimensions);
      std::bitset<Dimensions> studiedEntry  = TwoPowerD - studiedVertex.to_ulong() - 1;

      studiedEntry.set(normal,0);
      int currentRank = vertices[studiedVertex.to_ulong()].getAdjacentRanks( studiedEntry.to_ulong() );
      if (
        vertices[studiedVertex.to_ulong()].getState()!=GridVertex::State::HangingVertex
        or
        faceNumber>=Dimensions
      ) {
        isLocal |=  currentRank == _id;
        isLocal |=  splitTriggered.count(currentRank)>0;
        isLocal |=  splitting.count(currentRank)>0;
      }

      studiedEntry.set(normal,1);
      if (
        vertices[studiedVertex.to_ulong()].getState()!=GridVertex::State::HangingVertex
        or
        faceNumber<Dimensions
      ) {
        currentRank = vertices[studiedVertex.to_ulong()].getAdjacentRanks( studiedEntry.to_ulong() );
        isLocal |= currentRank == _id;
        isLocal |=  splitTriggered.count(currentRank)>0;
        isLocal |=  splitting.count(currentRank)>0;
      }
    }

    result[faceNumber] = isLocal;
  }
  return result;
}

bool peano4::grid::GridTraversalEventGenerator::isVertexAdjacentToLocalSpacetree(
  [[maybe_unused]] GridVertex                 vertex,
  [[maybe_unused]] const SplitSpecification&  splitTriggered,
  [[maybe_unused]] const std::set<int>&       splitting,
  [[maybe_unused]] const std::set< int >&     joinTriggered,
  [[maybe_unused]] const std::set< int >&     joining,
  [[maybe_unused]] bool        splittingIsConsideredLocal,
  [[maybe_unused]] bool        joiningIsConsideredLocal
) const {
/*
  if (vertex.getState()==GridVertex::State::HangingVertex) {
    return false;
  }
  else {
*/
    logTraceInWith3Arguments( "isVertexAdjacentToLocalSpacetree(...)", vertex.toString(), splittingIsConsideredLocal, joiningIsConsideredLocal );
    bool result = false;
    for(int i=0; i<TwoPowerD; i++) {
      assertion( splitTriggered.count( vertex.getAdjacentRanks(i) )<=1 );
      assertion( splitting.count( vertex.getAdjacentRanks(i) )<=1 );

      result |= vertex.getAdjacentRanks(i)==_id;

      result |= splitTriggered.count( vertex.getAdjacentRanks(i) )==1;
      result |= (splittingIsConsideredLocal and splitting.count( vertex.getAdjacentRanks(i) )==1);

      result |= (joiningIsConsideredLocal and joining.count( vertex.getAdjacentRanks(i) )==1);
    }
    logTraceOutWith1Argument( "isVertexAdjacentToLocalSpacetree(...)", result );
    return result;
//  }
}

bool peano4::grid::GridTraversalEventGenerator::isSpacetreeNodeLocal(
  [[maybe_unused]] GridVertex    vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] bool          splittingIsConsideredLocal,
  [[maybe_unused]] bool          joiningIsConsideredLocal
) const {
  bool isLocal = true;
  dfor2(k)
    isLocal &= (
      (vertices[kScalar].getState()==GridVertex::State::HangingVertex)
      or
      (
        vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1)==_id
      )
      or
      ( splitTriggered.count(vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1)) > 0)
      or
      (
        splittingIsConsideredLocal and splitting.count(vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1))>0
      )
      or
      (
        joiningIsConsideredLocal and joining.count(vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1))>0
      )
    );
  enddforx

  return isLocal;
}

std::bitset<TwoPowerD> peano4::grid::GridTraversalEventGenerator::areVerticesAdjacentToParallelDomainBoundary(
  [[maybe_unused]] GridVertex                                vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] bool calledByLeaveCell
) const {
  std::bitset<TwoPowerD> result;
  for (int j=0; j<TwoPowerD; j++) {
    tarch::la::Vector< TwoPowerD, int > adjacency = vertices[j].getBackupOfAdjacentRanks();
    bool oneLocal  = false;
    bool oneRemote = false;
    for (int i=0; i<TwoPowerD; i++ ) {
      bool valid = adjacency(i)>=0 or adjacency(i)==peano4::grid::Spacetree::RankOfPeriodicBoundaryCondition;
      bool local = adjacency(i)==_id or splitting.count(adjacency(i))>0 or splitTriggered.count(adjacency(i))>0;
      oneLocal  |= (valid and local);
      oneRemote |= (valid and not local);
    }
    result.set(j,oneLocal and oneRemote);
  }
  return result;
}

std::bitset<TwoTimesD> peano4::grid::GridTraversalEventGenerator::areFacesAdjacentToParallelDomainBoundary(
  [[maybe_unused]] GridVertex                                vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] bool calledByLeaveCell
) const {
  std::bitset<TwoTimesD> result;
  for (int j=0; j<TwoTimesD; j++) {
    // @todo
    tarch::la::Vector< TwoPowerD, int > adjacency = getAdjacentRanksOfFace(vertices, j, calledByLeaveCell);
    //tarch::la::Vector< TwoPowerD, int > adjacency = getAdjacentRanksOfFace(vertices, i, false);
    bool oneLocal  = false;
    bool oneRemote = false;
    for (int i=0; i<TwoPowerD; i++ ) {
      // @todo noch rueber auf vertices
      // So ganz stimmt es natuerlich net, weil das ist ja jetzt eins zu frueh
      bool valid = adjacency(i)>=0 or adjacency(i)==peano4::grid::Spacetree::RankOfPeriodicBoundaryCondition;
/*
      bool local = adjacency(i)==_id
                or (splitting.count(adjacency(i))>0 and calledByEnterCell)
                or splitTriggered.count(adjacency(i))>0;
      bool remote = adjacency(i)!=_id
                and not splitTriggered.count(adjacency(i))>0
                and (not splitting.count(adjacency(i))>0 or not calledByEnterCell);
*/
      bool local = adjacency(i)==_id
                or splitTriggered.count(adjacency(i))>0
//                or splitting.count(adjacency(i))>0;
                or (splitting.count(adjacency(i))>0 and not calledByLeaveCell);
      bool remote = adjacency(i)!=_id
                and not (splitTriggered.count(adjacency(i))>0)
                and (splitting.count(adjacency(i))==0 or calledByLeaveCell);
//      and not (splitTriggered.count(adjacency(i))>0 and calledByLeaveCell);
//      bool remote = not local;
      oneLocal  |= (valid and local);
      oneRemote |= (valid and remote);
    }
    result.set(j,oneLocal and oneRemote);
  }
  return result;
}

peano4::grid::GridTraversalEvent peano4::grid::GridTraversalEventGenerator::createGenericCellTraversalEvent(
  [[maybe_unused]] GridVertex                                coarseGridVertices[TwoPowerD],
  [[maybe_unused]] GridVertex                                fineGridVertices[TwoPowerD],
  [[maybe_unused]] const AutomatonState&                     state,
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] const tarch::la::Vector<Dimensions,int>&  relativePositionToFather,
  [[maybe_unused]] bool                                      spacetreeStateIsRunning
) const {
  #if Dimensions==2
  logTraceInWith7Arguments( "createGenericCellTraversalEvent(...)", state.toString(), relativePositionToFather, fineGridVertices[0].toString(), fineGridVertices[1].toString(), fineGridVertices[2].toString(), fineGridVertices[3].toString(), _id );
  #else
  logTraceInWith3Arguments( "createGenericCellTraversalEvent(...)", state.toString(), relativePositionToFather, _id );
  #endif
  GridTraversalEvent  event;
  event.setX( state.getX() + state.getH()*0.5 );
  event.setH( state.getH() );

  event.setHasBeenRefined( haveVerticesBeenRefined(fineGridVertices) );
  event.setWillBeRefined( willVerticesBeRefined(fineGridVertices) );
  event.setRelativePositionToFather( relativePositionToFather );

  event.setIsCellLocal( isSpacetreeNodeLocal( fineGridVertices,   splitTriggered, splitting, joinTriggered, joining, true, true) );
  event.setIsParentCellLocal( state.getLevel()>1 and isSpacetreeNodeLocal( coarseGridVertices, splitTriggered, splitting, joinTriggered, joining, true, true) );
  event.setIsFaceLocal( areFacesLocal( fineGridVertices, splitTriggered, splitting, joinTriggered, joining) );
  event.setIsVertexLocal( areVerticesLocal( fineGridVertices, splitTriggered, splitting, joinTriggered, joining) );

  event.setIsParentVertexLocal( areVerticesLocal( coarseGridVertices, splitTriggered, splitting, joinTriggered, joining) );

  event.setInvokingSpacetree( _id );
  event.setInvokingSpacetreeIsNotInvolvedInAnyDynamicLoadBalancing(
        spacetreeStateIsRunning and
        joinTriggered.empty() and
        joining.empty() and
        splitTriggered.empty() and
        splitting.empty()
  );

  event.setIsAdjacentCellLocal(0);
  dfor2(k) {
    if (fineGridVertices[kScalar].getState()!=GridVertex::State::HangingVertex) {
      dfor2(j) {
        tarch::la::Vector<Dimensions,int>  entryIn3x3Patch = k + j;
        int owningTree = fineGridVertices[kScalar].getAdjacentRanks(jScalar);
        event.setIsAdjacentCellLocal(
          peano4::utils::dLinearised(entryIn3x3Patch,3),
          owningTree == _id
        );
        event.setIsVertexParentOfSubtree( fineGridVertices[kScalar].getHasBeenParentOfSubtreeVertexInPreviousTreeSweep() );
      } enddforx
    }
  } enddforx

  logTraceOutWith2Arguments( "createGenericCellTraversalEvent(...)", event.toString(), _id );
  return event;
}

peano4::grid::GridTraversalEvent peano4::grid::GridTraversalEventGenerator::createLeaveCellTraversalEvent(
  [[maybe_unused]] GridVertex                                coarseGridVertices[TwoPowerD],
  [[maybe_unused]] GridVertex                                fineGridVertices[TwoPowerD],
  [[maybe_unused]] const AutomatonState&                     state,
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] const std::set< int >&                    hasSplit,
  [[maybe_unused]] const tarch::la::Vector<Dimensions,int>&  relativePositionToFather,
  [[maybe_unused]] bool                                      spacetreeStateIsRunning
) const {
  logTraceInWith3Arguments( "createLeaveCellTraversalEvent(...)", state.toString(), _id, relativePositionToFather );
  GridTraversalEvent  event = createGenericCellTraversalEvent(coarseGridVertices, fineGridVertices, state, splitTriggered, splitting, joinTriggered, joining, relativePositionToFather, spacetreeStateIsRunning);

  event.setIsVertexAdjacentToParallelDomainBoundary( areVerticesAdjacentToParallelDomainBoundary(fineGridVertices, splitTriggered, splitting, joinTriggered, joining, true) );
  event.setIsFaceAdjacentToParallelDomainBoundary( areFacesAdjacentToParallelDomainBoundary(fineGridVertices, splitTriggered, splitting, joinTriggered, joining, true));

  const std::bitset<Dimensions> coordinates = PeanoCurve::getFirstVertexIndex(state);
  for (int i=0; i<TwoPowerD; i++) {
    const std::bitset<Dimensions>           vertexIndex( coordinates ^ std::bitset<Dimensions>(i) );
    const int   stackNumber = PeanoCurve::getVertexWriteStackNumber(state,vertexIndex);

    event.setVertexDataFrom(i,vertexIndex.to_ulong());
    switch ( fineGridVertices[vertexIndex.to_ulong()].getState() ) {
      case GridVertex::State::HangingVertex:
        event.setVertexDataTo(i,TraversalObserver::CreateOrDestroyHangingGridEntity);
        break;
      case GridVertex::State::New:
      case GridVertex::State::Unrefined:
      case GridVertex::State::Refined:
      case GridVertex::State::RefinementTriggered:
      case GridVertex::State::Refining:
      case GridVertex::State::EraseTriggered:
      case GridVertex::State::Erasing:
        event.setVertexDataTo(i,stackNumber);
        break;
      case GridVertex::State::Delete:
        {
          if ( PeanoCurve::isInOutStack(stackNumber) ) {
            event.setVertexDataTo(i,TraversalObserver::CreateOrDestroyPersistentGridEntity);
          }
          else {
            event.setVertexDataTo(i,stackNumber);
          }
        }
        break;
    }

    if (
      PeanoCurve::isInOutStack(event.getVertexDataTo(i))
      and
      not event.getIsVertexLocal(vertexIndex.to_ulong())
    ) {
      logDebug(
        "createLeaveCellTraversalEvent(...)",
        "reset event entry " << i << " to no-data as we have " << event.toString()
        << ". vertex=" << fineGridVertices[vertexIndex.to_ulong()].toString()
      );
      event.setVertexDataTo(i,TraversalObserver::NoData);
    }
  }

  for (int i=0; i<2*Dimensions; i++) {
    int        faceIndex   = PeanoCurve::getFaceNumberAlongCurve(state,i);
    FaceType   type        = getFaceType(coarseGridVertices,relativePositionToFather,faceIndex);
    const int  stackNumber = PeanoCurve::getFaceWriteStackNumber(state,faceIndex);

    event.setFaceDataFrom(i,faceIndex);

    switch (type) {
      case FaceType::Hanging:
        event.setFaceDataTo(i,TraversalObserver::CreateOrDestroyHangingGridEntity);
        break;
      case FaceType::New:
      case FaceType::Persistent:
        event.setFaceDataTo(i,stackNumber);
        break;
      case FaceType::Delete:
        if ( PeanoCurve::isInOutStack(stackNumber) ) {
          event.setFaceDataTo(i,TraversalObserver::CreateOrDestroyPersistentGridEntity);
        }
        else {
          event.setFaceDataTo(i,stackNumber);
        }
        break;
    }

    if (
      PeanoCurve::isInOutStack(event.getFaceDataTo(i))
      and
      not event.getIsFaceLocal(faceIndex)
    ) {
      event.setFaceDataTo(i,TraversalObserver::NoData);
    }
  }

  {
    CellType type = getCellType(coarseGridVertices,relativePositionToFather);
    const int  stackNumber = PeanoCurve::getCellWriteStackNumber(state);

    switch (type) {
      case CellType::New:
        event.setCellData(stackNumber);
        break;
          case CellType::Persistent:
        event.setCellData(stackNumber);
        break;
      case CellType::Delete:
        event.setCellData(TraversalObserver::CreateOrDestroyPersistentGridEntity);
        break;
    }

    if (
      // always true here, but if I write it here explicitly, then it is consistent with faces/vertices
      PeanoCurve::isInOutStack(event.getCellData())
      and
      not event.getIsCellLocal()
    ) {
      event.setCellData(TraversalObserver::NoData);
    }
  }

  #if !defined(SharedMemoryParallelisation)
  for (int i=0; i<TwoPowerD; i++)
    assertion1( event.getIsVertexLocal(i) or not event.getIsVertexAdjacentToParallelDomainBoundary(i), event.toString() );
  for (int i=0; i<TwoTimesD; i++)
    assertion1( event.getIsFaceLocal(i) or not event.getIsFaceAdjacentToParallelDomainBoundary(i), event.toString() );
  #endif

  logTraceOutWith3Arguments( "createLeaveCellTraversalEvent(...)", state.toString(), event.toString(), _id );
  return event;
}

peano4::grid::GridTraversalEvent peano4::grid::GridTraversalEventGenerator::createPrunedEnterCellTraversalEvent( SpacetreeState spacetreeState, const GridTraversalEvent& event ) const {
  GridTraversalEvent result = event;

  if(
    spacetreeState==SpacetreeState::EmptyRun or
    spacetreeState==SpacetreeState::NewFromSplit or
    spacetreeState==SpacetreeState::Joining
  ) {
    result.setIsCellLocal(false);
    result.setIsParentCellLocal(false);
    result.setIsFaceLocal(0);
    result.setIsVertexLocal(0);
  }

/*
  if(
    spacetreeState==SpacetreeState::EmptyRun or
    spacetreeState==SpacetreeState::NewFromSplit
  ) {
    result.setVertexDataFrom(TraversalObserver::NoData);
    result.setFaceDataFrom(TraversalObserver::NoData);
    result.setCellData(TraversalObserver::NoData);
  }
*/

  return result;
}


peano4::grid::GridTraversalEvent peano4::grid::GridTraversalEventGenerator::createPrunedLeaveCellTraversalEvent( SpacetreeState spacetreeState, const GridTraversalEvent& event ) const {
  GridTraversalEvent result = event;

  if(
    spacetreeState==SpacetreeState::EmptyRun or
    spacetreeState==SpacetreeState::NewFromSplit or
    spacetreeState==SpacetreeState::Joining
  ) {
    result.setIsCellLocal(false);
    result.setIsParentCellLocal(false);
    result.setIsFaceLocal(0);
    result.setIsVertexLocal(0);
  }

  /*
  if(
    spacetreeState==SpacetreeState::EmptyRun or
//    spacetreeState==SpacetreeState::NewFromSplit or
    spacetreeState==SpacetreeState::Joining
  ) {
    result.setVertexDataTo(TraversalObserver::NoData);
    result.setFaceDataTo(TraversalObserver::NoData);
    result.setCellData(TraversalObserver::NoData);
  }
*/

  return result;
}

peano4::grid::GridTraversalEvent peano4::grid::GridTraversalEventGenerator::createEnterCellTraversalEvent(
  [[maybe_unused]] GridVertex                                coarseGridVertices[TwoPowerD],
  [[maybe_unused]] GridVertex                                fineGridVertices[TwoPowerD],
  [[maybe_unused]] const AutomatonState&                     state,
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining,
  [[maybe_unused]] const std::set< int >&                    hasSplit,
  [[maybe_unused]] const tarch::la::Vector<Dimensions,int>&  relativePositionToFather,
  [[maybe_unused]] bool                                      spacetreeStateIsRunning
) const {
  logTraceInWith7Arguments( "createEnterCellTraversalEvent(...)", state.toString(), _id, relativePositionToFather, coarseGridVertices[0].toString(), coarseGridVertices[1].toString(), coarseGridVertices[2].toString(), coarseGridVertices[3].toString() );
  GridTraversalEvent  event = createGenericCellTraversalEvent(coarseGridVertices, fineGridVertices, state, splitTriggered, splitting, joinTriggered, joining, relativePositionToFather, spacetreeStateIsRunning);

  event.setIsVertexAdjacentToParallelDomainBoundary( areVerticesAdjacentToParallelDomainBoundary(fineGridVertices, splitTriggered, splitting, joinTriggered, joining, false) );
  event.setIsFaceAdjacentToParallelDomainBoundary( areFacesAdjacentToParallelDomainBoundary(fineGridVertices, splitTriggered, splitting, joinTriggered, joining, false));

  #if !defined(SharedMemoryParallelisation)
  for (int i=0; i<TwoPowerD; i++)
    assertion1( event.getIsVertexLocal(i) or not event.getIsVertexAdjacentToParallelDomainBoundary(i), event.toString() );
  for (int i=0; i<TwoTimesD; i++)
    assertion1( event.getIsFaceLocal(i) or not event.getIsFaceAdjacentToParallelDomainBoundary(i), event.toString() );
  #endif

  const std::bitset<Dimensions> coordinates = PeanoCurve::getFirstVertexIndex(state);
  for (int i=0; i<TwoPowerD; i++) {
    const std::bitset<Dimensions>  vertexIndex( coordinates ^ std::bitset<Dimensions>(i) );
    const int  stackNumber    = PeanoCurve::getVertexReadStackNumber(state,vertexIndex);
    const int  vertexPosition = vertexIndex.to_ulong();

    switch ( fineGridVertices[vertexPosition].getState() ) {
      case GridVertex::State::HangingVertex:
        event.setVertexDataFrom(i,TraversalObserver::CreateOrDestroyHangingGridEntity);
        break;
      case GridVertex::State::New:
        {
          if ( PeanoCurve::isInOutStack(stackNumber) ) {
            event.setVertexDataFrom(i,TraversalObserver::CreateOrDestroyPersistentGridEntity);
          }
          else {
            event.setVertexDataFrom(i,stackNumber);
          }
        }
        break;
      case GridVertex::State::Unrefined:
      case GridVertex::State::Refined:
      case GridVertex::State::RefinementTriggered:
      case GridVertex::State::Refining:
      case GridVertex::State::EraseTriggered:
      case GridVertex::State::Erasing:
      case GridVertex::State::Delete:
        event.setVertexDataFrom(i,stackNumber);
        break;
    }
    event.setVertexDataTo(i,vertexIndex.to_ulong());

    bool mayResetToNoData =
      PeanoCurve::isInOutStack(event.getVertexDataFrom(i))
      and
      not event.getIsVertexLocal(vertexPosition);

    for (auto p: hasSplit) {
      mayResetToNoData &= not tarch::la::contains( fineGridVertices[vertexPosition].getAdjacentRanks(), p );
    }

    if (mayResetToNoData) {
      event.setVertexDataFrom(i,TraversalObserver::NoData);
    }

    logDebug(
      "createEnterCellTraversalEvent(...)",
      "vertex " << i << " on on tree " << _id << ": " <<
      event.toString() << ", mayResetToNoData=" << mayResetToNoData <<
      ", is-local=" << event.getIsVertexLocal(vertexPosition) << ", is in-out=" << PeanoCurve::isInOutStack(event.getVertexDataFrom(i))
      << ", vertex[0]=" << fineGridVertices[0].toString()
      << ", vertex[1]=" << fineGridVertices[1].toString()
      << ", vertex[2]=" << fineGridVertices[2].toString()
      << ", vertex[3]=" << fineGridVertices[3].toString()
      << ", vertex[4]=" << fineGridVertices[4].toString()
      << ", vertex[5]=" << fineGridVertices[5].toString()
      << ", vertex[6]=" << fineGridVertices[6].toString()
      << ", vertex[7]=" << fineGridVertices[7].toString()
    );
  }

  for (int i=0; i<2*Dimensions; i++) {
    int        faceIndex   = PeanoCurve::getFaceNumberAlongCurve(state,i);
    FaceType   type        = getFaceType(coarseGridVertices,relativePositionToFather,faceIndex);
    const int  stackNumber = PeanoCurve::getFaceReadStackNumber(state,faceIndex);

    switch (type) {
      case FaceType::New:
        {
          if ( PeanoCurve::isInOutStack(stackNumber) ) {
            event.setFaceDataFrom(i,TraversalObserver::CreateOrDestroyPersistentGridEntity);
          }
          else {
            event.setFaceDataFrom(i,stackNumber);
          }
        }
        break;
      case FaceType::Hanging:
        event.setFaceDataFrom(i,TraversalObserver::CreateOrDestroyHangingGridEntity);
        break;
      case FaceType::Persistent:
      case FaceType::Delete:
        event.setFaceDataFrom(i,stackNumber);
        break;
    }
    event.setFaceDataTo(i,faceIndex);

    bool mayResetToNoData =
      PeanoCurve::isInOutStack(event.getFaceDataFrom(i))
      and
      not event.getIsFaceLocal(faceIndex);

    for (auto p: hasSplit) {
      mayResetToNoData &= not tarch::la::contains( getAdjacentRanksOfFace(fineGridVertices, faceIndex, false), p );
    }

    if (mayResetToNoData) {
      event.setFaceDataFrom(i,TraversalObserver::NoData);
    }
  }

  {
    CellType type = getCellType(coarseGridVertices,relativePositionToFather);
    const int  stackNumber = PeanoCurve::getCellReadStackNumber(state);
    switch (type) {
      case CellType::New:
        event.setCellData(TraversalObserver::CreateOrDestroyPersistentGridEntity);
        break;
          case CellType::Persistent:
        event.setCellData(stackNumber);
        break;
      case CellType::Delete:
        event.setCellData(stackNumber);
        break;
    }
  }

  bool mayResetToNoData =
    // always true here, but if I write it here explicitly, then it is consistent with faces/vertices
    PeanoCurve::isInOutStack(event.getCellData())
    and
    not event.getIsCellLocal();

  for (auto p: hasSplit) {
    mayResetToNoData &= getTreeOwningSpacetreeNode(fineGridVertices, splitTriggered, splitting, joinTriggered, joining)!=p;
  }

  if (mayResetToNoData) {
    event.setCellData(TraversalObserver::NoData);
  }

  #if !defined(SharedMemoryParallelisation)
  for (int i=0; i<TwoPowerD; i++)
    assertion1( event.getIsVertexLocal(i) or not event.getIsVertexAdjacentToParallelDomainBoundary(i), event.toString() );
  for (int i=0; i<TwoTimesD; i++)
    assertion1( event.getIsFaceLocal(i) or not event.getIsFaceAdjacentToParallelDomainBoundary(i), event.toString() );
  #endif

  logTraceOutWith3Arguments( "createEnterCellTraversalEvent(...)", state.toString(), event.toString(), _id );
  return event;
}

peano4::grid::FaceType peano4::grid::GridTraversalEventGenerator::getFaceType(
  [[maybe_unused]] GridVertex                         coarseGridVertices[TwoPowerD],
  [[maybe_unused]] tarch::la::Vector<Dimensions,int>  positionOfCell,
  [[maybe_unused]] int                                faceNumber
) {
  logTraceInWith1Argument( "getFaceType(...)", faceNumber );

  bool allVerticesAreHanging         = true;
  bool allVerticesAreDeleteOrHanging = true;
  bool allVerticesAreNewOrHanging    = true;

  const int normal = faceNumber % Dimensions;
  for (int i=0; i<TwoPowerD; i++) {
    std::bitset<Dimensions> studiedVertex = i;
    studiedVertex.set(normal,faceNumber>=Dimensions);
    switch ( getVertexType( coarseGridVertices, positionOfCell + tarch::la::Vector<Dimensions,int>(studiedVertex) ) ) {
      case VertexType::Hanging:
        break;
      case VertexType::New:
        allVerticesAreDeleteOrHanging = false;
        allVerticesAreHanging         = false;
        break;
      case VertexType::Persistent:
        allVerticesAreHanging         = false;
        allVerticesAreDeleteOrHanging = false;
        allVerticesAreNewOrHanging    = false;
        break;
      case VertexType::Delete:
        allVerticesAreHanging         = false;
        allVerticesAreNewOrHanging    = false;
        break;
    }
  }

  FaceType  result = FaceType::Persistent;
  if ( allVerticesAreHanging ) {
    result = FaceType::Hanging;
  }
  else if ( allVerticesAreNewOrHanging ) {
    result = FaceType::New;
  }
  else if ( allVerticesAreDeleteOrHanging ) {
    result = FaceType::Delete;
  }

  logTraceOutWith4Arguments( "getFaceType(...)", toString(result), allVerticesAreDeleteOrHanging, allVerticesAreNewOrHanging, allVerticesAreHanging );
  return result;
}

peano4::grid::CellType peano4::grid::GridTraversalEventGenerator::getCellType(
  [[maybe_unused]] GridVertex                         coarseGridVertices[TwoPowerD],
  [[maybe_unused]] tarch::la::Vector<Dimensions,int>  positionOfCell
) {
  bool allVerticesAreDelete = true;
  bool allVerticesAreNew    = true;

  for (int i=0; i<TwoPowerD; i++) {
        std::bitset<Dimensions> vectorPosition = i;
    switch ( getVertexType( coarseGridVertices, positionOfCell + tarch::la::Vector<Dimensions,int>(vectorPosition) )) {
      case VertexType::Hanging:
        break;
      case VertexType::Persistent:
        allVerticesAreDelete = false;
        allVerticesAreNew    = false;
        break;
      case VertexType::New:
        allVerticesAreDelete = false;
        break;
      case VertexType::Delete:
        allVerticesAreNew    = false;
        break;
    }
  }

  assertion( not (allVerticesAreDelete and allVerticesAreNew) );

  if (allVerticesAreDelete) {
    logDebug( "getCellType(...)", "delete cell" );
    return CellType::Delete;
  }
  else if (allVerticesAreNew) {
    logDebug( "getCellType(...)", "create new cell" );
    return CellType::New;
  }
  else {
    logDebug( "getCellType(...)", "keep cell" );
    return CellType::Persistent;
  }
}

peano4::grid::VertexType peano4::grid::GridTraversalEventGenerator::getVertexType(
  [[maybe_unused]] GridVertex                         coarseGridVertices[TwoPowerD],
  [[maybe_unused]] tarch::la::Vector<Dimensions,int>  position,
  [[maybe_unused]] int                                dimension
) {
  if (dimension<0) {
    switch (coarseGridVertices[ peano4::utils::dLinearised(position,2) ].getState()) {
      case GridVertex::State::HangingVertex:
        return VertexType::Hanging;
      case GridVertex::State::Delete:
      case GridVertex::State::New:
      case GridVertex::State::Unrefined:
        return VertexType::Hanging;
      case GridVertex::State::Refined:
        return VertexType::Persistent;
      case GridVertex::State::RefinementTriggered:
        return VertexType::Hanging;
      case GridVertex::State::Refining:
        return VertexType::New;
      case GridVertex::State::EraseTriggered:
        return VertexType::Persistent;
      case GridVertex::State::Erasing:
        return VertexType::Delete;
    }
    assertion3(false, peano4::utils::dLinearised(position,2), position, coarseGridVertices[ peano4::utils::dLinearised(position,2) ].toString() );
  }
  else if ( position(dimension)==0 ) {
    position(dimension)=0;
    return getVertexType(coarseGridVertices,position,dimension-1);
  }
  else if ( position(dimension)==3 ) {
    position(dimension)=1;
    return getVertexType(coarseGridVertices,position,dimension-1);
  }

  logTraceInWith2Arguments( "getVertexType(...)", position, dimension );
  position(dimension)=0;
  peano4::grid::VertexType lhs = getVertexType(coarseGridVertices,position,dimension-1);

  position(dimension)=1;
  peano4::grid::VertexType rhs = getVertexType(coarseGridVertices,position,dimension-1);

  VertexType result = lhs;
  if (
    (lhs==VertexType::New and rhs==VertexType::Hanging) or (lhs==VertexType::Hanging and rhs==VertexType::New)
  ) {
    result = VertexType::New;
  }
  if (
    (lhs==VertexType::New and rhs==VertexType::Persistent) or (lhs==VertexType::Persistent and rhs==VertexType::New)
  ) {
    result = VertexType::Persistent;
  }
  if (
    (lhs==VertexType::New and rhs==VertexType::Delete) or (lhs==VertexType::Delete and rhs==VertexType::New)
  ) {
    result = VertexType::Persistent;
  }
  if (
    (lhs==VertexType::Hanging and rhs==VertexType::Persistent) or (lhs==VertexType::Persistent and rhs==VertexType::Hanging)
  ) {
    result = VertexType::Persistent;
  }
  if (
    (lhs==VertexType::Hanging and rhs==VertexType::Delete) or (lhs==VertexType::Delete and rhs==VertexType::Hanging)
  ) {
    result = VertexType::Delete;
  }
  if (
    (lhs==VertexType::Persistent and rhs==VertexType::Delete) or (lhs==VertexType::Delete and rhs==VertexType::Persistent)
  ) {
    result = VertexType::Persistent;
  }

  logTraceOutWith1Argument( "getVertexType(...)", toString(result) );
  return result;
}

int peano4::grid::GridTraversalEventGenerator::getTreeOwningSpacetreeNode(
  [[maybe_unused]] GridVertex            vertices[TwoPowerD],
  [[maybe_unused]] const SplitSpecification&                 splitTriggered,
  [[maybe_unused]] const std::set<int>&                      splitting,
  [[maybe_unused]] const std::set< int >&                    joinTriggered,
  [[maybe_unused]] const std::set< int >&                    joining
) const {
  const int NotSetYet = -1;
  int id     = NotSetYet;

  int weakId = NotSetYet;
  dfor2(k)
    if (
      vertices[kScalar].getState()!=GridVertex::State::HangingVertex
      and
      vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1)!=InvalidRank
    ) {
      weakId = vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1);
    }
    if (
      vertices[kScalar].getState()!=GridVertex::State::HangingVertex
      and
      isVertexAdjacentToLocalSpacetree(vertices[kScalar], splitTriggered, splitting, joinTriggered, joining, true, false)
    ) {
      assertion8(
        id==NotSetYet
        or
        vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1)==id,
        id, kScalar, vertices[kScalar].toString(),
        vertices[0].toString(), vertices[1].toString(), vertices[2].toString(), vertices[3].toString(),
        _id
      );
      id = vertices[kScalar].getAdjacentRanks(TwoPowerD-kScalar-1);
    }
  enddforx
  assertion1(id!=NotSetYet or not isSpacetreeNodeLocal(vertices, splitTriggered, splitting, joinTriggered, joining, false, false),id);

  if (id==NotSetYet) {
    id = weakId;
  }

  return id;
}

tarch::la::Vector< TwoPowerD, int >  peano4::grid::GridTraversalEventGenerator::getAdjacentRanksOfFace( GridVertex fineGridVertices[TwoPowerD], int faceNumber, bool calledByReceivingProcess ) const {
  tarch::la::Vector< TwoPowerD, int >  adjacentRanksOfFace(InvalidRank);

  int counter = 0;
  const int normal = faceNumber % Dimensions;
  dfore( i, 2, normal, faceNumber<Dimensions ? 0 : 1 ) {
    int currentVertex = peano4::utils::dLinearised(i,2);

    std::bitset<Dimensions> studiedEntry = TwoPowerD - currentVertex - 1;
    studiedEntry[normal] = 0;
    assertion3(studiedEntry.to_ullong()>=0,        studiedEntry,currentVertex,i);
    assertion3(studiedEntry.to_ullong()<TwoPowerD, studiedEntry,currentVertex,i);
    int rankEntry = calledByReceivingProcess ? fineGridVertices[currentVertex].getBackupOfAdjacentRanks(studiedEntry.to_ullong())
                                             : fineGridVertices[currentVertex].getAdjacentRanks(studiedEntry.to_ullong());

    if (
      ( fineGridVertices[currentVertex].getState()==GridVertex::State::HangingVertex )
      or
      ( calledByReceivingProcess and fineGridVertices[currentVertex].getState()==GridVertex::State::New )
      or
      ( not calledByReceivingProcess and fineGridVertices[currentVertex].getState()==GridVertex::State::Delete )
    ) {
      rankEntry = InvalidRank;
    };

    adjacentRanksOfFace(counter) =  rankEntry;
    counter++;

    studiedEntry[normal] = 1;
    assertion3(studiedEntry.to_ullong()>=0,        studiedEntry,currentVertex,i);
    assertion3(studiedEntry.to_ullong()<TwoPowerD, studiedEntry,currentVertex,i);
    rankEntry = calledByReceivingProcess ? fineGridVertices[currentVertex].getBackupOfAdjacentRanks(studiedEntry.to_ullong())
                                         : fineGridVertices[currentVertex].getAdjacentRanks(studiedEntry.to_ullong());

    if (
      ( fineGridVertices[currentVertex].getState()==GridVertex::State::HangingVertex )
      or
      ( calledByReceivingProcess and fineGridVertices[currentVertex].getState()==GridVertex::State::New )
      or
      ( not calledByReceivingProcess and fineGridVertices[currentVertex].getState()==GridVertex::State::Delete )
    ) {
      rankEntry = InvalidRank;
    };

    adjacentRanksOfFace(counter) =  rankEntry;
    counter++;
  }

  return adjacentRanksOfFace;
}
