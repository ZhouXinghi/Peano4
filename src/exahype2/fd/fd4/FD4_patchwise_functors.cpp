#include "FD4.h"

#include "exahype2/fd/LoopBody.h"
#include "exahype2/fd/fd4/LoopBody.h"

#include "tarch/multicore/otter.h"
#include "tarch/accelerator/accelerator.h"

void exahype2::fd::fd4::timeStep_patchwise_heap_functors(
  ::exahype2::CellData&   patchData,
  int                     numberOfGridCellsPerPatchPerAxis,
  int                     haloSize,
  int                     unknowns,
  int                     auxiliaryVariables,
  double                  KOSigma,
  bool                    evaluateFlux,
  bool                    evaluateDifferentialSource, //for ncp
  bool                    evaluateAlgebraicSource, //real source
  bool                    copyOldTimeStepAndScaleWithTimeStepSize,
  DifferentialSourceTermVariant                  variant,
  Flux                    flux,
  NonconservativeProduct  DifferentialSource,
  Source                  AlgebraicSource
) {
  static tarch::logging::Log _log( "exahype2::fd" );

  assertionEquals(haloSize,3);

  bool evaluateKODissipation = KOSigma>0.0;
  bool second_order_formulation = false;

  exahype2::enumerator::AoSLexicographicEnumerator QInEnumerator (1,numberOfGridCellsPerPatchPerAxis,haloSize,unknowns,auxiliaryVariables);
  exahype2::enumerator::AoSLexicographicEnumerator QOutEnumerator(1,numberOfGridCellsPerPatchPerAxis,0,       unknowns,0);

  logDebug( "timeStep_patchwise_heap_functors(...)", "size of temp flux field=" << QOutEnumerator.size() );
  double* tempFluxX =  evaluateFlux                                                ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempFluxY =  evaluateFlux                                                ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempFluxZ = (evaluateFlux and Dimensions==3)                             ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempDifferentialSourceX =  evaluateDifferentialSource                    ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempDifferentialSourceY =  evaluateDifferentialSource                    ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempDifferentialSourceZ = (evaluateDifferentialSource and Dimensions==3) ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempKODissipationX      =  evaluateKODissipation                         ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempKODissipationY      =  evaluateKODissipation                         ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;
  double* tempKODissipationZ      = (evaluateKODissipation      and Dimensions==3) ? tarch::allocateMemory<double>( QOutEnumerator.size(),tarch::MemoryLocation::Heap ) : nullptr;

  for (int patchIndex=0; patchIndex<patchData.numberOfCells; patchIndex++) {
    const double timeStepSize = copyOldTimeStepAndScaleWithTimeStepSize ? patchData.dt[patchIndex] : 1.0;

    // STEP 1: copy solution from old array to new array, or set all to zero
    if ( Dimensions==3 and copyOldTimeStepAndScaleWithTimeStepSize ) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        ::exahype2::fd::internal::copySolution_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchIndex,
          gridCellIndex3d(x,y,z),
          unknown,
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
    }
    else if (Dimensions==3 and not copyOldTimeStepAndScaleWithTimeStepSize){
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
      {
        for (int unknown=0; unknown<unknowns; unknown++) {
          ::exahype2::fd::internal::clearSolution_LoopBody(
            patchIndex,
            gridCellIndex3d(x,y,z),
            unknown,
            patchData.QOut[patchIndex],
            QOutEnumerator
          );
        }
      }
    }
    else if (Dimensions==2 and copyOldTimeStepAndScaleWithTimeStepSize) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        ::exahype2::fd::internal::copySolution_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchIndex,
          gridCellIndex2d(x,y),
          unknown,
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
    }
    else {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++){
        for (int unknown=0; unknown<unknowns; unknown++) {
          ::exahype2::fd::internal::clearSolution_LoopBody(
            patchIndex,
            gridCellIndex2d(x,y),
            unknown,
            patchData.QOut[patchIndex],
            QOutEnumerator
          );
        }
      }
    }

    // STEP 1.5: Calculate and update auxiliary variables
    /*
    if (Dimensions==2 and second_order_formulation) {

      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeAuxiliaryVariables_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          0, // normal
          QOutEnumerator
        );
      }

      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeAuxiliaryVariables_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          1, // normal
          QOutEnumerator
        );
      }

    } else if (Dimensions==3 and second_order_formulation) {

        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            0, // normal
            QOutEnumerator
          );
        }

        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            1, // normal
            QOutEnumerator
          );
        }

        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            2, // normal
            QOutEnumerator
          );
        }

    }
    */


    // Step 2: Add source terms
    if (Dimensions==2 and evaluateAlgebraicSource) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        ::exahype2::fd::internal::addAlgebraicSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          AlgebraicSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
    }
    else if (Dimensions==3 and evaluateAlgebraicSource) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        ::exahype2::fd::internal::addAlgebraicSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          AlgebraicSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
    }

    // STEP 3: Evaluate flux and update solution accordingly
    if (evaluateFlux and Dimensions==2) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeFlux_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          flux,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          0, // normal
          tempFluxX,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeFlux_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          flux,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          1, // normal
          tempFluxY,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        // @todo 2d missing
        internal::updateSolutionWithFlux_LoopBody(
          tempFluxX, tempFluxY, tempFluxZ,
          QOutEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          unknown,
          timeStepSize,
          *(patchData.QOut + patchIndex),
          QOutEnumerator
        );
      }
    }
    else if (evaluateFlux and Dimensions==3) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeFlux_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          flux,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          0, // normal
          tempFluxX,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeFlux_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          flux,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          1, // normal
          tempFluxY,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeFlux_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          flux,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          2, // normal
          tempFluxZ,
          QOutEnumerator
        );
      }

      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        internal::updateSolutionWithFlux_LoopBody(
          tempFluxX, tempFluxY, tempFluxZ,
          QOutEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          unknown,
          timeStepSize,
          *(patchData.QOut + patchIndex),
          QOutEnumerator
        );
      }
    }

    // STEP 4: Evaluate differential source and update solution accordingly
    if (evaluateDifferentialSource and Dimensions==2) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeDifferentialSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          DifferentialSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          0, // normal
          tempDifferentialSourceX,
          QOutEnumerator,
          variant
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeDifferentialSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          DifferentialSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          1, // normal
          tempDifferentialSourceY,
          QOutEnumerator,
          variant
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        internal::updateSolutionWithDifferentialSourceTerm_LoopBody(
          tempDifferentialSourceX, tempDifferentialSourceY, tempDifferentialSourceZ,
          QOutEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          unknown,
          timeStepSize,
          *(patchData.QOut + patchIndex),
          QOutEnumerator
        );
      }
    }
    else if (evaluateDifferentialSource and Dimensions==3) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeDifferentialSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          DifferentialSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          0, // normal
          tempDifferentialSourceX,
          QOutEnumerator,
          variant
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeDifferentialSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          DifferentialSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          1, // normal
          tempDifferentialSourceY,
          QOutEnumerator,
          variant
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
        internal::computeDifferentialSourceTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          DifferentialSource,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          patchData.t[patchIndex],
          timeStepSize,
          2, // normal
          tempDifferentialSourceZ,
          QOutEnumerator,
          variant
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        internal::updateSolutionWithDifferentialSourceTerm_LoopBody(
          tempDifferentialSourceX, tempDifferentialSourceY, tempDifferentialSourceZ,
          QOutEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex3d(x,y,z),
          unknown,
          timeStepSize,
          *(patchData.QOut + patchIndex),
          QOutEnumerator
        );
      }
    }

    // STEP 5: compute the Kreiss Oliger dissipation and add it to the solution
    if (evaluateKODissipation and Dimensions==2) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeKreissOligerDissipationTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          0, // normal
          tempKODissipationX,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeKreissOligerDissipationTerm_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          patchData.t[patchIndex],
          timeStepSize,
          1, // normal
          tempKODissipationY,
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        internal::updateSolutionWithKODissipationTerm_LoopBody(
          KOSigma, tempKODissipationX, tempKODissipationY, tempKODissipationZ,
          QOutEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          unknown,
          timeStepSize,
          *(patchData.QOut + patchIndex),
          QOutEnumerator
        );
      }
    }
    else if (evaluateKODissipation and Dimensions==3) {
        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeKreissOligerDissipationTerm_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            patchData.t[patchIndex],
            timeStepSize,
            0, // normal
            tempKODissipationX,
            QOutEnumerator
          );
        }
        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeKreissOligerDissipationTerm_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            patchData.t[patchIndex],
            timeStepSize,
            1, // normal
            tempKODissipationY,
            QOutEnumerator
          );
        }
        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeKreissOligerDissipationTerm_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            patchData.t[patchIndex],
            timeStepSize,
            2, // normal
            tempKODissipationZ,
            QOutEnumerator
          );
        }
        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
        for (int unknown=0; unknown<unknowns; unknown++) {
          internal::updateSolutionWithKODissipationTerm_LoopBody(
            KOSigma, tempKODissipationX, tempKODissipationY, tempKODissipationZ,
            QOutEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            unknown,
            timeStepSize,
            *(patchData.QOut + patchIndex),
            QOutEnumerator
          );
        }
    }
  }

  if (tempFluxX!=nullptr) tarch::freeMemory(tempFluxX, tarch::MemoryLocation::Heap);
  if (tempFluxY!=nullptr) tarch::freeMemory(tempFluxY, tarch::MemoryLocation::Heap);
  if (tempFluxZ!=nullptr) tarch::freeMemory(tempFluxZ, tarch::MemoryLocation::Heap);
  if (tempDifferentialSourceX!=nullptr) tarch::freeMemory(tempDifferentialSourceX, tarch::MemoryLocation::Heap);
  if (tempDifferentialSourceY!=nullptr) tarch::freeMemory(tempDifferentialSourceY, tarch::MemoryLocation::Heap);
  if (tempDifferentialSourceZ!=nullptr) tarch::freeMemory(tempDifferentialSourceZ, tarch::MemoryLocation::Heap);
  if (tempKODissipationX!=nullptr)      tarch::freeMemory(tempKODissipationX,      tarch::MemoryLocation::Heap);
  if (tempKODissipationY!=nullptr)      tarch::freeMemory(tempKODissipationY,      tarch::MemoryLocation::Heap);
  if (tempKODissipationZ!=nullptr)      tarch::freeMemory(tempKODissipationZ,      tarch::MemoryLocation::Heap);
}

void exahype2::fd::fd4::reconstruct_first_derivatives(
  ::exahype2::CellData&   patchData,
  int                     numberOfGridCellsPerPatchPerAxis,
  int                     haloSize,
  int                     unknowns,
  int                     auxiliaryVariables
) {
  static tarch::logging::Log _log( "exahype2::fd" );

  assertionEquals(haloSize,3);

  exahype2::enumerator::AoSLexicographicEnumerator QInEnumerator (1,numberOfGridCellsPerPatchPerAxis,haloSize,unknowns,auxiliaryVariables);
  exahype2::enumerator::AoSLexicographicEnumerator QOutEnumerator(1,numberOfGridCellsPerPatchPerAxis,0,       unknowns,0);
  exahype2::enumerator::AoSLexicographicEnumerator QOutEnumeratorWithAuxiliary(1,numberOfGridCellsPerPatchPerAxis,0,unknowns,auxiliaryVariables);

  for (int patchIndex=0; patchIndex<patchData.numberOfCells; patchIndex++) {
    if (Dimensions==2) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeAuxiliaryVariables_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          0, // normal
          QOutEnumerator
        );
      }
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++) {
        internal::computeAuxiliaryVariables_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          gridCellIndex2d(x,y),
          1, // normal
          QOutEnumerator
        );
      }
    } else if (Dimensions==3) {
        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            0, // normal
            QOutEnumerator
          );
        }

        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            1, // normal
            QOutEnumerator
          );
        }

        for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
        for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
        for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++) {
          internal::computeAuxiliaryVariables_LoopBody(
            patchData.QIn[patchIndex],
            QInEnumerator,
            patchData.cellCentre[patchIndex],
            patchData.cellSize[patchIndex],
            patchIndex,
            gridCellIndex3d(x,y,z),
            2, // normal
            QOutEnumerator
          );
        }
    }

    if ( Dimensions==3 ) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int z = 0; z < numberOfGridCellsPerPatchPerAxis; z++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        ::exahype2::fd::internal::copySolution_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchIndex,
          gridCellIndex3d(x,y,z),
          unknown,
          patchData.QOut[patchIndex],
          QOutEnumeratorWithAuxiliary
        );
      }
    } else if (Dimensions==2) {
      for (int x = 0; x < numberOfGridCellsPerPatchPerAxis; x++)
      for (int y = 0; y < numberOfGridCellsPerPatchPerAxis; y++)
      for (int unknown=0; unknown<unknowns; unknown++) {
        ::exahype2::fd::internal::copySolution_LoopBody(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchIndex,
          gridCellIndex2d(x,y),
          unknown,
          patchData.QOut[patchIndex],
          QOutEnumeratorWithAuxiliary
        );
      }
    }

  }
}
