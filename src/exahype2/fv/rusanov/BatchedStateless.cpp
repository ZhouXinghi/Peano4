// This file is part of the ExaHyPE2 project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org

namespace exahype2::fv::rusanov::internal {
  template <
    class SolverType,
    int  NumberOfVolumesPerAxisInPatch,
    int  HaloSize,
    int  NumberOfUnknowns,
    int  NumberOfAuxiliaryVariables,
    bool EvaluateFlux,
    bool EvaluateNonconservativeProduct,
    bool EvaluateSource,
    bool EvaluateMaximumEigenvalueAfterTimeStep,
    class TempDataEnumeratorType>
  KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovBatchedStateless(
    CellData&                                      patchData,
    peano4::utils::LoopPlacement                 loopParallelism,
    const enumerator::AoSLexicographicEnumerator& QInEnumerator,
    const enumerator::AoSLexicographicEnumerator& QOutEnumerator,
    const TempDataEnumeratorType&                  fluxEnumerator,
    const TempDataEnumeratorType&                  ncpEnumerator,
    const TempDataEnumeratorType&                  eigenvalueEnumerator,
    double*                                        tempFluxX,
    double*                                        tempFluxY,
    double*                                        tempFluxZ,
    double*                                        tempNonconservativeProductX,
    double*                                        tempNonconservativeProductY,
    double*                                        tempNonconservativeProductZ,
    double*                                        tempEigenvalueX,
    double*                                        tempEigenvalueY,
    double*                                        tempEigenvalueZ
  ) InlineMethod {
    static_assert(HaloSize == 1);

    static tarch::logging::Log _log("exahype2::fv::rusanov");
    logTraceIn("timeStepWithRusanovBatchedStateless()");

    // ====================================================
    // Copy solution over and evaluate source (if required)
    // ====================================================
    if constexpr (EvaluateSource) {
      simtDPlusOneForWithSchedulerInstructions(volume, internal::rangeOverVolumesTimesPatches(NumberOfVolumesPerAxisInPatch, patchData.numberOfCells), loopParallelism) {
        const int patchIndex = volume(Dimensions);
        loopbodies::copySolutionAndAddSourceTerm<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
      endSimtDfor
    } else {
      simtDPlusTwoForWithSchedulerInstructions(
        volume,
        internal::rangeOverVolumesTimesUnknownsPlusAuxiliaryVariablesTimesPatches(
          NumberOfVolumesPerAxisInPatch, NumberOfUnknowns, NumberOfAuxiliaryVariables, patchData.numberOfCells
        ),
        loopParallelism
      ) {
        const int unknown    = volume(Dimensions);
        const int patchIndex = volume(Dimensions + 1);
        loopbodies::copySolution(
          patchData.QIn[patchIndex], QInEnumerator, patchIndex, tarch::la::slice<Dimensions>(volume, 0), unknown, patchData.QOut[patchIndex], QOutEnumerator
        );
      }
      endSimtDfor
    }

    // ====================================================
    // Compute damping due to max eigenvalue
    // ====================================================
    simtDPlusOneForWithSchedulerInstructions(
      volume, internal::rangeOverVolumesTimesPatchesPlusHaloInXDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
    ) {
      volume(0) -= HaloSize;
      const int patchIndex = volume(Dimensions);
      loopbodies::computeMaxEigenvalue<SolverType>(
        patchData.QIn[patchIndex],
        QInEnumerator,
        patchData.cellCentre[patchIndex],
        patchData.cellSize[patchIndex],
        patchIndex,
        tarch::la::slice<Dimensions>(volume, 0),
        patchData.t[patchIndex],
        patchData.dt[patchIndex],
        0,
        tempEigenvalueX,
        eigenvalueEnumerator
      );
    }
    endSimtDfor

      simtDPlusOneForWithSchedulerInstructions(
        volume, internal::rangeOverVolumesTimesPatchesPlusHaloInYDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
      ) {
      volume(1) -= HaloSize;
      const int patchIndex = volume(Dimensions);
      loopbodies::computeMaxEigenvalue<SolverType>(
        patchData.QIn[patchIndex],
        QInEnumerator,
        patchData.cellCentre[patchIndex],
        patchData.cellSize[patchIndex],
        patchIndex,
        tarch::la::slice<Dimensions>(volume, 0),
        patchData.t[patchIndex],
        patchData.dt[patchIndex],
        1,
        tempEigenvalueY,
        eigenvalueEnumerator
      );
    }
    endSimtDfor

#if Dimensions == 3
      simtDPlusOneForWithSchedulerInstructions(
        volume, internal::rangeOverVolumesTimesPatchesPlusHaloInZDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
      ) {
      volume(2) -= HaloSize;
      const int patchIndex = volume(Dimensions);
      loopbodies::computeMaxEigenvalue<SolverType>(
        patchData.QIn[patchIndex],
        QInEnumerator,
        patchData.cellCentre[patchIndex],
        patchData.cellSize[patchIndex],
        patchIndex,
        tarch::la::slice<Dimensions>(volume, 0),
        patchData.t[patchIndex],
        patchData.dt[patchIndex],
        2,
        tempEigenvalueZ,
        eigenvalueEnumerator
      );
    }
    endSimtDfor
#endif

      simtDPlusTwoForWithSchedulerInstructions(
        volume,
        internal::rangeOverVolumesTimesUnknownsPlusAuxiliaryVariablesTimesPatches(
          NumberOfVolumesPerAxisInPatch, NumberOfUnknowns, NumberOfAuxiliaryVariables, patchData.numberOfCells
        ),
        loopParallelism
      ) {
      const int unknown    = volume(Dimensions);
      const int patchIndex = volume(Dimensions + 1);
      loopbodies::updateSolutionWithEigenvalueDamping(
        patchData.QIn[patchIndex],
        QInEnumerator,
        tempEigenvalueX,
        tempEigenvalueY,
        tempEigenvalueZ,
        eigenvalueEnumerator,
        patchData.cellCentre[patchIndex],
        patchData.cellSize[patchIndex],
        patchIndex,
        tarch::la::slice<Dimensions>(volume, 0),
        unknown,
        patchData.dt[patchIndex],
        patchData.QOut[patchIndex],
        QOutEnumerator
      );
    }
    endSimtDfor

      // ====================================================
      // Normal (conservative) flux
      // ====================================================
      if constexpr (EvaluateFlux) {
      simtDPlusOneForWithSchedulerInstructions(
        volume, internal::rangeOverVolumesTimesPatchesPlusHaloInXDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
      ) {
        volume(0) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          0, // normal
          tempFluxX,
          fluxEnumerator
        );
      }
      endSimtDfor

        simtDPlusOneForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesPatchesPlusHaloInYDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
        ) {
        volume(1) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          1, // normal
          tempFluxY,
          fluxEnumerator
        );
      }
      endSimtDfor

#if Dimensions == 3
        simtDPlusOneForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesPatchesPlusHaloInZDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
        ) {
        volume(2) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          2, // normal
          tempFluxZ,
          fluxEnumerator
        );
      }
      endSimtDfor
#endif

        simtDPlusTwoForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesUnknownsTimesPatches(NumberOfVolumesPerAxisInPatch, NumberOfUnknowns, patchData.numberOfCells), loopParallelism
        ) {
        const int unknown    = volume(Dimensions);
        const int patchIndex = volume(Dimensions + 1);
        loopbodies::updateSolutionWithFlux(
          tempFluxX,
          tempFluxY,
          tempFluxZ,
          fluxEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          unknown,
          patchData.dt[patchIndex],
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
      endSimtDfor
    }

    if constexpr (EvaluateNonconservativeProduct) {
      simtDPlusOneForWithSchedulerInstructions(
        volume, internal::rangeOverVolumesTimesPatchesPlusHaloInXDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
      ) {
        volume(0) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeNonconservativeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          0, // normal
          tempNonconservativeProductX,
          ncpEnumerator
        );
      }
      endSimtDfor

        simtDPlusOneForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesPatchesPlusHaloInYDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
        ) {
        volume(1) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeNonconservativeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          1, // normal
          tempNonconservativeProductY,
          ncpEnumerator
        );
      }
      endSimtDfor

#if Dimensions == 3
        simtDPlusOneForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesPatchesPlusHaloInZDirection(NumberOfVolumesPerAxisInPatch, HaloSize, patchData.numberOfCells), loopParallelism
        ) {
        volume(2) -= HaloSize;
        const int patchIndex = volume(Dimensions);
        loopbodies::computeNonconservativeFlux<SolverType>(
          patchData.QIn[patchIndex],
          QInEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          patchData.t[patchIndex],
          patchData.dt[patchIndex],
          2, // normal
          tempNonconservativeProductZ,
          ncpEnumerator
        );
      }
      endSimtDfor
#endif

        simtDPlusTwoForWithSchedulerInstructions(
          volume, internal::rangeOverVolumesTimesUnknownsTimesPatches(NumberOfVolumesPerAxisInPatch, NumberOfUnknowns, patchData.numberOfCells), loopParallelism
        ) {
        const int unknown    = volume(Dimensions);
        const int patchIndex = volume(Dimensions + 1);
        loopbodies::updateSolutionWithNonconservativeFlux(
          tempNonconservativeProductX,
          tempNonconservativeProductY,
          tempNonconservativeProductZ,
          ncpEnumerator,
          patchData.cellCentre[patchIndex],
          patchData.cellSize[patchIndex],
          patchIndex,
          tarch::la::slice<Dimensions>(volume, 0),
          unknown,
          patchData.dt[patchIndex],
          patchData.QOut[patchIndex],
          QOutEnumerator
        );
      }
      endSimtDfor
    }

    if constexpr (EvaluateMaximumEigenvalueAfterTimeStep) {
      simtForWithSchedulerInstructions(patchIndex, patchData.numberOfCells, loopParallelism) {
        double newMaxEigenvalue = 0.0;
        dfor(volume, NumberOfVolumesPerAxisInPatch) {
          newMaxEigenvalue = std::max(
            newMaxEigenvalue,
            loopbodies::reduceMaxEigenvalue<SolverType>(
              patchData.QOut[patchIndex],
              QOutEnumerator,
              patchData.cellCentre[patchIndex],
              patchData.cellSize[patchIndex],
              patchIndex,
              volume,
              patchData.t[patchIndex],
              patchData.dt[patchIndex]
            )
          );
        }
        patchData.maxEigenvalue[patchIndex] = newMaxEigenvalue;
      }
      endSimtFor
    }

    logTraceOut("timeStepWithRusanovBatchedStateless()");
  }
} // namespace exahype2::fv::rusanov::internal


template <
  class SolverType,
  int  NumberOfVolumesPerAxisInPatch,
  int  HaloSize,
  int  NumberOfUnknowns,
  int  NumberOfAuxiliaryVariables,
  bool EvaluateFlux,
  bool EvaluateNonconservativeProduct,
  bool EvaluateSource,
  bool EvaluateMaximumEigenvalueAfterTimeStep,
  class TempDataEnumeratorType>
void exahype2::fv::rusanov::timeStepWithRusanovBatchedCallStackStateless(
  CellData& patchData, tarch::timing::Measurement& measurement, peano4::utils::LoopPlacement loopParallelism
) {
  static tarch::logging::Log _log("exahype2::fv::rusanov");
  logTraceIn("timeStepWithRusanovBatchedCallStackStateless()");

  const enumerator::AoSLexicographicEnumerator QInEnumerator(1, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const enumerator::AoSLexicographicEnumerator QOutEnumerator(1, NumberOfVolumesPerAxisInPatch, 0, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const TempDataEnumeratorType                  fluxEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  ncpEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch + 1, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  eigenvalueEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, 1, 0);

  double tempFluxX[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempFluxY[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempFluxZ[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductX[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductY[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductZ[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueX[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueY[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueZ[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));

  tarch::timing::Watch watch("exahype2::fv::rusanov", "timeStepWithRusanovBatchedCallStackStateless", false, true);
  internal::timeStepWithRusanovBatchedStateless<
    SolverType,
    NumberOfVolumesPerAxisInPatch,
    HaloSize,
    NumberOfUnknowns,
    NumberOfAuxiliaryVariables,
    EvaluateFlux,
    EvaluateNonconservativeProduct,
    EvaluateSource,
    EvaluateMaximumEigenvalueAfterTimeStep,
    TempDataEnumeratorType>(
    patchData,
    loopParallelism,
    QInEnumerator,
    QOutEnumerator,
    fluxEnumerator,
    ncpEnumerator,
    eigenvalueEnumerator,
    tempFluxX,
    tempFluxY,
    tempFluxZ,
    tempNonconservativeProductX,
    tempNonconservativeProductY,
    tempNonconservativeProductZ,
    tempEigenvalueX,
    tempEigenvalueY,
    tempEigenvalueZ
  );
  watch.stop();
  measurement.setValue(watch.getCalendarTime());

  logTraceOut("timeStepWithRusanovBatchedCallStackStateless()");
}


template <
  class SolverType,
  int  NumberOfVolumesPerAxisInPatch,
  int  HaloSize,
  int  NumberOfUnknowns,
  int  NumberOfAuxiliaryVariables,
  bool EvaluateFlux,
  bool EvaluateNonconservativeProduct,
  bool EvaluateSource,
  bool EvaluateMaximumEigenvalueAfterTimeStep,
  class TempDataEnumeratorType>
void exahype2::fv::rusanov::timeStepWithRusanovBatchedCallStackStateless(CellData& patchData, peano4::utils::LoopPlacement loopParallelism) {
  static tarch::logging::Log _log("exahype2::fv::rusanov");
  logTraceIn("timeStepWithRusanovBatchedCallStackStateless()");

  const enumerator::AoSLexicographicEnumerator QInEnumerator(1, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const enumerator::AoSLexicographicEnumerator QOutEnumerator(1, NumberOfVolumesPerAxisInPatch, 0, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const TempDataEnumeratorType                  fluxEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  ncpEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch + 1, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  eigenvalueEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, 1, 0);

  double tempFluxX[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempFluxY[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempFluxZ[fluxEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductX[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductY[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempNonconservativeProductZ[ncpEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueX[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueY[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));
  double tempEigenvalueZ[eigenvalueEnumerator.size()] __attribute__((aligned(AlignmentOnHeap)));

  internal::timeStepWithRusanovBatchedStateless<
    SolverType,
    NumberOfVolumesPerAxisInPatch,
    HaloSize,
    NumberOfUnknowns,
    NumberOfAuxiliaryVariables,
    EvaluateFlux,
    EvaluateNonconservativeProduct,
    EvaluateSource,
    EvaluateMaximumEigenvalueAfterTimeStep,
    TempDataEnumeratorType>(
    patchData,
    loopParallelism,
    QInEnumerator,
    QOutEnumerator,
    fluxEnumerator,
    ncpEnumerator,
    eigenvalueEnumerator,
    tempFluxX,
    tempFluxY,
    tempFluxZ,
    tempNonconservativeProductX,
    tempNonconservativeProductY,
    tempNonconservativeProductZ,
    tempEigenvalueX,
    tempEigenvalueY,
    tempEigenvalueZ
  );

  logTraceOut("timeStepWithRusanovBatchedCallStackStateless()");
}


template <
  class SolverType,
  int  NumberOfVolumesPerAxisInPatch,
  int  HaloSize,
  int  NumberOfUnknowns,
  int  NumberOfAuxiliaryVariables,
  bool EvaluateFlux,
  bool EvaluateNonconservativeProduct,
  bool EvaluateSource,
  bool EvaluateMaximumEigenvalueAfterTimeStep,
  class TempDataEnumeratorType>
void exahype2::fv::rusanov::timeStepWithRusanovBatchedHeapStateless(CellData& patchData, tarch::timing::Measurement& measurement, peano4::utils::LoopPlacement loopParallelism) {
  static tarch::logging::Log _log("exahype2::fv::rusanov");
  logTraceIn("timeStepWithRusanovBatchedHeapStateless()");

  const enumerator::AoSLexicographicEnumerator QInEnumerator(1, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const enumerator::AoSLexicographicEnumerator QOutEnumerator(1, NumberOfVolumesPerAxisInPatch, 0, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const TempDataEnumeratorType                  fluxEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  ncpEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch + 1, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  eigenvalueEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, 1, 0);

  double* tempFluxX                   = EvaluateFlux ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempFluxY                   = EvaluateFlux ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempFluxZ                   = (EvaluateFlux and Dimensions == 3) ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductX = EvaluateNonconservativeProduct ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductY = EvaluateNonconservativeProduct ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductZ = (EvaluateNonconservativeProduct and Dimensions == 3)
                                          ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap)
                                          : nullptr;
  double* tempEigenvalueX             = tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap);
  double* tempEigenvalueY             = tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap);
  double* tempEigenvalueZ             = (Dimensions == 3) ? tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;

  tarch::timing::Watch watch("exahype2::fv::rusanov", "timeStepWithRusanovBatchedHeapStateless", false, true);
  internal::timeStepWithRusanovBatchedStateless<
    SolverType,
    NumberOfVolumesPerAxisInPatch,
    HaloSize,
    NumberOfUnknowns,
    NumberOfAuxiliaryVariables,
    EvaluateFlux,
    EvaluateNonconservativeProduct,
    EvaluateSource,
    EvaluateMaximumEigenvalueAfterTimeStep,
    TempDataEnumeratorType>(
    patchData,
    loopParallelism,
    QInEnumerator,
    QOutEnumerator,
    fluxEnumerator,
    ncpEnumerator,
    eigenvalueEnumerator,
    tempFluxX,
    tempFluxY,
    tempFluxZ,
    tempNonconservativeProductX,
    tempNonconservativeProductY,
    tempNonconservativeProductZ,
    tempEigenvalueX,
    tempEigenvalueY,
    tempEigenvalueZ
  );
  watch.stop();
  measurement.setValue(watch.getCalendarTime());

  if (tempFluxX != nullptr) {
    tarch::freeMemory(tempFluxX, tarch::MemoryLocation::Heap);
  }
  if (tempFluxY != nullptr) {
    tarch::freeMemory(tempFluxY, tarch::MemoryLocation::Heap);
  }
  if (tempFluxZ != nullptr) {
    tarch::freeMemory(tempFluxZ, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductX != nullptr) {
    tarch::freeMemory(tempNonconservativeProductX, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductY != nullptr) {
    tarch::freeMemory(tempNonconservativeProductY, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductZ != nullptr) {
    tarch::freeMemory(tempNonconservativeProductZ, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueX != nullptr) {
    tarch::freeMemory(tempEigenvalueX, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueY != nullptr) {
    tarch::freeMemory(tempEigenvalueY, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueZ != nullptr) {
    tarch::freeMemory(tempEigenvalueZ, tarch::MemoryLocation::Heap);
  }

  logTraceOut("timeStepWithRusanovBatchedHeapStateless()");
}


template <
  class SolverType,
  int  NumberOfVolumesPerAxisInPatch,
  int  HaloSize,
  int  NumberOfUnknowns,
  int  NumberOfAuxiliaryVariables,
  bool EvaluateFlux,
  bool EvaluateNonconservativeProduct,
  bool EvaluateSource,
  bool EvaluateMaximumEigenvalueAfterTimeStep,
  class TempDataEnumeratorType>
void exahype2::fv::rusanov::timeStepWithRusanovBatchedHeapStateless(CellData& patchData, peano4::utils::LoopPlacement loopParallelism) {
  static tarch::logging::Log _log("exahype2::fv::rusanov");
  logTraceIn("timeStepWithRusanovBatchedHeapStateless()");

  const enumerator::AoSLexicographicEnumerator QInEnumerator(1, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const enumerator::AoSLexicographicEnumerator QOutEnumerator(1, NumberOfVolumesPerAxisInPatch, 0, NumberOfUnknowns, NumberOfAuxiliaryVariables);
  const TempDataEnumeratorType                  fluxEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  ncpEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch + 1, HaloSize, NumberOfUnknowns, 0);
  const TempDataEnumeratorType                  eigenvalueEnumerator(patchData.numberOfCells, NumberOfVolumesPerAxisInPatch, HaloSize, 1, 0);

  double* tempFluxX                   = EvaluateFlux ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempFluxY                   = EvaluateFlux ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempFluxZ                   = (EvaluateFlux and Dimensions == 3) ? tarch::allocateMemory<double>(fluxEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductX = EvaluateNonconservativeProduct ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductY = EvaluateNonconservativeProduct ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;
  double* tempNonconservativeProductZ = (EvaluateNonconservativeProduct and Dimensions == 3)
                                          ? tarch::allocateMemory<double>(ncpEnumerator.size(), tarch::MemoryLocation::Heap)
                                          : nullptr;
  double* tempEigenvalueX             = tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap);
  double* tempEigenvalueY             = tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap);
  double* tempEigenvalueZ             = (Dimensions == 3) ? tarch::allocateMemory<double>(eigenvalueEnumerator.size(), tarch::MemoryLocation::Heap) : nullptr;

  internal::timeStepWithRusanovBatchedStateless<
    SolverType,
    NumberOfVolumesPerAxisInPatch,
    HaloSize,
    NumberOfUnknowns,
    NumberOfAuxiliaryVariables,
    EvaluateFlux,
    EvaluateNonconservativeProduct,
    EvaluateSource,
    EvaluateMaximumEigenvalueAfterTimeStep,
    TempDataEnumeratorType>(
    patchData,
    loopParallelism,
    QInEnumerator,
    QOutEnumerator,
    fluxEnumerator,
    ncpEnumerator,
    eigenvalueEnumerator,
    tempFluxX,
    tempFluxY,
    tempFluxZ,
    tempNonconservativeProductX,
    tempNonconservativeProductY,
    tempNonconservativeProductZ,
    tempEigenvalueX,
    tempEigenvalueY,
    tempEigenvalueZ
  );

  if (tempFluxX != nullptr) {
    tarch::freeMemory(tempFluxX, tarch::MemoryLocation::Heap);
  }
  if (tempFluxY != nullptr) {
    tarch::freeMemory(tempFluxY, tarch::MemoryLocation::Heap);
  }
  if (tempFluxZ != nullptr) {
    tarch::freeMemory(tempFluxZ, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductX != nullptr) {
    tarch::freeMemory(tempNonconservativeProductX, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductY != nullptr) {
    tarch::freeMemory(tempNonconservativeProductY, tarch::MemoryLocation::Heap);
  }
  if (tempNonconservativeProductZ != nullptr) {
    tarch::freeMemory(tempNonconservativeProductZ, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueX != nullptr) {
    tarch::freeMemory(tempEigenvalueX, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueY != nullptr) {
    tarch::freeMemory(tempEigenvalueY, tarch::MemoryLocation::Heap);
  }
  if (tempEigenvalueZ != nullptr) {
    tarch::freeMemory(tempEigenvalueZ, tarch::MemoryLocation::Heap);
  }

  logTraceOut("timeStepWithRusanovBatchedHeapStateless()");
}