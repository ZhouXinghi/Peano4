#pragma once

#include "GPUCellData.h"
#include "GPUCellDataPacked.h"
#include "GPUCellDataAsync.h"


namespace exahype2::fv::rusanov::omp {

template <
    class SolverType,
    int         NumberOfVolumesPerAxisInPatch,
    int         HaloSize,
    int         NumberOfUnknowns,
    int         NumberOfAuxiliaryVariables,
    bool        EvaluateFlux,
    bool        EvaluateNonconservativeProduct,
    bool        EvaluateSource,
    bool        EvaluateMaximumEigenvalueAfterTimeStep,
    class TempDataEnumeratorType,
    int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessOneHugeBuffer(int targetDevice, CellData& patchData) InlineMethod;

template <
    class SolverType,
    int         NumberOfVolumesPerAxisInPatch,
    int         HaloSize,
    int         NumberOfUnknowns,
    int         NumberOfAuxiliaryVariables,
    bool        EvaluateFlux,
    bool        EvaluateNonconservativeProduct,
    bool        EvaluateSource,
    bool        EvaluateMaximumEigenvalueAfterTimeStep,
    class TempDataEnumeratorType,
    int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessOneHugeBufferAsync(int targetDevice, CellData& patchData) InlineMethod;

template <
    class SolverType,
    int         NumberOfVolumesPerAxisInPatch,
    int         HaloSize,
    int         NumberOfUnknowns,
    int         NumberOfAuxiliaryVariables,
    bool        EvaluateFlux,
    bool        EvaluateNonconservativeProduct,
    bool        EvaluateSource,
    bool        EvaluateMaximumEigenvalueAfterTimeStep,
    class TempDataEnumeratorType,
    int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessOneHugeBufferPacked(int targetDevice, CellData& patchData) InlineMethod;

};


#include "PatchwiseStateless.cpph"
