// This file is part of the ExaHyPE2 project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once



namespace exahype2::fv::rusanov::omp {

void KernelTest() InlineMethod;

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
class TempDataEnumeratorType
> void timeStepWithRusanovBatchedHeapStatelessPacked(int targetDevice, CellData& patchData, tarch::timing::Measurement& measurement) InlineMethod;

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
class TempDataEnumeratorType,
int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessPacked(int targetDevice, CellData& patchData) InlineMethod;

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
class TempDataEnumeratorType,
int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessPackedOnTheFly(int targetDevice, CellData& patchData) InlineMethod;


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
class TempDataEnumeratorType,
int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessAllPacked(int targetDevice, CellData& patchData) InlineMethod;

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
class TempDataEnumeratorType
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessAllPacked(int targetDevice, CellData& patchData) InlineMethod;

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
class TempDataEnumeratorType,
int IterationsPerTransfer
>
KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStateless(int targetDevice, CellData& patchData) InlineMethod;

}// namespace exahype2::fv::rusanov::omp
 
namespace exahype2::fv::rusanov::omp {
class PackedDouble {
public:
#ifdef USE_HPC_EXT
    void printInfo()
    {
        std::printf("Using hpc-ext, mantissa == 32");
    }
    [[clang::truncate_mantissa(32)]] // using hpc-ext
#elif defined(USE_SOURCE_TO_SOURCE_TRANSFORM)
    void printInfo()
    {
        std::printf("Using source-to-source-transform, mantissa == 20");
    }
    [[clang::truncate_mantissa(20)]] // using source-to-source-transform
#endif
    double _d;

    PackedDouble() : _d(0.0) 
    {
        //std::printf("default constructor is called, _d = %f\n", _d);
    }
    PackedDouble(double other) 
    { 
        _d = other;
        //std::printf("constructor is called, _d = %f\n", _d);
    }
    operator double() const 
    { 
        return _d; 
    }
    PackedDouble& operator=(double other) 
    {
        _d = other;
        return *this;
    }
    void print() const 
    {
        std::printf("%f\n", _d);
    }
};
}

#include "KernelTest.cpph"
#include "BatchedStatelessPacked.cpph"
#include "PatchwiseStatelessPacked.cpph"
// #include "PatchwiseStatelessPackedOnTheFly.cpph"
#include "PatchwiseStatelessAllPacked.cpph"
// #include "PatchwiseStateless.cpph"
