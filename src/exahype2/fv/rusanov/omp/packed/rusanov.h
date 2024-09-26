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
        class TempDataEnumeratorType>
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
        class TempDataEnumeratorType>
      KeywordToAvoidDuplicateSymbolsForInlinedFunctions void timeStepWithRusanovPatchwiseHeapStatelessAllPacked(int targetDevice, CellData& patchData) InlineMethod;

}// namespace exahype2::fv::rusanov::omp
 
namespace exahype2::fv::rusanov::omp {
    class PackedDouble {
    public:
        [[clang::truncate_mantissa(20)]]
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
};

#include "KernelTest.cpph"
#include "BatchedStatelessPacked.cpph"
#include "PatchwiseStatelessPacked.cpph"
#include "PatchwiseStatelessAllPacked.cpph"
