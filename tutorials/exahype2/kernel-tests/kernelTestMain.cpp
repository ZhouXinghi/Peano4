#include <fenv.h>
#pragma float_control(precise, on)
#pragma STDC FENV_ACCESS ON

//#include "exahype2/CellData.h"
#include "exahype2/fv/rusanov/rusanov.h"
#include "peano4/peano.h"
#include "FVRusanovSolver.h"
#include "tarch/accelerator/accelerator.h"
#include "ticktock.h"
#include <iostream>

using namespace benchmarks::exahype2::kernelbenchmarks;
double                               validMaxEigenvalue = 0.0;
double*                              validOutcome       = nullptr;
static constexpr double              TimeStamp          = 0.5;
static constexpr double              TimeStepSize       = 1e-6;
static constexpr double              CellSize           = 0.1;
static constexpr double              CellOffset         = 4.0;
static constexpr int                 HaloSize           = 1;
::tarch::timing::Measurement         timingComputeKernel;
::tarch::multicore::BooleanSemaphore validateOutcomeSemaphore;
static_assert(Accuracy >= std::numeric_limits<double>::epsilon() || Accuracy == 0.0);

#if Dimensions == 2
static constexpr int NumberOfInputEntries = (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 2 * HaloSize)
                                            * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 2 * HaloSize) * FVRusanovSolver::NumberOfUnknowns;
static constexpr int NumberOfOutputEntries = (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 0) * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 0)
                                             * FVRusanovSolver::NumberOfUnknowns;
static constexpr int NumberOfFiniteVolumesPerPatch = FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch * FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch;
#else
static constexpr int NumberOfInputEntries = (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 2 * HaloSize)
                                            * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 2 * HaloSize)
                                            * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 2 * HaloSize) * FVRusanovSolver::NumberOfUnknowns;
static constexpr int NumberOfOutputEntries = (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 0) * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 0)
                                             * (FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch + 0) * FVRusanovSolver::NumberOfUnknowns;
static constexpr int NumberOfFiniteVolumesPerPatch = FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch * FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch
                                                     * FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch;
#endif

void initInputData(double* Q) {
  for (int i = 0; i < NumberOfInputEntries; i++) {
    Q[i] = std::sin(1.0 * i / (NumberOfInputEntries) * ::tarch::la::PI);
  }
}


int main(int argc, char** argv)
{
    //std::printf("%d\n", NumberOfInputEntries);
    //int result = ::peano4::initParallelEnvironment(&argc, &argv);
    //std::cout << result << std::endl;
    int patches = 512;
    int device = 0;
    static constexpr int  NumberOfFiniteVolumesPerAxisPerPatch = 8;
    ::exahype2::CellData patchData(patches, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
    
    auto initPatchData = [&] {
        for (int i = 0; i < patches; i++) {
            patchData.QIn[i]           = ::tarch::allocateMemory<double>(NumberOfInputEntries, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
            patchData.t[i]             = TimeStamp;
            patchData.dt[i]            = TimeStepSize;
            patchData.QOut[i]          = ::tarch::allocateMemory<double>(NumberOfOutputEntries, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
            patchData.cellCentre[i]    = ::tarch::la::Vector<Dimensions, double>(CellOffset + 0.5 * CellSize);
            patchData.cellSize[i]      = ::tarch::la::Vector<Dimensions, double>(CellSize);
            patchData.maxEigenvalue[i] = 0.0;
            initInputData(patchData.QIn[i]);
            std::memset(patchData.QOut[i], 0.0, NumberOfOutputEntries * sizeof(double));
        }
    };
    
    //::exahype2::fv::rusanov::timeStepWithRusanovBatchedHeapStateless<
        //FVRusanovSolver,
        //FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //HaloSize,
        //FVRusanovSolver::NumberOfUnknowns,
        //FVRusanovSolver::NumberOfAuxiliaryVariables,
        //EvaluateFlux,
        //EvaluateNonconservativeProduct,
        //EvaluateSource,
        //EvaluateMaximumEigenvalueAfterTimeStep,
        //::exahype2::enumerator::AoSLexicographicEnumerator>(patchData, timingComputeKernel, ::peano4::utils::LoopPlacement::Serial);



    initPatchData();
    TICK(timeStepWithRusanovBatchedHeapStateless)
    for (int i = 0; i < 50; ++i) {
        ::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStateless<
        //::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStatelessPacked<
            FVRusanovSolver,
            FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
            HaloSize,
            FVRusanovSolver::NumberOfUnknowns,
            FVRusanovSolver::NumberOfAuxiliaryVariables,
            EvaluateFlux,
            EvaluateNonconservativeProduct,
            EvaluateSource,
            EvaluateMaximumEigenvalueAfterTimeStep,
            ::exahype2::enumerator::AoSLexicographicEnumerator>
            (device, patchData, timingComputeKernel);
    }
    TOCK(timeStepWithRusanovBatchedHeapStateless)

    initPatchData();
    TICK(timeStepWithRusanovBatchedHeapStatelessPacked)
    for (int i = 0; i < 50; ++i) {
        ::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStatelessPacked<
        //::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStateless<
            FVRusanovSolver,
            FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
            HaloSize,
            FVRusanovSolver::NumberOfUnknowns,
            FVRusanovSolver::NumberOfAuxiliaryVariables,
            EvaluateFlux,
            EvaluateNonconservativeProduct,
            EvaluateSource,
            EvaluateMaximumEigenvalueAfterTimeStep,
            ::exahype2::enumerator::AoSLexicographicEnumerator
            //::exahype2::enumerator::SoALexicographicEnumerator
        >
            (device, patchData, timingComputeKernel);
    }
    TOCK(timeStepWithRusanovBatchedHeapStatelessPacked)

    //PackedDouble pd;
    //pd._d = 2.0;
    //for (int i = 0; i < 10; ++i) {
        //::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStatelessPacked<
        //::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStateless<
        //::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStateless<
            //FVRusanovSolver,
            //FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
            //HaloSize,
            //FVRusanovSolver::NumberOfUnknowns,
            //FVRusanovSolver::NumberOfAuxiliaryVariables,
            //EvaluateFlux,
            //EvaluateNonconservativeProduct,
            //EvaluateSource,
            //EvaluateMaximumEigenvalueAfterTimeStep,
            ////::exahype2::enumerator::AoSLexicographicEnumerator
            //::exahype2::enumerator::SoALexicographicEnumerator
        //>
            //(device, patchData, timingComputeKernel);
    //}

        //for (int i = 0; i < 5; ++i) {
            //std::printf("%f\n", patchData.maxEigenvalue[i]);
            //for (int j = 0; j < 5; ++j)
                //std::printf("%f\n", patchData.QOut[i][j]);
            //std::printf("\n");
        //}

        //for (int i = 0; i < patchData.numberOfCells; ++i) {
            //std::printf("%f\n", patchData.maxEigenvalue[i]);
            ////for (int j = 0; j < 10; ++j)
                ////std::printf("%f\n", patchData.QOut[i][j]);
        //}
    //}

    //::exahype2::fv::rusanov::omp::KernelTest();

    return 0;
}
