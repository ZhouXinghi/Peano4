// #include <fenv.h>

//#include "exahype2/CellData.h"
#include "exahype2/fv/rusanov/rusanov.h"
#include "peano4/peano.h"
#include "FVRusanovSolver.h"
#include "tarch/accelerator/accelerator.h"
#include "ticktock.h"
#include <cstdio>

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
    using ::exahype2::fv::rusanov::omp::PackedDouble; 
    PackedDouble{}.printInfo();

    // std::cout << "sizeof(double)\t" << sizeof(double) << "\n";
    // std::cout << "sizeof(PackedDouble)\t" << sizeof(PackedDouble) << "\n";
    // std::cout << "sizeof(PackedDouble*)\t" << sizeof(PackedDouble*) << "\n";
    // std::cout << "sizeof(PackedDouble{})\t" << sizeof(PackedDouble{}) << "\n";
    // PackedDouble d;
    // std::cout << "sizeof(d)\t" << sizeof(d) << "\n";
    //
    // struct Data {
    //     [[clang::pack]]
    //     bool b[8];
    // };
    // std::cout << "sizeof(Data)\t" << sizeof(Data) << "\n";
    // std::cout << "sizeof(Data{})\t" << sizeof(Data{}) << "\n";
    // Data d2;
    // std::cout << "sizeof(d2)\t" << sizeof(d2) << "\n";

    // std::cout << omp_get_num_devices() << '\n';
    // std::cout << omp_get_device_num() << '\n';
    // std::cout << omp_get_default_device() << '\n';
    // std::cout << omp_get_initial_device() << '\n';

    // CTest t;
    // t.data = 1;
    // std::cout << t.data << std::endl;

    // using ::exahype2::fv::rusanov::omp::PackedDouble; 
    // ::exahype2::fv::rusanov::omp::PackedDouble d;
    // std::cout << sizeof(d) << std::endl;
    // std::cout << sizeof(PackedDouble(0.0)) << std::endl;
    // std::cout << d << std::endl;
    // d.printInfo();

    // PackedDouble* p = new PackedDouble[1];
    // std::cout << p << "\t" << *p << std::endl;
    // std::cout << p + 1 << "\t" << *(p + 1) << std::endl;
    // std::cout << p + 2 << "\t" << *(p + 2) << std::endl;


    // constexpr int totalIterations = 20;
    constexpr int transferTimes = TotalIterations / IterationsPerTransfer; 
    //std::printf("%d\n", NumberOfInputEntries);
    //int result = ::peano4::initParallelEnvironment(&argc, &argv);
    //std::cout << result << std::endl;
    std::printf("===============Dimension = %d\n", Dimensions);
    std::printf("===============Patch-size = %d\n", FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch);
    std::printf("===============Number of InputEntries = %d\n", NumberOfInputEntries);
    std::printf("===============Number of OutputEntries = %d\n", NumberOfOutputEntries);

  const exahype2::enumerator::AoSLexicographicEnumerator QInEnumerator(1, FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch, HaloSize, FVRusanovSolver::NumberOfUnknowns, FVRusanovSolver::NumberOfAuxiliaryVariables);
  const exahype2::enumerator::AoSLexicographicEnumerator QOutEnumerator(1, FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch, 0, FVRusanovSolver::NumberOfUnknowns, FVRusanovSolver::NumberOfAuxiliaryVariables);

    std::printf("===============Number of InputEntries = %d\n", QInEnumerator.size());
    std::printf("===============Number of OutputEntries = %d\n", QOutEnumerator.size());

    for (int p = 0; p < NumberOfPatchesToStudy.size(); ++p) {
        int patches = NumberOfPatchesToStudy[p];

        std::printf("========Patches = %d\n", patches);

        // int patches = 512;
        int device = 0;
        // static constexpr int  NumberOfFiniteVolumesPerAxisPerPatch = 8;

        auto initPatchData = [&](::exahype2::CellData& patchData) {
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

        auto freePatchData = [&](::exahype2::CellData& patchData) {
            for (int i = 0; i < patches; i++) {
                ::tarch::freeMemory(patchData.QIn[i], ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
                ::tarch::freeMemory(patchData.QOut[i], ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
            }
        };

        auto printPatch = [&](const ::exahype2::CellData& patchData) {
            for (int i = 0; i < 1; ++i) {
                std::printf("%f\n", patchData.maxEigenvalue[i]);
                for (int j = 0; j < 5; ++j)
                    std::printf("%f\n", patchData.QOut[i][j]);
                std::printf("\n");
            }

            for (int i = patches - 1; i < patches; ++i) {
                std::printf("%f\n", patchData.maxEigenvalue[i]);
                for (int j = 0; j < 5; ++j)
                    std::printf("%f\n", patchData.QOut[i][j]);
                std::printf("\n");
            }
        };

        auto comparePatch = [&](const ::exahype2::CellData& patch_lhs, const ::exahype2::CellData& patch_rhs, double precision) {
            for (int i = 0; i < patches; ++i) {
                if (std::abs(patch_lhs.maxEigenvalue[i] - patch_rhs.maxEigenvalue[i]) > precision) {
                    return false;
                }
            }
            for (int i = 0; i < patches; ++i) {
                for (int j = 0; j < NumberOfOutputEntries; ++j) {
                    if (std::abs(patch_lhs.QOut[i][j] - patch_rhs.QOut[i][j]) > precision) {
                        return false;
                    }
                }
            }
            return true;
        };






        // initPatchData();
        // TICK(CPU)
        // for (int i = 0; i < TotalIterations; ++i) {
        //     ::exahype2::fv::rusanov::timeStepWithRusanovPatchwiseHeapStateless<
        //         FVRusanovSolver,
        //         FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //         HaloSize,
        //         FVRusanovSolver::NumberOfUnknowns,
        //         FVRusanovSolver::NumberOfAuxiliaryVariables,
        //         true, //EvaluateFlux,
        //         false, //EvaluateNonconservativeProduct,
        //         false, //EvaluateSource,
        //         true, //EvaluateMaximumEigenvalueAfterTimeStep,
        //         ::exahype2::enumerator::AoSLexicographicEnumerator>
        //         (patchData);
        // }
        // TOCK(CPU)
        // printPatch();
        // freePatchData();

        // This is for the program to run due to a bug. Still do not know the reason
        //
        ::exahype2::CellData patchData(patches, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
        for (int i = 0; i < 0; ++i)
            ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStateless<
                FVRusanovSolver,
                FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
                HaloSize,
                FVRusanovSolver::NumberOfUnknowns,
                FVRusanovSolver::NumberOfAuxiliaryVariables,
                true, //EvaluateFlux,
                true, //EvaluateNonconservativeProduct,
                false, //EvaluateSource,
                true, //EvaluateMaximumEigenvalueAfterTimeStep,
                ::exahype2::enumerator::AoSLexicographicEnumerator,
                IterationsPerTransfer
                >
                (device, patchData);

        ::exahype2::CellData patchData_baseline(patches, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
        initPatchData(patchData_baseline);
        TICK(GPUBaseline)
        for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
            ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStateless<
                FVRusanovSolver,
                FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
                HaloSize,
                FVRusanovSolver::NumberOfUnknowns,
                FVRusanovSolver::NumberOfAuxiliaryVariables,
                true, //EvaluateFlux,
                false, //EvaluateNonconservativeProduct,
                false, //EvaluateSource,
                true, //EvaluateMaximumEigenvalueAfterTimeStep,
                ::exahype2::enumerator::AoSLexicographicEnumerator,
                IterationsPerTransfer
                >
                (device, patchData_baseline);
        }
        TOCK(GPUBaseline)
        printPatch(patchData_baseline);


        // initPatchData();
        // TICK(GPUNoPackingOneHugeBufferAsync)
        // for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
        //     ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessOneHugeBufferAsync<
        //         FVRusanovSolver,
        //         FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //         HaloSize,
        //         FVRusanovSolver::NumberOfUnknowns,
        //         FVRusanovSolver::NumberOfAuxiliaryVariables,
        //         true, //EvaluateFlux,
        //         false, //EvaluateNonconservativeProduct,
        //         false, //EvaluateSource,
        //         true, //EvaluateMaximumEigenvalueAfterTimeStep,
        //         ::exahype2::enumerator::AoSLexicographicEnumerator,
        //         IterationsPerTransfer
        //         >
        //         (device, patchData);
        // }
        // TOCK(GPUNoPackingOneHugeBufferAsync)
        // printPatch();
        // freePatchData();

        ::exahype2::CellData patchData_OneHugeBuffer(patches, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
        initPatchData(patchData_OneHugeBuffer);
        TICK(GPUOneHugeBuffer)
        for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
            ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessOneHugeBuffer<
                FVRusanovSolver,
                FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
                HaloSize,
                FVRusanovSolver::NumberOfUnknowns,
                FVRusanovSolver::NumberOfAuxiliaryVariables,
                true, //EvaluateFlux,
                false, //EvaluateNonconservativeProduct,
                false, //EvaluateSource,
                true, //EvaluateMaximumEigenvalueAfterTimeStep,
                ::exahype2::enumerator::AoSLexicographicEnumerator,
                IterationsPerTransfer
                >
                (device, patchData_OneHugeBuffer);
        }
        TOCK(GPUOneHugeBuffer)
        printPatch(patchData_OneHugeBuffer);


        ::exahype2::CellData patchData_OneHugeBufferPacked(patches, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory, device);
        initPatchData(patchData_OneHugeBufferPacked);
        TICK(GPUOneHugeBufferPacked)
        for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
            ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessOneHugeBufferPacked<
                FVRusanovSolver,
                FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
                HaloSize,
                FVRusanovSolver::NumberOfUnknowns,
                FVRusanovSolver::NumberOfAuxiliaryVariables,
                true, //EvaluateFlux,
                false, //EvaluateNonconservativeProduct,
                false, //EvaluateSource,
                true, //EvaluateMaximumEigenvalueAfterTimeStep,
                ::exahype2::enumerator::AoSLexicographicEnumerator,
                IterationsPerTransfer
                >
                (device, patchData_OneHugeBufferPacked);
        }
        TOCK(GPUOneHugeBufferPacked)
        printPatch(patchData_OneHugeBufferPacked);
        
        if (comparePatch(patchData_baseline, patchData_OneHugeBuffer, 1e-20)) {
            std::cout << "OneHugeBuffer correct!\n";
        } else {
            std::cout << "OneHugeBuffer Error!\n";
        }
        if (comparePatch(patchData_baseline, patchData_OneHugeBufferPacked, 1e-6)) {
            std::cout << "OneHugeBufferPacked correct!\n";
        } else {
            std::cout << "OneHugeBufferPacked Error!\n";
        }

        freePatchData(patchData_baseline);
        freePatchData(patchData_OneHugeBuffer);
        freePatchData(patchData_OneHugeBufferPacked);
        

        // initPatchData();
        // TICK(PatchWiseGPUPacked)
        // for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
        //     ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessPacked<
        //         FVRusanovSolver,
        //         FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //         HaloSize,
        //         FVRusanovSolver::NumberOfUnknowns,
        //         FVRusanovSolver::NumberOfAuxiliaryVariables,
        //         true, //EvaluateFlux,
        //         false, //EvaluateNonconservativeProduct,
        //         false, //EvaluateSource,
        //         true, //EvaluateMaximumEigenvalueAfterTimeStep,
        //         ::exahype2::enumerator::AoSLexicographicEnumerator,
        //         IterationsPerTransfer
        //         >
        //         (device, patchData);
        //         // (device, patchData, timingComputeKernel);
        // }
        // TOCK(PatchWiseGPUPacked)
        // printPatch();
        // freePatchData();

        // initPatchData();
        // TICK(PatchWiseGPUPacked)
        // for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
        //     ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessPackedOnTheFly<
        //         FVRusanovSolver,
        //         FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //         HaloSize,
        //         FVRusanovSolver::NumberOfUnknowns,
        //         FVRusanovSolver::NumberOfAuxiliaryVariables,
        //         true, //EvaluateFlux,
        //         false, //EvaluateNonconservativeProduct,
        //         false, //EvaluateSource,
        //         true, //EvaluateMaximumEigenvalueAfterTimeStep,
        //         ::exahype2::enumerator::AoSLexicographicEnumerator,
        //         IterationsPerTransfer
        //         >
        //         (device, patchData);
        //         // (device, patchData, timingComputeKernel);
        // }
        // TOCK(PatchWiseGPUPacked)
        // printPatch();
        // freePatchData();


        // initPatchData();
        // TICK(PatchWiseGPUAllPacked)
        // for (int i = 0; i < TotalIterations / IterationsPerTransfer; ++i) {
        //     ::exahype2::fv::rusanov::omp::timeStepWithRusanovPatchwiseHeapStatelessAllPacked<
        //     // ::exahype2::fv::rusanov::omp::timeStepWithRusanovBatchedHeapStatelessPacked<
        //         FVRusanovSolver,
        //         FVRusanovSolver::NumberOfFiniteVolumesPerAxisPerPatch,
        //         HaloSize,
        //         FVRusanovSolver::NumberOfUnknowns,
        //         FVRusanovSolver::NumberOfAuxiliaryVariables,
        //         true, //EvaluateFlux,
        //         false, //EvaluateNonconservativeProduct,
        //         false, //EvaluateSource,
        //         true, //EvaluateMaximumEigenvalueAfterTimeStep,
        //         ::exahype2::enumerator::AoSLexicographicEnumerator,
        //         //::exahype2::enumerator::SoALexicographicEnumerator
        //         IterationsPerTransfer
        //         >
        //         (device, patchData);
        //         // (device, patchData, timingComputeKernel);
        // }
        // TOCK(PatchWiseGPUAllPacked)
        // printPatch();
        // freePatchData();


        


        printf("\n");
    }

    //::exahype2::fv::rusanov::omp::KernelTest();

    return 0;
}
