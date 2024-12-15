#pragma once

#include "exahype2/enumerator/AoSLexicographicEnumerator.h"
#include "peano4/utils/Globals.h"
#include "tarch/accelerator/accelerator.h"
#include "tarch/la/Vector.h"

#include <cstring>
#include <iostream>
#include <omp.h>

namespace exahype2 {

struct CellData;

namespace fv::rusanov::omp {

struct GPUCellDataAsync {
    int dep;
    int                                      targetDevice;

    int                                      numberOfCells;

    ::tarch::la::Vector<Dimensions, double>* cellCentre;
    ::tarch::la::Vector<Dimensions, double>* cellSize;
    double*                                  t;
    double*                                  dt;
    double*                                  maxEigenvalue;

    // Points to a huge buffer that stores all the data
    double* QInCopyInOneHugeBuffer;
    double* QOutCopyInOneHugeBuffer;

    // QIn and QOut contain pointers that point to different positions in the huge buffer 
    double**                                 QIn;
    double**                                 QOut;

    GPUCellDataAsync() = delete;

    // allocate memory on gpu and copy the data from cpu (hostCellData object) onto gpu 
    GPUCellDataAsync(
        int _targetDevice, 
        int _numberOfCells,
        ::tarch::la::Vector<Dimensions, double>* _cellCentre,
        ::tarch::la::Vector<Dimensions, double>* _cellSize,
        double*                                  _t,
        double*                                  _dt,

        double**                                 _QIn,
        double**                                 _QOut,
        double* _maxEigenvalue,

        const enumerator::AoSLexicographicEnumerator& inEnumerator,
        const enumerator::AoSLexicographicEnumerator& outEnumerator
    );

    // Copying GPU Data is forbidden.
    GPUCellDataAsync(const GPUCellDataAsync&) = delete;
    GPUCellDataAsync& operator=(const GPUCellDataAsync&) = delete;

    // For convenience I deleted move constructor and move assignment. It can be implemented in future if needed.
    GPUCellDataAsync(GPUCellDataAsync&&) = delete;
    GPUCellDataAsync& operator=(GPUCellDataAsync&&) = delete;

    ~GPUCellDataAsync();
    

    // copy the data from gpu back to cpu into the hostCellData object and release memory on gpu
    void copyToHost(
        const enumerator::AoSLexicographicEnumerator& outEnumerator,
        bool copyMaxEigenvalue
    ) const;
};

};

};


exahype2::fv::rusanov::omp::GPUCellDataAsync::~GPUCellDataAsync()
{
    #pragma omp target exit data map(delete: cellCentre[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: cellSize[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: t[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: dt[0 : numberOfCells]) device(targetDevice)
    // #pragma omp target exit data map(delete: maxEigenvalue[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: QIn[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: QOut[0 : numberOfCells]) device(targetDevice)
    #pragma omp target exit data map(delete: QInCopyInOneHugeBuffer[0 : numberOfCells]) device(targetDevice)
    // #pragma omp target exit data map(delete: QOutCopyInOneHugeBuffer[0 : numberOfCells]) device(targetDevice)
    delete[] QInCopyInOneHugeBuffer;
    delete[] QOutCopyInOneHugeBuffer;
}

exahype2::fv::rusanov::omp::GPUCellDataAsync::GPUCellDataAsync(
    int _targetDevice, 
    int _numberOfCells,
    ::tarch::la::Vector<Dimensions, double>* _cellCentre,
    ::tarch::la::Vector<Dimensions, double>* _cellSize,
    double*                                  _t,
    double*                                  _dt,

    double**                                 _QIn,
    double**                                 _QOut,
    double* _maxEigenvalue,

    const enumerator::AoSLexicographicEnumerator& inEnumerator,
    const enumerator::AoSLexicographicEnumerator& outEnumerator
) : 
targetDevice(_targetDevice),
numberOfCells(_numberOfCells)
{
    cellCentre = _cellCentre;
    cellSize = _cellSize;
    t = _t;
    dt = _dt;
    maxEigenvalue = _maxEigenvalue;
    #pragma omp target enter data map(to: cellCentre[0 : numberOfCells]) device(targetDevice)           
    #pragma omp target enter data map(to: cellSize[0 : numberOfCells]) device(targetDevice)             
    #pragma omp target enter data map(to: t[0 : numberOfCells]) device(targetDevice)                    
    #pragma omp target enter data map(to: dt[0 : numberOfCells]) device(targetDevice)                   
    #pragma omp target enter data map(alloc: maxEigenvalue[0 : numberOfCells]) device(targetDevice)     


    int hostDevice = omp_get_initial_device();

    QOut = _QOut;

    #pragma omp target enter data map(alloc: QIn[0 : numberOfCells]) device(targetDevice)               
    #pragma omp target enter data map(alloc: QOut[0 : numberOfCells]) device(targetDevice)              

    QInCopyInOneHugeBuffer = new double[numberOfCells * inEnumerator.size()];
    QOutCopyInOneHugeBuffer = new double[numberOfCells * outEnumerator.size()];

    // copy the host QIn into one huge buffer
    for (int i = 0; i < numberOfCells; ++i) {
        std::memcpy(&QInCopyInOneHugeBuffer[i * inEnumerator.size()], _QIn[i], inEnumerator.size() * sizeof(double));
    }


    #pragma omp target enter data map(to: QInCopyInOneHugeBuffer[0 : numberOfCells * inEnumerator.size()]) device(targetDevice)       
    #pragma omp target enter data map(alloc: QOutCopyInOneHugeBuffer[0 : numberOfCells * outEnumerator.size()]) device(targetDevice)  

    #pragma omp target teams distribute parallel for simd firstprivate(inEnumerator, outEnumerator, QIn, QOut, QInCopyInOneHugeBuffer, QOutCopyInOneHugeBuffer) device(targetDevice) depend(out: QInCopyInOneHugeBuffer) nowait
    for (int i = 0; i < numberOfCells; ++i) {
        QIn[i] = &QInCopyInOneHugeBuffer[i * inEnumerator.size()];
        QOut[i] = &QOutCopyInOneHugeBuffer[i * outEnumerator.size()];
    }
}


void exahype2::fv::rusanov::omp::GPUCellDataAsync::copyToHost(
    const enumerator::AoSLexicographicEnumerator& outEnumerator,
    bool copyMaxEigenvalue
) const
{
    int hostDevice = omp_get_initial_device();

    // TODO: make memcpy async
    #pragma omp target exit data map(from: QOutCopyInOneHugeBuffer[0 : numberOfCells * outEnumerator.size()]) device(targetDevice) depend(in: QOutCopyInOneHugeBuffer) nowait

    if (copyMaxEigenvalue) {
        #pragma omp target exit data map(from: maxEigenvalue[0 : numberOfCells]) device(targetDevice) depend(in: QOutCopyInOneHugeBuffer) nowait
    }

    for (int i = 0; i < numberOfCells; ++i) {
        std::memcpy(QOut[i], &QOutCopyInOneHugeBuffer[i * outEnumerator.size()], outEnumerator.size() * sizeof(double));
    }

}


