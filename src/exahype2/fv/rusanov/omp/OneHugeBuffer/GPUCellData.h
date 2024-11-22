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

struct GPUCellData {
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

    // QIn and QOut contain pointers that point to somewhere in the huge buffer 
    double**                                 QIn;
    double**                                 QOut;

    GPUCellData() = delete;

    // allocate memory on gpu and copy the data from cpu (hostCellData object) onto gpu 
    GPUCellData(
        const CellData& hostCellData,
        const enumerator::AoSLexicographicEnumerator& inEnumerator,
        const enumerator::AoSLexicographicEnumerator& outEnumerator,
        int _targetDevice
    );

    // Copying GPU Data is forbidden.
    GPUCellData(const GPUCellData&) = delete;
    GPUCellData& operator=(const GPUCellData&) = delete;

    // For convenience I deleted move constructor and move assignment. It can be implemented in future if needed.
    GPUCellData(GPUCellData&&) = delete;
    GPUCellData& operator=(GPUCellData&&) = delete;

    ~GPUCellData();
    

    // copy the data from gpu back to cpu into the hostCellData object and release memory on gpu
    void copyToHost(
        CellData& hostCellData, 
        const enumerator::AoSLexicographicEnumerator& outEnumerator,
        bool copyMaxEigenvalue
    ) const;
};

};

};


exahype2::fv::rusanov::omp::GPUCellData::~GPUCellData()
{
    omp_target_free(cellCentre, targetDevice);
    omp_target_free(cellSize, targetDevice);
    omp_target_free(t, targetDevice);
    omp_target_free(dt, targetDevice);
    omp_target_free(maxEigenvalue, targetDevice);
    omp_target_free(QIn, targetDevice);
    omp_target_free(QOut, targetDevice);
    omp_target_free(QInCopyInOneHugeBuffer, targetDevice);
    omp_target_free(QOutCopyInOneHugeBuffer, targetDevice);
}

exahype2::fv::rusanov::omp::GPUCellData::GPUCellData(
    const CellData& hostCellData,
    const enumerator::AoSLexicographicEnumerator& inEnumerator,
    const enumerator::AoSLexicographicEnumerator& outEnumerator,
    int _targetDevice
) : targetDevice(_targetDevice)
{
    numberOfCells = hostCellData.numberOfCells;
    cellCentre = (::tarch::la::Vector<Dimensions, double>*)omp_target_alloc(numberOfCells * sizeof(::tarch::la::Vector<Dimensions, double>), targetDevice);
    cellSize = (::tarch::la::Vector<Dimensions, double>*)omp_target_alloc(numberOfCells * sizeof(::tarch::la::Vector<Dimensions, double>), targetDevice);
    t = (double*)omp_target_alloc(numberOfCells * sizeof(double), targetDevice);
    dt = (double*)omp_target_alloc(numberOfCells * sizeof(double), targetDevice);
    maxEigenvalue = (double*)omp_target_alloc(numberOfCells * sizeof(double), targetDevice);

    int hostDevice = omp_get_initial_device();

    // TODO: make memcpy async
    omp_target_memcpy(cellCentre, hostCellData.cellCentre, numberOfCells * sizeof(::tarch::la::Vector<Dimensions, double>), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(cellSize, hostCellData.cellSize, numberOfCells  * sizeof(::tarch::la::Vector<Dimensions, double>), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(t, hostCellData.t, numberOfCells * sizeof(double), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(dt, hostCellData.dt, numberOfCells * sizeof(double), 0, 0, targetDevice, hostDevice);

    QIn = (double**)omp_target_alloc(numberOfCells * sizeof(double*), targetDevice);
    QOut = (double**)omp_target_alloc(numberOfCells * sizeof(double*), targetDevice);
    QInCopyInOneHugeBuffer = (double*)omp_target_alloc(numberOfCells * inEnumerator.size() * sizeof(double), targetDevice);
    QOutCopyInOneHugeBuffer = (double*)omp_target_alloc(numberOfCells * outEnumerator.size() *  sizeof(double), targetDevice);

    // copy the host QIn into one huge buffer
    double* tmpQIn = new double[numberOfCells * inEnumerator.size()];
    for (int i = 0; i < numberOfCells; ++i) {
        std::memcpy(&tmpQIn[i * inEnumerator.size()], hostCellData.QIn[i], inEnumerator.size() * sizeof(double));
    }

    #pragma omp target teams distribute parallel for simd firstprivate(inEnumerator, outEnumerator, QIn, QOut, QInCopyInOneHugeBuffer, QOutCopyInOneHugeBuffer) device(targetDevice)
    for (int i = 0; i < numberOfCells; ++i) {
        QIn[i] = &QInCopyInOneHugeBuffer[i * inEnumerator.size()];
        QOut[i] = &QOutCopyInOneHugeBuffer[i * outEnumerator.size()];
    }

    // TODO: make memcpy async
    omp_target_memcpy(QInCopyInOneHugeBuffer, tmpQIn, numberOfCells * inEnumerator.size() * sizeof(double), 0, 0, targetDevice, hostDevice);
    delete[] tmpQIn;
}


void exahype2::fv::rusanov::omp::GPUCellData::copyToHost(
    CellData& hostCellData, 
    const enumerator::AoSLexicographicEnumerator& outEnumerator,
    bool copyMaxEigenvalue
) const
{
    double* tmpQOut = new double[numberOfCells * outEnumerator.size()];

    int hostDevice = omp_get_initial_device();

    // TODO: make memcpy async
    omp_target_memcpy(tmpQOut, QOutCopyInOneHugeBuffer, numberOfCells * outEnumerator.size() * sizeof(double), 0, 0, hostDevice, targetDevice);
    if (copyMaxEigenvalue) {
        omp_target_memcpy(hostCellData.maxEigenvalue, maxEigenvalue, numberOfCells * sizeof(double), 0, 0, hostDevice, targetDevice);
    }

    for (int i = 0; i < numberOfCells; ++i) {
        std::memcpy(hostCellData.QOut[i], &tmpQOut[i * outEnumerator.size()], outEnumerator.size() * sizeof(double));
    }

    delete[] tmpQOut;

}


