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

struct PackedDouble {
public:
#ifdef USE_HPC_EXT
    void printInfo()
    {
        std::printf("Using hpc-ext, mantissa == 32\n");
    }
    [[clang::truncate_mantissa(32)]] // using hpc-ext
#elif defined(USE_SOURCE_TO_SOURCE_TRANSFORM)
    void printInfo()
    {
        std::printf("Using source-to-source-transform, mantissa == 20\n");
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

class ToolFunc {
public:
    static void mempack(PackedDouble* dst, const double* src, size_t count)
    {
        for (size_t i = 0; i < count; ++i) {
            dst[i]._d = src[i];
        }
    }
    static void memunpack(double* dst, const PackedDouble* src, size_t count)
    {
        for (size_t i = 0; i < count; ++i) {
            dst[i] = src[i]._d;
        }
    }
};

struct GPUCellDataPacked {
    int                                      targetDevice;

    int                                      numberOfCells;

    ::tarch::la::Vector<Dimensions, double>* cellCentre;
    ::tarch::la::Vector<Dimensions, double>* cellSize;
    double*                                  t;
    double*                                  dt;
    double*                                  maxEigenvalue;

    // Points to a huge buffer that stores all the data
    PackedDouble* QInCopyInOneHugeBufferPacked;
    PackedDouble* QOutCopyInOneHugeBufferPacked;
    double* QInCopyInOneHugeBuffer;
    double* QOutCopyInOneHugeBuffer;

    // QIn and QOut contain pointers that point to somewhere in the huge buffer 
    double**                                 QIn;
    double**                                 QOut;

    GPUCellDataPacked() = delete;

    // allocate memory on gpu and copy the data from cpu (hostCellData object) onto gpu 
    GPUCellDataPacked(
        const CellData& hostCellData,
        const enumerator::AoSLexicographicEnumerator& inEnumerator,
        const enumerator::AoSLexicographicEnumerator& outEnumerator,
        int _targetDevice
    );

    // Copying GPU Data is forbidden.
    GPUCellDataPacked(const GPUCellDataPacked&) = delete;
    GPUCellDataPacked& operator=(const GPUCellDataPacked&) = delete;

    // For convenience I deleted move constructor and move assignment. It can be implemented in future if needed.
    GPUCellDataPacked(GPUCellDataPacked&&) = delete;
    GPUCellDataPacked& operator=(GPUCellDataPacked&&) = delete;

    ~GPUCellDataPacked();
    

    // copy the data from gpu back to cpu into the hostCellData object and release memory on gpu
    void copyToHost(
        CellData& hostCellData, 
        const enumerator::AoSLexicographicEnumerator& outEnumerator,
        bool copyMaxEigenvalue
    ) const;
};

};

};


exahype2::fv::rusanov::omp::GPUCellDataPacked::~GPUCellDataPacked()
{
    omp_target_free(cellCentre, targetDevice);
    omp_target_free(cellSize, targetDevice);
    omp_target_free(t, targetDevice);
    omp_target_free(dt, targetDevice);
    omp_target_free(maxEigenvalue, targetDevice);
    omp_target_free(QIn, targetDevice);
    omp_target_free(QOut, targetDevice);
    omp_target_free(QInCopyInOneHugeBufferPacked, targetDevice);
    omp_target_free(QOutCopyInOneHugeBufferPacked, targetDevice);
    omp_target_free(QInCopyInOneHugeBuffer, targetDevice);
    omp_target_free(QOutCopyInOneHugeBuffer, targetDevice);
}

exahype2::fv::rusanov::omp::GPUCellDataPacked::GPUCellDataPacked(
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
    QInCopyInOneHugeBufferPacked = (PackedDouble*)omp_target_alloc(numberOfCells * inEnumerator.size() * sizeof(PackedDouble), targetDevice);
    QOutCopyInOneHugeBufferPacked = (PackedDouble*)omp_target_alloc(numberOfCells * outEnumerator.size() *  sizeof(PackedDouble), targetDevice);
    QInCopyInOneHugeBuffer = (double*)omp_target_alloc(numberOfCells * inEnumerator.size() * sizeof(double), targetDevice);
    QOutCopyInOneHugeBuffer = (double*)omp_target_alloc(numberOfCells * outEnumerator.size() *  sizeof(double), targetDevice);

    // copy and pack the host QIn into one huge buffer
    PackedDouble* tmpQIn = new PackedDouble[numberOfCells * inEnumerator.size()];
    #pragma omp parallel for simd
    for (int i = 0; i < numberOfCells; ++i) {
        // std::memcpy(&tmpQIn[i * inEnumerator.size()], hostCellData.QIn[i], inEnumerator.size() * sizeof(double));
        ToolFunc::mempack(&tmpQIn[i * inEnumerator.size()], hostCellData.QIn[i], inEnumerator.size());
    }

    #pragma omp target teams distribute parallel for simd firstprivate(inEnumerator, outEnumerator, QIn, QOut, QInCopyInOneHugeBufferPacked, QOutCopyInOneHugeBufferPacked)  device(targetDevice)
    for (int i = 0; i < numberOfCells; ++i) {
        QIn[i] = &QInCopyInOneHugeBuffer[i * inEnumerator.size()];
        QOut[i] = &QOutCopyInOneHugeBuffer[i * outEnumerator.size()];
    }

    // copy the packed data onto gpu
    // TODO: make memcpy async
    omp_target_memcpy(QInCopyInOneHugeBufferPacked, tmpQIn, numberOfCells * inEnumerator.size() * sizeof(PackedDouble), 0, 0, targetDevice, hostDevice);
    delete[] tmpQIn;

    // unpack QIn on gpu
    #pragma omp target teams distribute parallel for simd device(targetDevice) firstprivate(inEnumerator, QInCopyInOneHugeBuffer, QInCopyInOneHugeBufferPacked)
    for (int i = 0; i < numberOfCells * inEnumerator.size(); ++i) {
        QInCopyInOneHugeBuffer[i] = QInCopyInOneHugeBufferPacked[i]._d;
    }
}


void exahype2::fv::rusanov::omp::GPUCellDataPacked::copyToHost(
    CellData& hostCellData, 
    const enumerator::AoSLexicographicEnumerator& outEnumerator,
    bool copyMaxEigenvalue
) const 
{
    // pack QOut on gpu
    #pragma omp target teams distribute parallel for simd device(targetDevice) firstprivate(outEnumerator, QOutCopyInOneHugeBuffer, QOutCopyInOneHugeBufferPacked)
    for (int i = 0; i < numberOfCells * outEnumerator.size(); ++i) {
        QOutCopyInOneHugeBufferPacked[i]._d = QOutCopyInOneHugeBuffer[i];
    }

    PackedDouble* tmpQOut = new PackedDouble[numberOfCells * outEnumerator.size()];

    int hostDevice = omp_get_initial_device();

    // copy the packed data onto cpu
    // TODO: make memcpy async
    omp_target_memcpy(tmpQOut, QOutCopyInOneHugeBufferPacked, numberOfCells * outEnumerator.size() * sizeof(PackedDouble), 0, 0, hostDevice, targetDevice);
    if (copyMaxEigenvalue) {
        omp_target_memcpy(hostCellData.maxEigenvalue, maxEigenvalue, numberOfCells * sizeof(double), 0, 0, hostDevice, targetDevice);
    }

    // copy and unpack the QOut into seperate host buffer
    #pragma omp parallel for simd
    for (int i = 0; i < numberOfCells; ++i) {
        // std::memcpy(hostCellData.QOut[i], &tmpQOut[i * outEnumerator.size()], outEnumerator.size() * sizeof(double));
        ToolFunc::memunpack(hostCellData.QOut[i], &tmpQOut[i * outEnumerator.size()], outEnumerator.size());
    }

    delete[] tmpQOut;

}


