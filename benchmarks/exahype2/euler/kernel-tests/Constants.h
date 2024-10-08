// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#pragma once

#include <limits>
#include <string>

#include <bitset>
#include "tarch/la/Vector.h" 
#include "peano4/utils/Globals.h" 


namespace benchmarks{
namespace exahype2{
namespace kernelbenchmarks{

  const tarch::la::Vector<Dimensions,double> DomainOffset = {0.0,0.0};
  const tarch::la::Vector<Dimensions,double> DomainSize = {1.0,1.0};
  constexpr auto MinTerminalTime = 0.1;
  constexpr auto MaxTerminalTime = 0.1;
  constexpr auto FirstPlotTimeStamp = 0.0;
  constexpr auto TimeInBetweenPlots = 0.0;
  constexpr auto PlotterPrecision = 5;
  const std::bitset<2> PeriodicBC = 0;
  const std::string BuildInformation = "python3  /home/zhouxing/Documents/Masterarbeit/Peano/benchmarks/exahype2/euler/kernel-benchmarks/kernel-benchmarks-fv-rusanov.py -d 2 -ps 8 -p 512";
  const std::string ConfigureInformation = "$ ./configure --with-multithreading=omp --with-gpu=omp --enable-exahype --enable-loadbalancing --enable-blockstructured CC=clang CXX=clang++ 'CXXFLAGS=-std=c++17 -O0 -g -w -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89' 'LDFLAGS=-O0 -g -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89'";
  constexpr double Accuracy = 0.0;
  constexpr int NumberOfSamples = 10;
  const ::tarch::la::Vector<1, int> NumberOfPatchesToStudy = {512};
  constexpr int NumberOfLaunchingThreads = 1;
  constexpr bool EnableFPE = false;
  constexpr bool EvaluateFlux = true;
  constexpr bool EvaluateNonconservativeProduct = false;
  constexpr bool EvaluateSource = false;
  constexpr bool EvaluateMaximumEigenvalueAfterTimeStep = false;
  constexpr bool EvaluateHostKernels = true;
  constexpr bool EvaluateDeviceKernels = true;
  #define GAMMA 1.0

}
}
}

