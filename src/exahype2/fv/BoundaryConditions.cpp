// This file is part of the ExaHyPE2 project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#include "BoundaryConditions.h"

#include "peano4/utils/Globals.h"
#include "peano4/utils/Loop.h"
#include "PatchUtils.h"


void exahype2::fv::applyBoundaryConditions(
  std::function< void(
    const double* __restrict__                   Qinside,
    double * __restrict__                        Qoutside,
    const tarch::la::Vector<Dimensions,double>&  faceCentre,
    const tarch::la::Vector<Dimensions,double>&  volumeH,
    double                                       t,
    double                                       dt,
    int                                          normal
  ) >   boundaryCondition,
  const tarch::la::Vector<Dimensions,double>&  faceCentre,
  const tarch::la::Vector<Dimensions,double>&  patchSize,
  double                                       t,
  double                                       dt,
  int                                          numberOfVolumesPerAxisInPatch,
  int                                          overlap,
  int                                          unknowns,
  int                                          faceNumber,
  double* __restrict__                         Q
) {
  static tarch::logging::Log _log( "exahype2::fv" );

  auto serialisePatchIndex = [&](tarch::la::Vector<Dimensions,int> overlapCell) {{
    int base   = 1;
    int result = 0;
    for (int d=0; d<Dimensions; d++) {{
      result += overlapCell(d) * base;
      if (d==faceNumber % Dimensions) {{
        base *= 2*overlap;
      }}
      else {{
        base *= numberOfVolumesPerAxisInPatch;
      }}
    }}
    return result;
  }};

  logTraceInWith4Arguments( "applyBoundaryConditions(...)", faceCentre, patchSize, numberOfVolumesPerAxisInPatch, faceNumber);

  tarch::la::Vector<Dimensions,double> volumeH    = exahype2::fv::getVolumeSize(patchSize, numberOfVolumesPerAxisInPatch);
  tarch::la::Vector<Dimensions,double> faceOffset = faceCentre - 0.5 * patchSize;
  faceOffset(faceNumber%Dimensions) += 0.5 * patchSize(faceNumber%Dimensions);

  dfore(volume,numberOfVolumesPerAxisInPatch,faceNumber % Dimensions,0) {
    tarch::la::Vector<Dimensions,int> insideVolume  = volume;
    tarch::la::Vector<Dimensions,int> outsideVolume = volume;
    tarch::la::Vector<Dimensions,double> x          = faceOffset + tarch::la::multiplyComponents( tarch::la::convertScalar<double>(volume)+tarch::la::Vector<Dimensions,double>(0.5), volumeH);

    x(faceNumber%Dimensions) -= 0.5 * volumeH(faceNumber%Dimensions);

    for (int layer=0; layer<overlap; layer++){
    if (faceNumber<Dimensions) {
      insideVolume(faceNumber % Dimensions)  = overlap-layer;
      outsideVolume(faceNumber % Dimensions) = overlap-1-layer;
  
      int insideVolumeSerialised  = serialisePatchIndex(insideVolume);
      int outsideVolumeSerialised = serialisePatchIndex(outsideVolume);

      logDebug(
        "applyBoundaryConditions(...)",
        insideVolume << "->" << outsideVolume << " (" << insideVolumeSerialised << "->" << outsideVolumeSerialised << "): " <<
        "(" << *(Q + insideVolumeSerialised * unknowns) << ",...)^T->" <<
        "(" << *(Q + outsideVolumeSerialised * unknowns) << ",...)^T"
      );

      boundaryCondition(
        Q + insideVolumeSerialised * unknowns,
        Q + outsideVolumeSerialised * unknowns,
        x, volumeH, t, dt, faceNumber
      );
    }
    else {
      insideVolume(faceNumber % Dimensions)  = layer+overlap-1;
      outsideVolume(faceNumber % Dimensions) = layer+overlap;

      int insideVolumeSerialised  = serialisePatchIndex(insideVolume);
      int outsideVolumeSerialised = serialisePatchIndex(outsideVolume);

      logDebug(
        "applyBoundaryConditions(...)",
        insideVolume << "->" << outsideVolume << " (" << insideVolumeSerialised << "->" << outsideVolumeSerialised << "): " <<
        "(" << *(Q + insideVolumeSerialised * unknowns) << ",...)^T->" <<
        "(" << *(Q + outsideVolumeSerialised * unknowns) << ",...)^T"
      );

      boundaryCondition(
        Q + insideVolumeSerialised * unknowns,
        Q + outsideVolumeSerialised * unknowns,
        x, volumeH, t, dt, faceNumber
      );
    }
    }
  }

  logTraceOut( "applyBoundaryConditions(...)" );
}


