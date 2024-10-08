// This file is part of the ExaHyPE2 project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once



#include "tarch/la/Vector.h"
#include "peano4/utils/Globals.h"


#include <functional>


namespace exahype2 {
  namespace fv {
    /**
     * Apply boundary conditions. Works only for a halo size of 1.
     *
     *
     * @param faceNumber Is usually taken from marker.getSelectedFaceNumber() and
     *  is thus a number between 0 and 2d-1.
     *
     */
    void applyBoundaryConditions(
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
    );
  }
}


