/**
 *
 * This code is taken from the original ExaHyPE project written by
 * colleagues from the University of Trento.
 *
 */

#pragma once

#include "tarch/la/Vector.h"

#include "peano4/utils/Globals.h"

#ifdef IncludeTwoPunctures
#include "libtwopunctures/TP_bindding.h"
#endif


namespace applications {
  namespace exahype2 {
    namespace ccz4 {
      void gaugeWave(
        double * __restrict__ Q, // Q[64+0],
        const tarch::la::Vector<Dimensions,double>&  x,
        double t
      );
      void diagonal_gaugeWave(
        double * __restrict__ Q, // Q[64+0],
        const tarch::la::Vector<Dimensions,double>&  x,
        double t
      );
      void linearWave(
        double * __restrict__ Q, // Q[64+0],
        const tarch::la::Vector<Dimensions,double>&  X,
        double t
      );
      void flat(
        double * __restrict__ Q, // Q[64+0],
        const tarch::la::Vector<Dimensions,double>&  X,
        double t
      );

      #ifdef IncludeTwoPunctures
      void ApplyTwoPunctures(
        double * __restrict__ Q, // Q[64+0],
        const tarch::la::Vector<Dimensions,double>&  X,
        double t,
	TP::TwoPunctures* tp,
        bool low_res
      );
      #endif
    }
  }
}

