// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

namespace tarch {
  namespace plotter {
    /**
     * Transform inVal vom little-endian into big-endian.
     *
     * Paraview always requires all data to be stored as big-endian. However,
     * many architectures used by Peano, e.g., are little-endian ones. Intel
     * x86 architecture for example always is little-endian.
     */
    template <class InType>
    InType byteSwap(const InType inVal);

    /**
     * Conditional swap
     *
     * This operation forwards to byteSwap if the system is a little-endian
     * one. Otherwise, it reduces to identity.
     */
    template <class InType>
    InType byteSwapForParaviewBinaryFiles(const InType inVal);
  }
}

#include "tarch/plotter/ByteSwap.cpph"
