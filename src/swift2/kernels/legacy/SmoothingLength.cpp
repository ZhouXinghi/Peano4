//
// This file is delibarately empty. However, it includes the header and we
// therefore can see directly when we build the libraries if the particle
// self interaction header in itself is well-defined and consistent. We do
// not have to wait for a user to include it to see any bugs or issues
// there. Whenever that happens, it is often difficult to reconstruct what
// exactly is going on.
//

#define HYDRO_DIMENSION 2
#define HYDRO_GAMMA_5_3
#define QUARTIC_SPLINE_KERNEL

#include "SmoothingLength.h"