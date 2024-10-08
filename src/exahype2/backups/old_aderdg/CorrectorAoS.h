// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "config.h"

#include "tarch/la/Vector.h"
#include "tarch/accelerator/accelerator.h"

#include <functional>


namespace exahype2 {
  namespace aderdg {
    
    /**
     * @brief Allow user to modify solution at the given coordinates and time and
     * compute the max eigenvalue afterwards.
     * 
     * Writes to a single spatial degree of freedom.
     *
     * @param[in] adjustSolution
     * @param[inout] UOut
     * @param[in] nodes quadrature nodes; size: (order+1)
     * @param[in] cellCentre
     * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
     * @param[in] t time stamp
     * @param[in] nodesPerAxis nodes/Lagrange basis functions per coordinate axis (order+1)
     * @param[in] unknowns the number of PDE unknowns that we evolve
     * @param[in] strideQ
     * @param[in] scalarIndex
     */ 
    // geht net durch NVC++ muessen wir angucken, sobald es geht
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp declare target
    #endif
*/
    double corrector_adjustSolution_computeMaxEigenvalue_body_AoS(
        std::function< void(
          double * __restrict__                       Q,
          const tarch::la::Vector<Dimensions,double>& x,
          double                                      t
        ) >                                         adjustSolution,
        std::function< double(
          const double * const __restrict__           Q,
          const tarch::la::Vector<Dimensions,double>& x,
          double                                      t,
          const int                                   normal
        ) >                                         maxAbsoluteEigenvalue,
        double * __restrict__                       UOut,
        const double * const __restrict__           nodes,
        const tarch::la::Vector<Dimensions,double>& cellCentre,
        const double                                dx,
        const double                                t,
        const int                                   nodesPerAxis,
        const int                                   unknowns,
        const int                                   strideQ,
        const int                                   callMaxEigenvalue,
        const int                                   scalarIndex);
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp end declare target
    #endif
*/

    /**
     * @brief Add space-time volume flux contributions to the solution.
     * 
     * Writes to a single spatial degree of freedom.
     *
     * @note Directly run after Picard iterations in order to not store predictor.
     *
     * @param[in] flux
     * @param[inout] UOut
     * @param[in] QIn
     * @param[inout] FAux
     * @param[in] nodes quadrature nodes; size: (order+1)
     * @param[in] weights quadrature weights; size: (order+1)
     * @param[in] Kxi stiffness matrix; size: (order+1)*(order+1)
     * @param[in] cellCentre
     * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
     * @param[in] t time stamp
     * @param[in] dt time step size
     * @param[in] nodesPerAxis nodes/Lagrange basis functions per coordinate axis (order+1)
     * @param[in] unknowns the number of PDE unknowns that we evolve
     * @param[in] strideQ
     * @param[in] scalarIndex
     */
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp declare target
    #endif
*/
    void corrector_addFluxContributions_body_AoS(
        std::function< void(
          const double * const __restrict__                 Q,
          const tarch::la::Vector<Dimensions,double>& x,
          double                                      t,
          int                                         normal,
          double * __restrict__                       F
        ) >                                         flux,
        double* __restrict__                        UOut, 
        const double* __restrict__                  QIn,
        double* __restrict__                        FAux, // must be allocated per thread as size is runtime parameter
        const double* __restrict__                  nodes,
        const double* __restrict__                  weights,
        const double* __restrict__                  Kxi,
        const tarch::la::Vector<Dimensions,double>& cellCentre,
        const double                                dx,
        const double                                t,
        const double                                dt,
        const int                                   nodesPerAxis,
        const int                                   unknowns,
        const int                                   strideQ,
        const int                                   scalarIndex);
/*
     #if defined(GPUOffloadingOMP)
     #pragma omp end declare target
     #endif
*/
    
    /**
     * @brief Add source contributions to the solution.
     *
     * Writes to a single spatial degree of freedom.
     *
     * @note Directly run after Picard iterations in order to not store predictor.
     *
     * @param[in] algebraicSource
     * @param[inout] UOut
     * @param[inout] SAux
     * @param[in] QIn
     * @param[in] nodes quadrature nodes; size: (order+1)
     * @param[in] weights quadrature weights; size: (order+1)
     * @param[in] cellCentre
     * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
     * @param[in] t time stamp
     * @param[in] dt time step size
     * @param[in] nodesPerAxis nodes/Lagrange basis functions per coordinate axis (order+1)
     * @param[in] unknowns the number of PDE unknowns that we evolve
     * @param[in] strideQ
     * @param[in] scalarIndex
     */
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp declare target
    #endif
*/
    void corrector_addSourceContributions_body_AoS(
        std::function< void(
          const double * const __restrict__                 Q,
          const tarch::la::Vector<Dimensions,double>& x,
          double                                      t,
          double * __restrict__                       S
        ) >                                         algebraicSource,
        double* __restrict__                        UOut,
        double* __restrict__                        SAux,
        const double* __restrict__                  QIn,
        const double* __restrict__                  nodes,
        const double* __restrict__                  weights,
        const tarch::la::Vector<Dimensions,double>& cellCentre,
        const double                                dx,
        const double                                t,
        const double                                dt,
        const int                                   nodesPerAxis,
        const int                                   unknowns,
        const int                                   strideQ,
        const int                                   scalarIndex);
/*
     #if defined(GPUOffloadingOMP)
     #pragma omp end declare target
     #endif
*/
    
     /**
      * @brief Add nonconservative product contributions to the solution.
      *
      * Writes to a single spatial degree of freedom.
      *
      * @note Directly run after Picard iterations in order to not store predictor.
      *
      * @param[in] nonconservativeProduct
      * @param[inout] UOut
      * @param[in] QIn
      * @param[inout] gradQAux
      * @param[inout] SAux
      * @param[in] nodes quadrature nodes; size: (order+1)
      * @param[in] weights quadrature weights; size: (order+1)
      * @param[in] dudx derivative operator; size: (order+1)*(order+1)
      * @param[in] cellCentre
      * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
      * @param[in] t time stamp
      * @param[in] dt time step size
      * @param[in] nodesPerAxis nodes/Lagrange basis functions per coordinate axis (order+1)
      * @param[in] unknowns the number of PDE unknowns that we evolve
      * @param[in] strideQ
      * @param[in] scalarIndex
      */
/*
     #if defined(GPUOffloadingOMP)
     #pragma omp declare target
     #endif
*/
     void corrector_addNcpContributions_body_AoS(
       std::function< void(
         const double * const __restrict__                 Q,
         double * __restrict__                       dQ_or_deltaQ,
         const tarch::la::Vector<Dimensions,double>& x,
         double                                      t,
         int                                         normal,
         double * __restrict__                       BgradQ
       ) >                                         nonconservativeProduct,
       double* __restrict__                        UOut,
       double* __restrict__                        gradQAux,
       double* __restrict__                        SAux,
       const double* __restrict__                  QIn,
       const double* __restrict__                  nodes,
       const double* __restrict__                  weights,
       const double* __restrict__                  dudx,
       const tarch::la::Vector<Dimensions,double>& cellCentre,
       const double                                dx,
       const double                                t,
       const double                                dt,
       const int                                   nodesPerAxis,
       const int                                   unknowns,
       const int                                   strideQ,
       const int                                   scalarIndex);
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp end declare target
    #endif
*/
    
    /**
     * @brief Add Riemann flux contributions to the solution.
     * 
     * Writes to a single spatial degree of freedom.
     *
     * @param[inout] UOut
     * @param[in] riemannResultIn
     * @param[in] weights quadrature weights; size: (order+1)
     * @param[in] FCoeff
     * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
     * @param[in] nodesPerAxis nodes/Lagrange basis functions per coordinate axis (order+1)
     * @param[in] unknowns the number of PDE unknowns that we evolve
     * @param[in] strideQ
     * @param[in] scalarIndex
     */
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp declare target
    #endif
*/
    GPUCallableMethod void corrector_addRiemannContributions_body_AoS(
       double * __restrict__             UOut,
       const double * const __restrict__ riemannResultIn,
       const double * const __restrict__ weights,
       const double * const __restrict__ FLCoeff,
       const double                      dx,
       const double                      dt,
       const int                         nodesPerAxis,
       const int                         unknowns,
       const int                         strideQ,
       const int                         strideF,
       const int                         scalarIndexCell);
/*
    #if defined(GPUOffloadingOMP)
    #pragma omp end declare target
    #endif
*/
    
    // CPU launchers
    /**
      * @brief Add cell-local contributions to new solution or update vector.
      *
      * @param[in] flux
      * @param[in] algebraicSource
      * @param[in] nonconservativeProduct
      * @param[inout] UOut
      * @param[in] QIn
      * @param[in] weights quadrature weights; size: (order+1)
      * @param[in] nodes quadrature nodes; size: (order+1)
      * @param[in] Kxi stiffness matrix; size: (order+1)*(order+1)
      * @param[in] dudx derivative operator; size: (order+1)*(order+1)
      * @param[in] cellCentre
      * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
      * @param[in] t time stamp
      * @param[in] dt time step size
      * @param[in] order the DG approximation order, which corresponds to order+1 DG nodes/Lagrange basis functions per coordinate axis
      * @param[in] unknowns the number of PDE unknowns that we evolve
      * @param[in] auxiliaryVariables other quantities such as material parameters that we do not evolve
      * @param[in] callFlux
      * @param[in] callSource
      * @param[in] callNonconservativeProduct
      */
    void corrector_addCellContributions_loop_AoS(
      std::function< void(
        const double * const __restrict__           Q,
        const tarch::la::Vector<Dimensions,double>& x,
        double                                      t,
        int                                         normal,
        double * __restrict__                       F
      ) >   flux,
      std::function< void(
        const double * const __restrict__           Q,
        const tarch::la::Vector<Dimensions,double>& x,
        double                                      t,
        double * __restrict__                       S
      ) >   algebraicSource,
      std::function< void(
        const double * const __restrict__           Q,
        double * __restrict__                       dQ_or_deltaQ,
        const tarch::la::Vector<Dimensions,double>& x,
        double                                      t,
        int                                         normal,
        double * __restrict__                       BgradQ
      ) >                                         nonconservativeProduct,
      double * __restrict__                       UOut, 
      const double * const __restrict__           QIn, 
      const double * const __restrict__           nodes,
      const double * const __restrict__           weights,
      const double * const __restrict__           Kxi,
      const double * const __restrict__           dudx, 
      const tarch::la::Vector<Dimensions,double>& cellCentre,
      const double                                dx,
      const double                                t,
      const double                                dt,
      const int                                   order,
      const int                                   unknowns,
      const int                                   auxiliaryVariables,
      const bool                                  callFlux,
      const bool                                  callSource,
      const bool                                  callNonconservativeProduct);
   
    /** 
     * @brief Add cell-local contributions to new solution or update vector.
     *
     * @note This routine must be called AFTER corrector_addCellContributions_loop_AoS !
     *       Otherwise calling adjustSolution and maxAbsoluteEigenvalue will not have the intended effect.
     * 
     * @return maximum eigenvalue (absolute value) computed at the quadrature nodes in the interior of the cell
     *         after the full DG update and after all solution adjustments.
     *
     * @param[inout] UOut
     * @param[in] riemannResultIn
     * @param[in] weights quadrature weights; size: (order+1)
     * @param[in] FLCoeff values of basis functions evaluated at x=0.0 (left); size: order+1
     * @param[in] dx cell spacing (we assume the same spacing in all coordinate directions)
     * @param[in] dt time step size
     * @param[in] order the DG approximation order, which corresponds to order+1 DG nodes/Lagrange basis functions per coordinate axis
     * @param[in] unknowns the number of PDE unknowns that we evolve
     * @param[in] auxiliaryVariables other quantities such as material parameters that we do not evolve
     */
     double corrector_addRiemannContributions_loop_AoS(
       std::function< void(
         double * __restrict__                       Q,
         const tarch::la::Vector<Dimensions,double>& x,
         double                                      t
       ) >                                         adjustSolution,
       std::function< double(
         const double * const __restrict__           Q,
         const tarch::la::Vector<Dimensions,double>& x,
         double                                      t,
         const int                                   normal
       ) >                                         maxAbsoluteEigenvalue,
       double * __restrict__                       UOut,
       const double * const __restrict__           riemannResultIn,
       const double * const __restrict__           nodes,
       const double * const __restrict__           weights,
       const double * const __restrict__           FLCoeff,
       const tarch::la::Vector<Dimensions,double>& cellCentre,
       const double                                dx,
       const double                                t,
       const double                                dt,
       const int                                   order,
       const int                                   unknowns,
       const int                                   auxiliaryVariables,
       const bool                                  callMaxEigenvalue);
  
  }
}

