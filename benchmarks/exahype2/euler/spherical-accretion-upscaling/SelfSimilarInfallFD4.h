//
// ExaHyPE2 solver file
// Generated by Peano's Python API
// www.peano-framework.org
//
// This is generated. If you change fundamental properties, you will have to 
// generate this file. Backup your manual changes before you do so.
//
#pragma once

#include "AbstractSelfSimilarInfallFD4.h"
#include "tarch/logging/Log.h"

#include "spherical-accretion/MassAccumulator.h"
#include "spherical-accretion/GravityModel.h"

#include "EulerKernels.h"

namespace benchmarks {
  namespace exahype2 {
    namespace euler {
      namespace sphericalaccretionupscaling {
        class SelfSimilarInfallFD4;
      }
    }
  }
}


class benchmarks::exahype2::euler::sphericalaccretionupscaling::SelfSimilarInfallFD4: public benchmarks::exahype2::euler::sphericalaccretionupscaling::AbstractSelfSimilarInfallFD4 {
  private:
    static tarch::logging::Log   _log;

    /**
     * Use class from application folders.
     */
    ::applications::exahype2::euler::sphericalaccretion::MassAccumulator  _accumulator;

    /**
     * After each time step, we backup the total mass in this static variable.
     * It has to be static, as we want to use it in stateless terms later on.
     */
    static double TotalMassInPreviousTimeStep;

  public:
    static constexpr double BaseDensity         = 0.1;
    static constexpr double Gamma               = 5.0/3.0;
    static constexpr double pInitial            = 1e-6;
    static constexpr double aInitial            = 0.001;
    static constexpr double AdditionalMass      = 0.15;
    static constexpr double InitialTopHatRadius = 0.05;

    SelfSimilarInfallFD4();


    /**
     * Add density of one voxel to global bookkeeping
     */
    void addDensity(
      const tarch::la::Vector<Dimensions,double>&  volumeCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       density
    );

    
    void finishTimeStep() override;


    /**
     * Refinement criterion
     *
     * ExaHypE2 is guided by a maximum and minimum mesh (patch) size.
     * All (dynamic) AMR is constrained by these values, i.e. if your
     * mesh is coarser than the maximum mesh size, ExaHyPE 2 will
     * automatically refine. If you try to refine further than the
     * minimum mesh size, ExaHyPE 2 will ignore any refinement.
     *
     * Consequently, you are fine if you work with a regular mesh:
     * You set the maximum mesh size, and you leave everything else
     * to Peano 4/ExaHyPE 2. If you want to have an adaptive mesh,
     * use this routine to implement the refinement pattern.
     *
     * @param Q This is the (current) solution. The data is not set
     *  to a valid value throughout grid construction. That is: If
     *  t equals 0.0, you cannot assume that Q is properly initialised.
     *  Therefore, Q usually is only evaluated by dynamic AMR codes
     *  which make the solution follow
     */
    ::exahype2::RefinementCommand refinementCriterion(
      const double * __restrict__ Q, // Q[5+0],
      const tarch::la::Vector<Dimensions,double>&  volumeCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t
    ) override;
    


    
    void initialCondition(
      double * __restrict__ Q,
      const tarch::la::Vector<Dimensions,double>&  volumeCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      bool                                         gridIsConstructed
    ) override;
    

    
    
    virtual void boundaryConditions(
      const double * __restrict__ Qinside, // Qinside[5+0]
      double * __restrict__ Qoutside, // Qoutside[5+0]
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t,
      int                                          normal
    )  override;
    

    
  public:
    
    void sourceTerm(
      const double * __restrict__ Q,
      const tarch::la::Vector<Dimensions,double>&  gridCellCentre,
      const tarch::la::Vector<Dimensions,double>&  gridCellH,
      double                                       t,
      double                                       dt,
      double * __restrict__ S
    ) override;
    

    
    /**
     * Determine max eigenvalue over Jacobian in a given point with solution values
     * (states) Q. All parameters are in.
     *
     * @return Max eigenvalue. Result has to be positive, so we are actually speaking
     *   about the maximum absolute eigenvalue.
     */
    virtual double maxEigenvalue(
      const double * __restrict__ Q, // Q[5+0],
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  gridCellH,
      double                                       t,
      double                                       dt,
      int                                          normal
    ) override;
    


    
    virtual void flux(
      const double * __restrict__ Q, // Q[5+0],
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  gridCellH,
      double                                       t,
      double                                       dt,
      int                                          normal,
      double * __restrict__ F // F[5]
    ) override;
    
     
     
    
    
   
    


    


    
     
   
        
};


