//
// ExaHyPE2 solver file
// Generated by Peano's Python API
// www.peano-framework.org
//
// This is generated. If you change fundamental properties, you will have to 
// generate this file. Backup your manual changes before you do so.
//
#ifndef _EXAMPLES_EXAHYPE2_FINITEVOLUMES_LOH1_
#define _EXAMPLES_EXAHYPE2_FINITEVOLUMES_LOH1_


#include "AbstractLOH1.h"

#include "tarch/logging/Log.h"


namespace examples{
namespace exahype2{
namespace loh1{
  class LOH1;
}}}



class examples::exahype2::loh1::LOH1: public examples::exahype2::loh1::AbstractLOH1 {
  private:
    static tarch::logging::Log   _log;
    
    static constexpr int Unknowns = 13;
    
    struct VariablesShortcuts {
      const int v = 0; 
      const int sigma = 3;
      const int rho = 9;
      const int cp = 10;
      const int cs = 11;
      const int alpha = 12;
    } s;

    void prescribeGaussianWave(
        const tarch::la::Vector<Dimensions,double>&  x,
    		double Q[]);

  public:
    virtual void adjustSolution(
      double * __restrict__                        Q, // [9+4],
      const tarch::la::Vector<Dimensions,double>&  volumeCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t,
      double                                       dt
    ) override;

    double maxEigenvalue(
      const double * __restrict__                  Q, // [9+4],
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t,
      int                                          normal
    ) override;

    void boundaryConditions(
      const double* __restrict__                   Qinside, // [13]
      double* __restrict__                         Qoutside, // [13]
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t,
      int                                          normal
    ) override;

    void nonconservativeProduct(
      const double * __restrict__                  Q, // [9+4],
      const double * __restrict__                  deltaQ, // [9+4],
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t,
      int                                          normal,
      double * __restrict__ BgradQ // BgradQ[13]
     ) override;

};


#endif