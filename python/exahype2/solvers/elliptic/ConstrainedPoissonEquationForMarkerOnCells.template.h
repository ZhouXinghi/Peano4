// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#pragma once


#include "exahype2/RefinementControl.h"
#include "exahype2/Solver.h"

#include "tarch/la/Vector.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include "peano4/utils/Globals.h"
#include "peano4/datamanagement/FaceEnumerator.h"

#include "Constants.h"

#include "facedata/{{MARKER}}.h"
#include "celldata/{{MARKER}}.h"


{% for item in NAMESPACE -%}
  namespace {{ item }} {

{%- endfor %}
  class {{CLASSNAME}};

{% for item in NAMESPACE -%}
  }
{%- endfor %}



class {{NAMESPACE | join("::")}}::{{CLASSNAME}}: public ::exahype2::Solver {
  private:
    static constexpr double RelaxationCoefficient = {{RELAXATION_COEFFICIENT}};
    static constexpr double Rhs                   = {{RHS}};
    static constexpr double MinValue              = {{MIN_VALUE}};
    static constexpr double MaxValue              = {{MAX_VALUE}};
    static constexpr double MinThreshold          = {{MIN_THRESHOLD}};
    static constexpr double MaxThreshold          = {{MAX_THRESHOLD}};

  public:
    /**
     * We do not impose any constraints on the patch/cell size.
     */
    static constexpr double MaxAdmissibleCellH  = std::numeric_limits<double>::max();
    static constexpr double MinAdmissibleCellH  = 0.0;


    {{CLASSNAME}}();

    /**
     * Alias for periodic boundary conditions.
     */
    static std::bitset<Dimensions> PeriodicBC;

    /**
     * Infinity, as we do not restrict anything.
     */
    double getMinTimeStamp(bool ofCurrentlyRunningGridSweep=false) const final;
    double getMaxTimeStamp(bool ofCurrentlyRunningGridSweep=false) const final;


    double getMinTimeStepSize() const final;
    double getMaxTimeStepSize() const final;

    /**
     * Nop
     */
    void startGridConstructionStep() final;
    void finishGridConstructionStep() final;
    void startGridInitialisationStep() final;
    void finishGridInitialisationStep() final;
    void startTimeStep(
      double globalMinTimeStamp,
      double globalMaxTimeStamp,
      double globalMinTimeStepSize,
      double globalMaxTimeStepSize
    ) final;
    void finishTimeStep() final;
    void startPlottingStep(
      double globalMinTimeStamp,
      double globalMaxTimeStamp,
      double globalMinTimeStepSize,
      double globalMaxTimeStepSize
    ) final;
    void finishPlottingStep() final;

    virtual double getMaxMeshSize() const override final;
    virtual double getMinMeshSize() const override final;

    virtual bool mayPlot() const override;

    virtual void startSimulation() override;
    virtual void finishSimulation() override;

    /**
     * By construction, we always are in the first and last sweep of time
     * step. Logically, each individual sweep is a time step.
     */
    bool isFirstGridSweepOfTimeStep() const;
    bool isLastGridSweepOfTimeStep() const;

    /**
     * Ensure that the value is actually within the bounds specified. Also
     * ensure that the value is in line with the marker value.
     */
    void constrainValue( benchmarks::exahype2::euler::sphericalaccretionupscaling::celldata::{{MARKER}}& marker );


    /**
     * Compute a new state for this cell. This is the state using the current
     * cell value as it got delivered by the Poisson solver. It is not yet
     * subject to any a posteriori update.
     */
    void computeNewMarkerState(
      benchmarks::exahype2::euler::sphericalaccretionupscaling::celldata::{{MARKER}}&                                                cellMarker,
      const peano4::datamanagement::FaceEnumerator<benchmarks::exahype2::euler::sphericalaccretionupscaling::facedata::{{MARKER}}>&  faceMarker
    );
  protected:
    static tarch::logging::Log  _log;
};



{# Empty line here #}