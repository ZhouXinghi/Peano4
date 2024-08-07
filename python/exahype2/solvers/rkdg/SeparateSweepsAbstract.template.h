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
#include "tarch/NonCriticalAssertions.h"

#include "peano4/utils/Globals.h"

#include "Constants.h"

{% if USE_VARIABLE_SHORTCUT is sameas True %}
#include "VariableShortcuts.h"
{% endif %}

{{SOLVER_INCLUDES}}

{% for item in NAMESPACE -%}
  namespace {{ item }} {

{%- endfor %}
  class {{CLASSNAME}};

{% for item in NAMESPACE -%}
  }
{%- endfor %}

class {{NAMESPACE | join("::")}}::{{CLASSNAME}}: public ::exahype2::Solver {
  public:
    enum class SolverState {
      GridConstruction,
      GridInitialisation,
      Suspended,
      Plotting {% for item in range(RK_STEPS) -%}
      ,
      ProjectOnFacesAndComputeVolumeIntegralOfStep{{item}},
      SolveRiemannProblemAndAddToVolumeIntegralOfStep{{item}}
      {%- endfor %}
    };

    static std::string toString(SolverState);

    static constexpr int  DGOrder = {{DG_ORDER}};

    /**
     * Maximum mesh size that this solver wants/is allowed to support. Peano 4
     * has to ensure that none of its ExaHyPE 2 solvers has to work on a mesh
     * which is coarser than its MaxH value.
     */
    static constexpr double MaxAdmissibleCellH  = {{MAX_CELL_H}};

    static constexpr int    NumberOfUnknowns           = {{NUMBER_OF_UNKNOWNS}};
    static constexpr int    NumberOfAuxiliaryVariables = {{NUMBER_OF_AUXILIARY_VARIABLES}};

    /**
     * Minimum mesh size that this solver wants to support. Peano 4 tries
     * to ensure that none of its ExaHyPE 2 solvers has to work on a mesh
     * which is finer than its MinH value. However, if there are multiple
     * solvers, the interplay of these solvers might imply that some solvers
     * work with fine mesh resolutions even though they have specified a
     * coarser MinH.
     */
    static constexpr double MinAdmissibleCellH  = {{MIN_CELL_H}};

{{BASIS_DECLARATIONS | indent(4,True) }}

    {{SOLVER_CONSTANTS}}

    {{CLASSNAME}}();

    /**
     * Alias for periodic boundary conditions.
     */
    static std::bitset<Dimensions> PeriodicBC;

    double getMinTimeStamp(bool ofCurrentlyRunningGridSweep=false) const final;
    double getMaxTimeStamp(bool ofCurrentlyRunningGridSweep=false) const final;
    double getMinTimeStepSize() const final;
    double getMaxTimeStepSize() const final;

    /**
     * Evaluate refinement criterion over cell
     *
     * @param Q Vector of unknowns
     * @param x Centre of underlying cell
     * @param h Size of underlying cell
     * @param t Time
     */
    virtual ::exahype2::RefinementCommand refinementCriterion(
      const double * __restrict__                  Q,    // Q[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      const tarch::la::Vector<Dimensions,double>&  volumeCentre,
      const tarch::la::Vector<Dimensions,double>&  volumeH,
      double                                       t
    ) {% if REFINEMENT_CRITERION_IMPLEMENTATION=="<user-defined>" %}= 0{% else %} final {% endif %};

    /**
     *
     * @param x     Position of point that is to be initialised within computational domain.
     * @param h     Size of cell within which the point is placed.
     * @param index Index of point within cell. This one can be used to look up which shape
     *              function is tied to the point. You just have to scale it with h. By
     *              translating the shape function's centre to x, you can reconstruct where
     *              the shape function is placed within the computational domain.
     * @param gridIsConstructed Boolean holds if grid is completely constructed.
     *              initialCondition() is invoked both by the grid construction and the
     *              initialisation, and you might want to do different things depending on
     *              the context. For most applications, they initialise data if and only if
     *              gridIsConstructed holds, but others might also initialise data throughout
     *              the construction to guide AMR criteria.
     *
     */
    virtual void initialCondition(
      double * __restrict__                        Q,     // Q[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      const tarch::la::Vector<Dimensions,double>&  x,
      const tarch::la::Vector<Dimensions,double>&  h,
      const tarch::la::Vector<Dimensions,int>&     index,
      bool                                         gridIsConstructed
    ) {% if INITIAL_CONDITIONS_IMPLEMENTATION=="<user-defined>" %}= 0{% else %} final {% endif %};

    {% if BOUNDARY_CONDITIONS_IMPLEMENTATION=="<none>" %}
    {% else %}
    /**
     * Apply boundary conditions. You can overwrite both the inside and
     * outside values though most BCs only modify the outside ones. Please
     * note that the boundary conditions you set here are after that subject
     * to the Riemann solver, i.e. flux and eigenvalues.
     */
    virtual void boundaryConditions(
      const double * __restrict__                  Qinside, // Qinside[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      double * __restrict__                        Qoutside, // Qoutside[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      const tarch::la::Vector<Dimensions,double>&  faceCentre,
      double                                       t,
      int                                          normal
    ) {% if BOUNDARY_CONDITIONS_IMPLEMENTATION=="<user-defined>" %}= 0{% else %} final{% endif %};
    {% endif %}


    virtual void suspendSolversForOneGridSweep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void startGridConstructionStep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void finishGridConstructionStep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void startGridInitialisationStep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void finishGridInitialisationStep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void startTimeStep(
      double globalMinTimeStamp,
      double globalMaxTimeStamp,
      double globalMinTimeStepSize,
      double globalMaxTimeStepSize
    ) override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void finishTimeStep() override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void startPlottingStep(
      double globalMinTimeStamp,
      double globalMaxTimeStamp,
      double globalMinTimeStepSize,
      double globalMaxTimeStepSize
    ) override;

    /**
     * If you hook into this routine, ensure the abstract base class
     * operation is still invoked.
     */
    virtual void finishPlottingStep() override;

    /**
     * @param currentTimeStep If you set this to false, you'll get the
     *   quantity from the preceding time step. This is important for
     *   local time stepping with fixed subcycling, as they have to
     *   know the sizes from the last time step throughout the traversal,
     *   where the current patch size might just be re-evaluated.
     *
     * @return Actually observed sizes, not the admissible quantities
     */
    double getMaxCellSize(bool currentTimeStep = true) const;
    double getMinCellSize(bool currentTimeStep = true) const;

    /**
     * Within the DG context, mesh is an alias for cell.
     */
    virtual double getMaxMeshSize() const override final;
    virtual double getMinMeshSize() const override final;

    /**
     * Update the global solver state, i.e. inform the solver about some
     * updated global quantities.
     *
     * @see setTimeStepSize(double)
     */
    void update(double timeStepSize, double timeStamp, double cellSize);

    SolverState  getSolverState() const;

    /**
     * Always holds.
     */
    virtual bool mayPlot() const override;

    /**
     * This predicate is always true, as we work with a single-sweep
     * implementation, i.e. each grid sweep realises one time step.
     */
    bool isFirstGridSweepOfTimeStep() const;
    bool isLastGridSweepOfTimeStep() const;

    /**
     * Feel free to overwrite in user code, but ensure the superclass
     * implementation is still invoked, too.
     */
    virtual void startSimulation() override;

    /**
     * Feel free to overwrite in user code, but ensure the superclass
     * implementation is still invoked, too.
     */
    virtual void finishSimulation() override;

    {% if STATELESS_PDE_TERMS %}
    /**
     * By default each cell can use stateless PDE terms, i.e. the static
     * versions of the operators. Therefore, each cell can go to the GPU
     * or inline aggressively. You can disable this behaviour for selected
     * cells via this callback.
     */
    {% else %}
    /**
     * Warning: This operation has no semantics, as you have disabled
     * stateless PDE terms. The routine is here to allow you to write one
     * code version for both a GPU-enabled solver and a solver without
     * GPU support.
     */
    {% endif %}
    virtual bool cellCanUseStatelessPDETerms(
      const tarch::la::Vector<Dimensions,double>&  cellCentre,
      const tarch::la::Vector<Dimensions,double>&  cellH,
      double                                       t,
      double                                       dt
    ) const;

  protected:
    static tarch::logging::Log  _log;

    SolverState  _solverState;

    double     _minTimeStamp;
    double     _maxTimeStamp;

    double     _minTimeStampThisTimeStep;
    double     _maxTimeStampThisTimeStep;

    double     _localMinTimeStampThisTimeStep;
    double     _localMaxTimeStampThisTimeStep;

    double     _minCellH;
    double     _maxCellH;

    double     _minCellHThisTimeStep;
    double     _maxCellHThisTimeStep;

    double     _minTimeStepSize;
    double     _maxTimeStepSize;

    double     _minTimeStepSizeThisTimeStep;
    double     _maxTimeStepSizeThisTimeStep;

    int        _cellUpdates;

    tarch::multicore::BooleanSemaphore  _semaphore;

  {{ABSTRACT_SOLVER_USER_DECLARATIONS}}
};
