// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************

#include <iomanip>

#include "config.h"
#include "Constants.h"
#include "{{MAIN_NAME}}.h"
#include "exahype2/UnitTests.h"
#include "exahype2/UserInterface.h"
#include "observers/CreateGrid.h"
#include "observers/CreateGridAndConvergeLoadBalancing.h"
#include "observers/CreateGridButPostponeRefinement.h"
#include "observers/InitGrid.h"
#include "observers/PlotSolution.h"
#include "observers/TimeStep.h"
#include "peano4/peano.h"
#include "peano4/UnitTests.h"
#include "repositories/DataRepository.h"
#include "repositories/SolverRepository.h"
#include "repositories/StepRepository.h"
#include "tarch/UnitTests.h"
#include "peano4/grid/Spacetree.h"
#include "peano4/parallel/SpacetreeSet.h"
#include "tarch/logging/Log.h"
#include "tarch/logging/LogFilter.h"
#include "tarch/logging/Statistics.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/multicore.h"
#include "tarch/multicore/otter.h"
#include "tarch/NonCriticalAssertions.h"
#include "tarch/tests/TreeTestCaseCollection.h"
#include "tarch/timing/Measurement.h"
#include "tarch/timing/Watch.h"
#include "toolbox/blockstructured/UnitTests.h"
#include "toolbox/loadbalancing/loadbalancing.h"


{% if FENV_ARGS is defined and FENV_ARGS != "" -%}
#include <fenv.h>
#pragma float_control(precise, on)
#pragma STDC FENV_ACCESS ON
{% endif -%}

using namespace {{NAMESPACE | join("::")}};

tarch::logging::Log _log("::");


tarch::timing::Measurement timePerMeshSweepMeasurement;
tarch::timing::Measurement gridConstructionMeasurement;
tarch::timing::Measurement timeStepMeasurement;
tarch::timing::Measurement plotMeasurement;


/**
 * Decide which step to run next
 *
 * ## Control of the parallel grid construction
 *
 *
 * @return continues to run
 */
bool selectNextAlgorithmicStep() {
  static bool                                  gridConstructed                       = false;
  static bool                                  gridInitialised                       = false;
  static bool                                  gridBalanced                          = false;
  static double                                nextMaxPlotTimeStamp                  = FirstPlotTimeStamp;
  static double                                nextMinPlotTimeStamp                  = FirstPlotTimeStamp;
  static bool                                  haveJustWrittenSnapshot               = false;
  static bool                                  haveReceivedNoncriticialAssertion     = false;
  static bool                                  addGridSweepWithoutGridRefinementNext = false;
  static tarch::la::Vector<Dimensions, double> minH = tarch::la::Vector<Dimensions, double>(
    std::numeric_limits<double>::max()
  );
  static int globalNumberOfTrees = 0;
  bool       continueToSolve     = true;

  if (tarch::hasNonCriticalAssertionBeenViolated() and not haveReceivedNoncriticialAssertion) {
    peano4::parallel::Node::getInstance().setNextProgramStep(
      repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::PlotSolution)
    );
    haveReceivedNoncriticialAssertion = true;
    logError(
      "selectNextAlgorithmicStep()", "non-critical assertion has been triggered in code. Dump final state and terminate"
    );
  } else if (tarch::hasNonCriticalAssertionBeenViolated()) {
    continueToSolve = false;
  } else if (gridConstructed and not gridBalanced) {
    if (not repositories::loadBalancer.isEnabled(true) and not repositories::loadBalancer.hasSplitRecently()) {
      logInfo("selectNextAlgorithmicStep()", "all ranks have switched off their load balancing");
      gridBalanced = true;
    } else {
      logInfo(
        "selectNextAlgorithmicStep()", "wait for load balancing to become stable: " << repositories::loadBalancer
      );
    }

    peano4::parallel::Node::getInstance().setNextProgramStep(repositories::StepRepository::toProgramStep(
      repositories::StepRepository::Steps::CreateGridAndConvergeLoadBalancing
    ));
  } else if (gridBalanced and not gridInitialised) {
    peano4::parallel::Node::getInstance().setNextProgramStep(
      repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::InitGrid)
    );

    gridInitialised = true;
  } else if (not gridConstructed) {
    if (tarch::la::max(peano4::parallel::SpacetreeSet::getInstance().getGridStatistics().getMinH()) < tarch::la::max(minH)) {
      minH = peano4::parallel::SpacetreeSet::getInstance().getGridStatistics().getMinH();
      logDebug(
        "selectNextAlgorithmicStep()", "mesh has refined, so reset minH=" << minH << " and postpone further refinement"
      );
      addGridSweepWithoutGridRefinementNext = true;
    } else if (repositories::loadBalancer.getGlobalNumberOfTrees() > globalNumberOfTrees) {
      logInfo("selectNextAlgorithmicStep()", "mesh has rebalanced recently, so postpone further refinement)");
      addGridSweepWithoutGridRefinementNext = true;
      globalNumberOfTrees                   = repositories::loadBalancer.getGlobalNumberOfTrees();
    }
    else if (
      peano4::parallel::SpacetreeSet::getInstance().getGridStatistics().getStationarySweeps()>5
      and
      // ensure that a proper creation has been ran before, so the mesh had the opportunity
      // to refine further if it has not done yet
      repositories::StepRepository::toStepEnum( peano4::parallel::Node::getInstance().getCurrentProgramStep() ) == repositories::StepRepository::Steps::CreateGrid
    ) {
      logInfo(
        "selectNextAlgorithmicStep()", "grid has been stationary for quite some time. Terminate grid construction"
      );
      addGridSweepWithoutGridRefinementNext = false;
      gridConstructed                       = true;
    } else {
      logInfo(
        "selectNextAlgorithmicStep()",
        "mesh rebalancing seems to be stationary, so study whether to refine mesh further in next sweep: "
          << peano4::parallel::SpacetreeSet::getInstance().getGridStatistics().toString()
      );
      addGridSweepWithoutGridRefinementNext = false;
      globalNumberOfTrees                   = repositories::loadBalancer.getGlobalNumberOfTrees();
    }

    // Actual grid traversal choice
    if (addGridSweepWithoutGridRefinementNext) {
      peano4::parallel::Node::getInstance().setNextProgramStep(repositories::StepRepository::toProgramStep(
        repositories::StepRepository::Steps::CreateGridButPostponeRefinement
      ));
    } else {
      peano4::parallel::Node::getInstance().setNextProgramStep(
        repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::CreateGrid)
      );
    }

    continueToSolve = true;
  } else {
    if (TimeInBetweenPlots > 0.0 and repositories::getMinTimeStamp() < MinTerminalTime and repositories::getMaxTimeStamp() < MaxTerminalTime and (repositories::getMinTimeStamp() >= nextMinPlotTimeStamp or repositories::getMaxTimeStamp() >= nextMaxPlotTimeStamp) and repositories::mayPlot()) {
      if (repositories::getMinTimeStamp() >= nextMinPlotTimeStamp) {
        nextMinPlotTimeStamp += TimeInBetweenPlots;
      }
      if (repositories::getMaxTimeStamp() >= nextMaxPlotTimeStamp) {
        nextMaxPlotTimeStamp += TimeInBetweenPlots;
      }

      if (nextMinPlotTimeStamp < repositories::getMinTimeStamp()) {
        logWarning(
          "selectNextAlgorithmicStep()",
          "code is asked to plot every dt="
            << TimeInBetweenPlots << ", but this seems to be less than the time step size of the solvers. "
            << "So postpone next plot to t=" << (repositories::getMinTimeStamp() + TimeInBetweenPlots)
        );
        nextMinPlotTimeStamp = repositories::getMinTimeStamp() + TimeInBetweenPlots;
      } else if (nextMaxPlotTimeStamp < repositories::getMaxTimeStamp()) {
        logWarning(
          "selectNextAlgorithmicStep()",
          "code is asked to plot every dt="
            << TimeInBetweenPlots << ", but this seems to be less than the time step size of the solvers. "
            << "So postpone next plot to t=" << (repositories::getMaxTimeStamp() + TimeInBetweenPlots)
        );
        nextMaxPlotTimeStamp = repositories::getMaxTimeStamp() + TimeInBetweenPlots;
      }

      nextMaxPlotTimeStamp = std::max(nextMaxPlotTimeStamp, nextMinPlotTimeStamp);

      peano4::parallel::Node::getInstance().setNextProgramStep(
        repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::PlotSolution)
      );
      haveJustWrittenSnapshot = true;
      continueToSolve         = true;
    } else if (repositories::getMinTimeStamp() < MinTerminalTime and repositories::getMaxTimeStamp() < MaxTerminalTime) {
      peano4::parallel::Node::getInstance().setNextProgramStep(
        repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::TimeStep)
      );
      continueToSolve         = true;
      haveJustWrittenSnapshot = false;
    } else {
      if (not haveJustWrittenSnapshot and TimeInBetweenPlots > 0.0 and repositories::mayPlot()) {
        peano4::parallel::Node::getInstance().setNextProgramStep(
          repositories::StepRepository::toProgramStep(repositories::StepRepository::Steps::PlotSolution)
        );
        continueToSolve         = true; // don't want to terminate immediately
        haveJustWrittenSnapshot = true;
        nextMinPlotTimeStamp    = std::numeric_limits<double>::max();
        nextMaxPlotTimeStamp    = std::numeric_limits<double>::max();
      } else if (not haveJustWrittenSnapshot and TimeInBetweenPlots > 0.0 and not repositories::mayPlot()) {
        continueToSolve = true; // don't want to terminate immediately but to wait for incomplete time steps to complete
      } else {
        continueToSolve = false;
      }
    }
  }

  return continueToSolve;
}


void step() {
  int  stepIdentifier = peano4::parallel::Node::getInstance().getCurrentProgramStep();
  auto stepName       = repositories::StepRepository::toStepEnum(stepIdentifier);

  static tarch::logging::Log _log("");
#if PeanoDebug > 0
#else
  if (tarch::mpi::Rank::getInstance().isGlobalMaster())
#endif
  logInfo("step()", "run " << repositories::StepRepository::toString(stepName));

  static tarch::timing::Watch watch("::", "step()", false);

  static int creepingNumberOfLocalCells = 0;

  switch (stepName) {
  case repositories::StepRepository::Steps::CreateGridButPostponeRefinement: {
    tarch::logging::LogFilter::getInstance().switchProgramPhase("create-grid-but-postpone-refinement");

    OTTER_PHASE_SWITCH(otter::phase::create_grid_no_refine);

    repositories::startGridConstructionStep();

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::CreateGridButPostponeRefinement observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    gridConstructionMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishGridConstructionStep();
  } break;
  case repositories::StepRepository::Steps::CreateGrid: {
    tarch::logging::LogFilter::getInstance().switchProgramPhase("create-grid");

    OTTER_PHASE_SWITCH(otter::phase::create_grid);

    repositories::startGridConstructionStep();

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::CreateGrid observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    gridConstructionMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishGridConstructionStep();

    // We always overestimate so give the convergence the opportunity to catch up. The constant
    // here is a magic one.
    creepingNumberOfLocalCells = ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree()
                                 + tarch::multicore::Core::getInstance().getNumberOfThreads() * 3;
  } break;
  case repositories::StepRepository::Steps::CreateGridAndConvergeLoadBalancing: {
    if (creepingNumberOfLocalCells < ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree() - 1) {
      logInfo(
        "step()",
        "it seems the grid has just refined before we switched to the phase where we make the load balancing converge. Wait for a few iterations more to give load balancing chance to catch up"
      );
      creepingNumberOfLocalCells = ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree()
                                   + tarch::multicore::Core::getInstance().getNumberOfThreads() * 3;
    }

    tarch::logging::LogFilter::getInstance().switchProgramPhase("create-grid-and-converge-load-balancing");

    OTTER_PHASE_SWITCH(otter::phase::create_grid_converge);

    // The smaller here corresponds to the -1 below
    if (::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree() < 0 and repositories::loadBalancer.isEnabled(false)) {
      logInfo("step()", "rank is degenerated so disable load balancing temporarily");
      repositories::loadBalancer.enable(false);
    }
    if (
          ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree() >= creepingNumberOfLocalCells
          and
          repositories::loadBalancer.isEnabled(false)
        ) {
      logInfo(
        "step()",
        "grid construction and decomposition on this rank seem to be stable as we have around "
          << creepingNumberOfLocalCells << " local cells in the heaviest tree. Disable load balancing temporarily"
      );
      repositories::loadBalancer.enable(false);
    }

    repositories::startGridConstructionStep();

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::CreateGridButPostponeRefinement observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    gridConstructionMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishGridConstructionStep();

    if (
          ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree() <= creepingNumberOfLocalCells
          and
          not repositories::loadBalancer.hasSplitRecently()
          and
          repositories::loadBalancer.isEnabled(false)
        ) {
      logInfo(
        "step()",
        "have to decrement local cell counter "
          << creepingNumberOfLocalCells << " as maximum weight is "
          << ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree()
      );
      creepingNumberOfLocalCells = (creepingNumberOfLocalCells
                                    + ::toolbox::loadbalancing::getWeightOfHeaviestLocalSpacetree())
                                   / 2;
    }
  } break;
  case repositories::StepRepository::Steps::InitGrid: {
    tarch::logging::LogFilter::getInstance().switchProgramPhase("init-grid");
    repositories::loadBalancer.enable(false);

    OTTER_PHASE_SWITCH(otter::phase::init);

    repositories::startGridInitialisationStep();

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::InitGrid observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    gridConstructionMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishGridInitialisationStep();
  } break;
  case repositories::StepRepository::Steps::PlotSolution: {
    tarch::logging::LogFilter::getInstance().switchProgramPhase("plot-solution");
    const double minTimeStamp    = repositories::getMinTimeStamp();
    const double maxTimeStamp    = repositories::getMaxTimeStamp();
    const double minTimeStepSize = repositories::getMinTimeStepSize();
    const double maxTimeStepSize = repositories::getMaxTimeStepSize();

    OTTER_PHASE_SWITCH(otter::phase::plot);

    repositories::startPlottingStep(minTimeStamp, maxTimeStamp, minTimeStepSize, maxTimeStepSize);

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::PlotSolution observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    plotMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishPlottingStep();
  } break;
  case repositories::StepRepository::Steps::TimeStep: {
    tarch::logging::LogFilter::getInstance().switchProgramPhase("time-step");
    if (repositories::loadBalancer.isEnabled(false)) {
      logInfo("step()", "disable load balancing throughout initialisation (to be removed in later releases)");
      repositories::loadBalancer.enable(false);
    }

    const double minTimeStamp    = repositories::getMinTimeStamp();
    const double maxTimeStamp    = repositories::getMaxTimeStamp();
    const double minTimeStepSize = repositories::getMinTimeStepSize();
    const double maxTimeStepSize = repositories::getMaxTimeStepSize();
    const double minMeshSize     = repositories::getMinMeshSize();
    const double maxMeshSize     = repositories::getMaxMeshSize();

    OTTER_PHASE_SWITCH(otter::phase::timestep);

    repositories::startTimeStep(minTimeStamp, maxTimeStamp, minTimeStepSize, maxTimeStepSize);

    OTTER_DEFINE_TASK(step_task, OTTER_NULL_TASK, otter_no_add_to_pool, otter::label::step);
    OTTER_TASK_START(step_task);

    observers::TimeStep observer;
    watch.start();
    peano4::parallel::SpacetreeSet::getInstance().traverse(observer);
    watch.stop();
    timeStepMeasurement.setValue(watch.getCalendarTime());

    OTTER_TASK_WAIT_IMPLICIT(children);
    OTTER_TASK_END(step_task);
    OTTER_TASK_WAIT_IMPLICIT(descendants);

    repositories::finishTimeStep();
  } break;
  case repositories::StepRepository::Steps::Undef:
    assertion(false);
    break;
  }
}

int main(int argc, char** argv) {
  constexpr int ExitCodeSuccess          = 0;
  constexpr int ExitCodeUnitTestsFailed  = 1;
  constexpr int ExitCodeInvalidArguments = 2;
  constexpr int ExitCodeInvalidBuild     = 3;

  static tarch::timing::Watch watch("::", "main()", false);

  peano4::initParallelEnvironment(&argc, &argv);

  // Do this early, so people can use logInfo properly.
  repositories::initLogFilters();

  tarch::initNonCriticalAssertionEnvironment();
  tarch::multicore::initSmartMPI();
  peano4::fillLookupTables();

  peano4::initSingletons(DomainOffset, DomainSize, PeriodicBC);

  repositories::initSharedMemoryAndGPUEnvironment();

  if (tarch::mpi::Rank::getInstance().getNumberOfRanks() > 1 and tarch::multicore::Core::getInstance().getNumberOfThreads() <= 1) {
    logError("main()", "MPI runs without multithreading are not supported currently.");
    return ExitCodeInvalidBuild;
  }

  repositories::DataRepository::initDatatypes();

  {% if FENV_ARGS is defined and FENV_ARGS != "" -%}
  feenableexcept({{FENV_ARGS}} );
  {% endif -%}

#if PeanoDebug >= 2
  tarch::tests::TreeTestCaseCollection* unitTests = new tarch::tests::TreeTestCaseCollection();
  unitTests->addTestCase(peano4::getUnitTests());
  unitTests->addTestCase(tarch::getUnitTests());
  unitTests->addTestCase(toolbox::blockstructured::getUnitTests());
  unitTests->addTestCase(exahype2::getUnitTests());
  unitTests->run();
  if (unitTests->getNumberOfErrors() != 0) {
    logError("main()", "unit tests failed. Quit.");
    tarch::mpi::Rank::abort(ExitCodeUnitTestsFailed);
  }
  delete unitTests;
#endif

  repositories::startSimulation();

  tarch::logging::Statistics::getInstance().clear();

  OTTER_INITIALISE();

#if defined(SharedOMP)
#pragma omp parallel
  {
#pragma omp master
    {
#endif

#if defined(UseSmartMPI)
      const bool isGlobalMaster = tarch::mpi::Rank::getInstance().isGlobalMaster() and smartmpi::isComputeRank();
      const bool isPeanoComputeNode = not tarch::mpi::Rank::getInstance().isGlobalMaster() and smartmpi::isComputeRank();
#else
      const bool isGlobalMaster     = tarch::mpi::Rank::getInstance().isGlobalMaster();
      const bool isPeanoComputeNode = not tarch::mpi::Rank::getInstance().isGlobalMaster();
#endif

      if (isGlobalMaster) {
        while (selectNextAlgorithmicStep()) {
          watch.start();
          step();
          watch.stop();

          timePerMeshSweepMeasurement.setValue(watch.getCalendarTime());
          logInfo(
            "main()",
            "time per mesh sweep (current/average): " << std::fixed << std::setprecision(2) << watch.getCalendarTime() <<
            "s / " << timePerMeshSweepMeasurement.getValue() << "s"
          );
        }

        logInfo("main()", "terminated successfully");
        logInfo(
          "main()",
          "initial grid construction:   " << gridConstructionMeasurement.getAccumulatedValue() <<
          "s\t" << gridConstructionMeasurement.toString()
        );
        logInfo(
          "main()",
          "plotting:                    " << plotMeasurement.getAccumulatedValue() <<
          "s\t" << plotMeasurement.toString()
        );
        logInfo(
          "main()",
          "time stepping:               " << timeStepMeasurement.getAccumulatedValue() <<
          "s\t" << timeStepMeasurement.toString()
        );
        logInfo(
          "main()",
          "average time per mesh sweep: " << timePerMeshSweepMeasurement.getValue() <<
          "s\t" << timePerMeshSweepMeasurement.toString()
        );
      } else if (isPeanoComputeNode) {
        while (peano4::parallel::Node::getInstance().continueToRun()) {
          step();
        }
      }
#if defined(UseSmartMPI)
      else {
        while (smartmpi::continueToRun()) {
          smartmpi::tick();
        }
      }
#endif

#if defined(SharedOMP)
    }
  }
#endif

  tarch::logging::Statistics::getInstance().writeToCSV();

  OTTER_FINALISE();

  repositories::finishSimulation();

  peano4::shutdownSingletons();
  repositories::DataRepository::shutdownDatatypes();
  tarch::multicore::shutdownSmartMPI();
  tarch::shutdownNonCriticalAssertionEnvironment();
  peano4::shutdownParallelEnvironment();

  return ExitCodeSuccess;
}
