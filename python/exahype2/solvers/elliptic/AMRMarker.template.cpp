// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#include "{{CLASSNAME}}.h"
#include "exahype2/RefinementControl.h"


tarch::logging::Log   {{NAMESPACE | join("::")}}::{{CLASSNAME}}::_log( "{{NAMESPACE | join("::")}}::{{CLASSNAME}}" );


std::bitset<Dimensions> {{NAMESPACE | join("::")}}::{{CLASSNAME}}::{{CLASSNAME}}::PeriodicBC = {{NAMESPACE | join("::")}}::PeriodicBC;


{{NAMESPACE | join("::")}}::{{CLASSNAME}}::{{CLASSNAME}}() {
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMinTimeStamp(bool ofCurrentlyRunningGridSweep) const {
  return std::numeric_limits<double>::max();
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMaxTimeStamp(bool ofCurrentlyRunningGridSweep) const {
  return 0.0;
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMinTimeStepSize() const {
  return 0.0;
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMaxTimeStepSize() const {
  return std::numeric_limits<double>::max();
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMaxMeshSize() const {
  return std::numeric_limits<double>::max();
}


double {{NAMESPACE | join("::")}}::{{CLASSNAME}}::getMinMeshSize() const {
  return 0.0;
}


bool {{NAMESPACE | join("::")}}::{{CLASSNAME}}::mayPlot() const {
  return true;
}


bool {{NAMESPACE | join("::")}}::{{CLASSNAME}}::isFirstGridSweepOfTimeStep() const {
  return true;
}


bool {{NAMESPACE | join("::")}}::{{CLASSNAME}}::isLastGridSweepOfTimeStep() const {
  return true;
}


void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::startGridConstructionStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::finishGridConstructionStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::startGridInitialisationStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::finishGridInitialisationStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::startTimeStep(
  double globalMinTimeStamp,
  double globalMaxTimeStamp,
  double globalMinTimeStepSize,
  double globalMaxTimeStepSize
) {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::finishTimeStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::startPlottingStep(
  double globalMinTimeStamp,
  double globalMaxTimeStamp,
  double globalMinTimeStepSize,
  double globalMaxTimeStepSize
) {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::finishPlottingStep() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::startSimulation() {}
void {{NAMESPACE | join("::")}}::{{CLASSNAME}}::finishSimulation() {}
