  // **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#include "{{CLASSNAME}}.h"

tarch::logging::Log   {{NAMESPACE | join("::")}}::{{CLASSNAME}}::_log( "{{NAMESPACE | join("::")}}::{{CLASSNAME}}" );


{% if ADMISSIBILITY_IMPLEMENTATION=="<user-defined>" %}
bool {{NAMESPACE | join("::")}}::{{CLASSNAME}}::isPhysicallyAdmissible(
  const double* const                         Q,
  const tarch::la::Vector<Dimensions,double>& x,
  const tarch::la::Vector<Dimensions,double>& h,
  const double timeStamp
){
  return true;
}
{% endif %}