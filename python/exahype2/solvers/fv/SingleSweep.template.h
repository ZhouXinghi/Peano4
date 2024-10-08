// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#pragma once

#include "Abstract{{CLASSNAME}}.h"
#include "tarch/logging/Log.h"

{% for item in NAMESPACE -%}
  namespace {{ item }} {

{%- endfor %}
  class {{CLASSNAME}};

{% for item in NAMESPACE -%}
  }
{%- endfor %}

class {{NAMESPACE | join("::")}}::{{CLASSNAME}}: public {{NAMESPACE | join("::")}}::Abstract{{CLASSNAME}} {
  private:
    static tarch::logging::Log   _log;

  public:
    {% if REFINEMENT_CRITERION_IMPLEMENTATION=="<user-defined>" %}
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
     *  to a valid value throughout grid construction. In this case,
     *  it is nullptr.
     */
    virtual ::exahype2::RefinementCommand refinementCriterion(
      [[maybe_unused]] const double* __restrict__                   Q, // Q[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& volumeCentre,
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& volumeH,
      [[maybe_unused]] double                                       t
    ) override;
    {% endif %}

    {% if INITIAL_CONDITIONS_IMPLEMENTATION=="<user-defined>" %}
    virtual void initialCondition(
      [[maybe_unused]] double* __restrict__                         Q, // Q[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& volumeCentre,
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& volumeH,
      [[maybe_unused]] bool                                         gridIsConstructed
    ) override;
    {% endif %}

    {% if BOUNDARY_CONDITIONS_IMPLEMENTATION=="<user-defined>" %}
    virtual void boundaryConditions(
      [[maybe_unused]] const double* __restrict__                   Qinside, // Qinside[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      [[maybe_unused]] double* __restrict__                         Qoutside, // Qoutside[{{NUMBER_OF_UNKNOWNS}}+{{NUMBER_OF_AUXILIARY_VARIABLES}}]
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& faceCentre,
      [[maybe_unused]] const tarch::la::Vector<Dimensions, double>& volumeH,
      [[maybe_unused]] double                                       t,
      [[maybe_unused]] int                                          normal
    ) override;
    {% endif %}

    {{SOLVER_USER_DECLARATIONS}}
};
