//
// Solver file of Peano's PETSc add-on
// Generated by Peano's Python API
// www.peano-framework.org
//
// This is generated. Don't change it! Every rerun of the Python API will
// overwrite your changes.
//
#pragma once


#include "Abstract{{CLASSNAME}}.h"


{% for item in NAMESPACE -%}
  namespace {{ item }} {

{%- endfor %}
  class {{CLASSNAME}};

{% for item in NAMESPACE -%}
  }
{%- endfor %}


class {{NAMESPACE | join("::")}}::{{CLASSNAME}}: public {{NAMESPACE | join("::")}}::Abstract{{CLASSNAME}} {
  public:

    /**
     * Default constructor
     *
     * @todo Please add your documentation here.
     */
    {{CLASSNAME}}();

    /**
     * Destructor
     *
     * Has to be virtual, as there is a superclass with virtual functions.
     *
     * @todo Please add your documentation here.
     */
    virtual ~{{CLASSNAME}}();

    /**
     * Initialise a degree of freedom
     *
     * This routine will be called only for interior DoFs. See the correlation
     * to getCellDoFType() in the superclass.
     *
     * @todo Please add your documentation here
     */
    virtual void initNode(
      const tarch::la::Vector<Dimensions, double>&  x,
      const tarch::la::Vector<Dimensions, double>&  h,
      {% if CELL_UNKNOWNS==1 %}
      double&                                       value,
      double&                                       rhs,
      double&                                       exactSol
      {% else %}
      tarch::la::Vector< {{CELL_UNKNOWNS}}, double >&  value,
      tarch::la::Vector< {{CELL_UNKNOWNS}}, double >&  rhs,
      tarch::la::Vector< {{CELL_UNKNOWNS}}, double >&  exactSol
      {% endif %}
    ) override;


    {% if CELL_CELL_RHS_MATRIX==[] %}
    virtual tarch::la::Matrix< DoFsPerCell, DoFsPerCell, double > getLhsMatrix(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    ) override;
    {% endif %}


    {% if CELL_CELL_RHS_MATRIX==[] %}
    virtual tarch::la::Matrix< DoFsPerCell, DoFsPerCell, double > getRhsMatrix(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    ) override;
    {% endif %}


    {% if CELL_TO_FACE_MATRIX==[] %}
    virtual tarch::la::Matrix< NodesPerFace*FaceUnknowns*FacesPerCell, DoFsPerCell, double > getProjectionOfCellDataOntoFace(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    ) override;
    {% endif %}


    {% if FACE_TO_CELL_MATRIX==[] %}
    virtual tarch::la::Matrix< DoFsPerCell, NodesPerFace*FaceUnknowns*FacesPerCell, double > getProjectionOfRiemannSolutionOntoCell(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    ) override;
    {% endif %}


    {% if FACE_FACE_RIEMANN_PROBLEM_MATRIX==[] %}
    virtual tarch::la::Matrix< NodesPerFace*FaceUnknowns, 2*NodesPerFace*FaceUnknowns, double > getRiemannSolver(
      const tarch::la::Vector<Dimensions, double>&  faceCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    ) override;
    {% endif %}
};


