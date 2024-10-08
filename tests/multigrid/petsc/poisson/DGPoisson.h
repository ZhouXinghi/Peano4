//
// Solver file of Peano's PETSc add-on
// Generated by Peano's Python API
// www.peano-framework.org
//
// This is generated. Don't change it! Every rerun of the Python API will
// overwrite your changes.
//
#pragma once


#include "AbstractDGPoisson.h"


namespace benchmarks {namespace multigrid {namespace petsc {namespace poisson {
  class DGPoisson;

}}}}


class benchmarks::multigrid::petsc::poisson::DGPoisson: public benchmarks::multigrid::petsc::poisson::AbstractDGPoisson {
  public:

    /**
     * Default constructor
     *
     * @todo Please add your documentation here.
     */
    DGPoisson();

    /**
     * Destructor
     *
     * Has to be virtual, as there is a superclass with virtual functions.
     *
     * @todo Please add your documentation here.
     */
    virtual ~DGPoisson();

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
      
      double&                                       value,
      double&                                       rhs,
      double&                                       exactSol
      
    ) override;
    

    /*
    Here we provide some custom scaling for this function. 

    Please ensure we passed "None" to CELL_TO_FACE_MATRIX_SCALING

    */
    tarch::la::Matrix< 
      benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::NodesPerFace*benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::FaceUnknowns*benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::FacesPerCell,
      benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::DoFsPerCell,
      double > getProjectionOfCellDataOntoFace(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
    );
    
    tarch::la::Matrix< benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::DoFsPerCell, 
      benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::NodesPerFace*benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::FaceUnknowns*benchmarks::multigrid::petsc::poisson::AbstractDGPoisson::FacesPerCell,
      double > getProjectionOfRiemannSolutionOntoCell(
      const tarch::la::Vector<Dimensions, double>&  cellCentre,
      const tarch::la::Vector<Dimensions, double>&  cellSize
      );

};

