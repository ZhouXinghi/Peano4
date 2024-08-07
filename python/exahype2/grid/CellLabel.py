# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from peano4.solversteps.ActionSet import ActionSet

import dastgen2
import peano4.datamodel.DaStGen2


class UpdateCellLabel(ActionSet):
  """! Update the cell label within a sweep
  
  Update the cell label initialises cell labels correctly when a cell is 
  created, and it also resets the update flag.
  
  It is important that this label update happens before the actual cell
  update. This is ensured via a correct invocation order in 
  FV.add_actions_to_perform_time_step(). 
  
  Cell labels exist everywhere and are not masked out via the LoadStoreCompute
  flag. That is, we store them persistently in each and every cell on all 
  resolution levels all the time - even if the underlying solver des not 
  exist in such a cell. In this case, the cell label should be false all the
  time.
  
  """
  def get_attribute_name(solver_name):
    return solver_name + "CellLabel"

  """
  
   SetLabels is an action set which is automatically merged into all ExaHyPE2
   steps by the project. There's nothing for the user solvers to be done here.
  
  """
  def __init__(self, solver_name):
    super(UpdateCellLabel,self).__init__(descend_invocation_order=1,parallel=False)
    self._solver_name = solver_name
    pass


  def get_constructor_body(self):
    return ""

    
  def get_destructor_body(self):
    return ""


  #def get_body_of_getGridControlEvents(self):
  #  return self.__Template_GridControlEvents.format(**self.d)


  def get_action_set_name(self):
    return __name__.replace(".py", "").replace(".", "_")


  def user_should_modify_template(self):
    return False


  def get_body_of_operation(self,operation_name):
    result = "\n"
    if operation_name==ActionSet.OPERATION_CREATE_CELL:
      result += """
  fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setSemaphoreNumber(  ::exahype2::EnclaveBookkeeping::NoEnclaveTaskNumber );
  fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setAMDCriterionEvaluatedThroughoutGridConstruction( false );
  if ( tarch::la::equals(repositories::getMinTimeStamp(),0.0) ) {
    fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setTimeStamp(  0.0 );
    fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setTimeStepSize( 0.0 );
  }
  else {
    fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setTimeStamp(  coarseGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.getTimeStamp() );
    fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setTimeStepSize( coarseGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.getTimeStepSize() );
  }
  fineGridCell""" + UpdateCellLabel.get_attribute_name(self._solver_name) + """.setHasUpdated( false );
"""
    return result


  def get_attributes(self):
    return ""


  def get_includes(self):
    return """
#include "Constants.h"
#include "exahype2/EnclaveBookkeeping.h"
#include "repositories/SolverRepository.h"
"""    




def create_cell_label(solver_name):
  """!

   Build up the DaStGen2 data structure that we need per cell to maintain
   per-cell data per solver.
   
   :: Attributes per cell
   
   - SemaphoreNumber is a number that links to a task that has been spawn in the
     background. Consult the enclave tasking and its bookkeeping for details.
   - TimeStamp holds the local time stamp of the current cell. Even for global 
     time stepping, it is important that this time stamp is stored per cell, as 
     some postprocessing (particles, e.g.) rely on the field.
   - TimeStepSize (similar to arguments around TimeStamp).
   - AMDCriterionEvaluatedThroughoutGridConstruction
     
   solver_name: string
     Name of the solver
     
  """
  result = peano4.datamodel.DaStGen2( UpdateCellLabel.get_attribute_name( solver_name ) )
  
  result.data.add_attribute( dastgen2.attributes.Integer("SemaphoreNumber") )
  result.data.add_attribute( dastgen2.attributes.Double("TimeStamp") )
  result.data.add_attribute( dastgen2.attributes.Double("TimeStepSize") )
  result.data.add_attribute( dastgen2.attributes.Boolean("HasUpdated") )
  result.data.add_attribute( dastgen2.attributes.Boolean("AMDCriterionEvaluatedThroughoutGridConstruction") )
    
  result.peano4_mpi_and_storage_aspect.merge_implementation = """
  switch (context) {
    case ::peano4::grid::TraversalObserver::SendReceiveContext::BoundaryExchange:
      assertionMsg( false, "cells should never be subject of a boundary exchange" );
      break;
    case ::peano4::grid::TraversalObserver::SendReceiveContext::ForkDomain:
      *this = neighbour;
      break;
    default:
      assertionMsg( false, "not implemented yet" );
  }
"""

  return result
