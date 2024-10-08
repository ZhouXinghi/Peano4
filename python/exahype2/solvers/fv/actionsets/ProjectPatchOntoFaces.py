# This file is part of the ExaHyPE2 project. For conditions of distribution and 
# use, please see the copyright notice at www.peano-framework.org
import peano4
import jinja2


import peano4.toolbox.blockstructured.ProjectPatchOntoFaces


class ProjectPatchOntoFaces( peano4.toolbox.blockstructured.ProjectPatchOntoFaces ):
  """! Project patch data onto faces so we can reconstruct the halo layers later on 
  
   This class is very similar to peano4.toolbox.blockstructured.ProjectPatchOntoFaces, 
   but there are some delicate differences.
   
   The big difference is that the projection also sets the isUpdated() flag and 
   the UpdatedTimeStamp(). As the projection writes to these two updated 
   records, it is important that you roll it over afterwards. This is done via 
   the mapping RollOverUpdatedFace.
   
   It is important to study this action set in combination with DynamicAMR. In 
   the documentation of the latter I explain why we need the guard 
   
   not marker.hasBeenRefined()
   
   if we want to support dynamic coarsening. See exahype2.solvers.fv.FV.add_actions_to_init_grid() 
   for details on the ordering of the action sets.
      
  """
  def __init__(self,solver, guard):
    peano4.toolbox.blockstructured.ProjectPatchOntoFaces.__init__(
      self,
      solver._patch,
      solver._patch_overlap_update,
      "not marker.hasBeenRefined() and " + guard,
      solver._get_default_includes() + solver.user_action_set_includes, True)

    self.d[ "FACE_METADATA_ACCESSOR" ] = "fineGridFaces"  + solver._face_label.name
    self.d[ "CELL_METADATA_ACCESSOR" ] = "fineGridCell""" + solver._cell_label.name

    self.__Template_TouchCellLastTime_Core += """
  for (int d=0; d<Dimensions; d++) {{
    {FACE_METADATA_ACCESSOR}(d).setUpdated(1,true);
    {FACE_METADATA_ACCESSOR}(d).setUpdatedTimeStamp(1,{CELL_METADATA_ACCESSOR}.getTimeStamp());
    {FACE_METADATA_ACCESSOR}(d+Dimensions).setUpdated(0,true);
    {FACE_METADATA_ACCESSOR}(d+Dimensions).setUpdatedTimeStamp(0,{CELL_METADATA_ACCESSOR}.getTimeStamp());
    logDebug( "touchCellLastTime(...)", "update {FACE_METADATA_ACCESSOR}(" << d << ")(1)" );
    logDebug( "touchCellLastTime(...)", "update {FACE_METADATA_ACCESSOR}(" << (d+Dimensions) << ")(0)" );
  }}  
"""


  def get_action_set_name(self):
    return __name__.replace(".py", "").replace(".", "_")

