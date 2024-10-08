# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from peano4.solversteps.ActionSet import ActionSet



class ProjectPatchOntoFaces(ActionSet):
  """!
  
  This class assumes that you have NxNxN patch within your block. It also 
  assumes that you have an 2MxNxN patch on your faces. M<N. The code interprets
  these face-associated data as overlap with the patch, hooks into 
  touchCellLastTime, and projects the cell's patch data onto the overlap 
  (simple copy). 
  
  And then runs over your grid and maps the patch data on the respective entries of the 
  face patch which is an auxiliary data structure.
  
  
  # Attributes
  
  patch: peano4.datamodel.Patch
    This is the input data structure, i.e. the patch from which we extract
    the boundary data.
    
  patch_overlap: peano4.datamodel.Patch. 
    Consult remark above about how the dimensions of this overlap patch have 
    to match. 
    
  """
  
  
  def __init__(self,patch,patch_overlap,guard,additional_includes, add_assertions = True):
    super(ProjectPatchOntoFaces,self).__init__(descend_invocation_order=1,parallel=False)
    self.d = {}
    if patch_overlap.dim[0] % 2 != 0:
      print( "Error: Patch associated to face has to have even number of cells. Otherwise, it is not a symmetric overlap." )
      assert( patch_overlap.dim[0] % 2 == 0 )
    if patch.dim[0] != patch.dim[1]:
      print( "Error: Can only handle square patches." )
      assert( patch.dim[0] == patch.dim[1] )
    if patch_overlap.dim[1] != patch.dim[0]:
      print( "Error: Patch of overlap and patch of cell have to match" )
      assert( patch_overlap.dim[1] == patch.dim[0] )
      
    assert isinstance(patch.no_of_unknowns,int)
    assert isinstance(patch.dim[0],int)
    assert isinstance(patch_overlap.dim[0],int)
      
    self.d[ "UNKNOWNS" ]           = str(patch.no_of_unknowns)
    self.d[ "DOFS_PER_AXIS" ]      = str(patch.dim[0])
    self.d[ "OVERLAP" ]            = str(int(patch_overlap.dim[0]/2))
    self.d[ "FACES_ACCESSOR" ]     = "fineGridFaces"  + patch_overlap.name
    self.d[ "CELL_ACCESSOR" ]      = "fineGridCell" + patch.name
    self.d[ "GUARD" ]              = guard
    if add_assertions:
      self.d[ "ASSERTION_PREFIX" ] = "false"
    else:
      self.d[ "ASSERTION_PREFIX" ] = "true"
      
    
    self.additional_includes = additional_includes


  @property
  def guard(self):
      return self.d[ "GUARD" ]
  
  @guard.setter
  def guard(self,value):
    self.d[ "GUARD" ]              = value

  def get_constructor_body(self):
    return ""

  
  def get_destructor_body(self):
    return ""


  def get_body_of_getGridControlEvents(self):
    return "  return std::vector< peano4::grid::GridControlEvent >();\n" 


  def get_action_set_name(self):
    return __name__.replace(".py", "").replace(".", "_")


  def user_should_modify_template(self):
    return False


  __Template_TouchCellLastTime_Preamble = """
  if ( {GUARD} ) {{
    logTraceIn( "touchCellLastTime(...)---ProjectPatchOntoFaces" );
"""


  __Template_TouchCellLastTime_Core = """
 // @todo Might want to use routine from Projection.h in toolbox.
    for(int d=0; d<Dimensions; d++) {{
      /**
       * d-loop over all dimensions except d. The vector k's entry d is set
       * to 0. We start with the left/bottom face, i.e. the one closer to the 
       * coordinate system's origin.
       */
      dfore(k,{DOFS_PER_AXIS},d,0) {{
        for (int i=0; i<{OVERLAP}; i++) {{
          tarch::la::Vector<Dimensions,int> patchCell   = k;
          tarch::la::Vector<Dimensions,int> overlapCell = k;
          patchCell(d)   = i;
          overlapCell(d) = i+{OVERLAP};
          
          int patchCellSerialised   = peano4::utils::dLinearised(patchCell,{DOFS_PER_AXIS});
          int overlapCellSerialised = toolbox::blockstructured::serialiseVoxelIndexInOverlap(overlapCell,{DOFS_PER_AXIS},{OVERLAP},d);
          for (int j=0; j<{UNKNOWNS}; j++) {{
            assertion7( 
              {ASSERTION_PREFIX}
              or
              {CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j]=={CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j], j,i,k,d, patchCell, overlapCell,
              marker.toString() 
            );
            {FACES_ACCESSOR}(d).value[overlapCellSerialised*{UNKNOWNS}+j] = 
              {CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j];
          }}
  
          patchCell(d)   = i+{DOFS_PER_AXIS}-{OVERLAP};
          overlapCell(d) = i;
          
          patchCellSerialised   = peano4::utils::dLinearised(patchCell,{DOFS_PER_AXIS});
          overlapCellSerialised = toolbox::blockstructured::serialiseVoxelIndexInOverlap(overlapCell,{DOFS_PER_AXIS},{OVERLAP},d);
          for (int j=0; j<{UNKNOWNS}; j++) {{
            assertion7( 
              {ASSERTION_PREFIX}
              or
              {CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j]=={CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j], j,i,k,d, patchCell, overlapCell, 
              marker.toString() 
            );
            {FACES_ACCESSOR}(d+Dimensions).value[overlapCellSerialised*{UNKNOWNS}+j] = 
              {CELL_ACCESSOR}.value[patchCellSerialised*{UNKNOWNS}+j];
          }}
        }}
      }}
    }}
"""


  __Template_TouchCellLastTime_Epilogue = """
    logTraceOut( "touchCellLastTime(...)---ProjectPatchOntoFaces" );
  }}
  else {{
    logTraceInWith1Argument( "touchCellLastTime(...)---ProjectPatchOntoFaces [skip]", marker.toString() );
    logTraceOut( "touchCellLastTime(...)---ProjectPatchOntoFaces [skip]" );
  }}
"""


  def get_body_of_operation(self,operation_name):
    result = "\n"
    if operation_name==ActionSet.OPERATION_TOUCH_CELL_LAST_TIME:
      result  = self.__Template_TouchCellLastTime_Preamble.format(**self.d)
      result += self.__Template_TouchCellLastTime_Core.format(**self.d)
      result += self.__Template_TouchCellLastTime_Epilogue.format(**self.d)
      pass 
    return result


  def get_attributes(self):
    return ""


  def get_includes(self):
    return """
#include "peano4/utils/Loop.h"
#include "toolbox/blockstructured/Enumeration.h"
""" + self.additional_includes
