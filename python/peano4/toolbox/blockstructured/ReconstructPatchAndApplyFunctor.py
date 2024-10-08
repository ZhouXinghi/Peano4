# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from enum import IntEnum

from peano4.solversteps.ActionSet import ActionSet

class ReconstructedArrayMemoryLocation(IntEnum):
  """
   All arrays are held on the call stack. Might not work with PGI compilers or
   Clang, but it does work with GCC and Intel. This mode is the preferred one,
   as all memory frees are implicitly done.
  """
  CallStack = 0,
  """
   Create data on the heap. Done via a plain new. You have to delete the
   reconstructed data yourself.
  """
  Heap = 1,
  HeapWithoutDelete = 2,
  """
    Use tarch::allocateMemory()
  """  
  HeapThroughTarch = 3,
  HeapThroughTarchWithoutDelete = 4,
  """
    Use tarch::allocateMemory()
  """
  ManagedSharedAcceleratorDeviceMemoryThroughTarch = 5,
  ManagedSharedAcceleratorDeviceMemoryThroughTarchWithoutDelete = 6


class ReconstructPatchAndApplyFunctor(ActionSet):
  """!
  
  This class assumes that you have NxNxN patch within your block. It also assumes that
  you have an 2MxNxN patch on your faces. M<N. The code interprets these face-associated
  data as overlap with the patch, hooks into touchCellLastTime, and projects the cell's
  patch data onto the overlap (simple copy). See ProjectPatchOntoFaces for details.

  If the face patches are auxiliary patches holding copies from the cell patches, then
  we can reconstruct N+2M x N+2M x N+2M patches per cell whenever we hit a cell: We
  create a new temporary array, and copy data from both the cell patch and the faces
  into this array (gather operation). After that, we can launch the passed functor
  giving it access to the temporary, large array plus the original patch data.

  If you want to swap the functor, please replace self.functor_implementation, 
  and ensure that all external references (text replacements) are already
  resolved. Consult _add_action_set_entries_to_dictionary() for details.
  
  """

  def __init__(self,
               patch,
               patch_overlap,
               functor_implementation,
               reconstructed_array_memory_location=ReconstructedArrayMemoryLocation.Heap,
               guard="true", 
               add_assertions_to_halo_exchange = True
               ):
    """!
    
  patch: peano4.datamodel.Patch
    Patch which is to be used

  patch_overlap: peano4.datamodel.Patch
    Consult remark above about how the dimensions of this overlap 
    patch have to match

  ## Functor_implementation 

  The functor implementation is a plain C/C++ 

  - If you want to use brackets in your implementation, please use double brackets {{ }} as 
    the template system otherwise gets confused.
  - The following C++ variables are defined:

    oldQWithHalo
    newQ

    Both are plain double pointers.

  ## Data validation

  I use quite a lot of validitiy checks for the copied data via comparison to its self. This
  way, I can at least spot nans. I use the dictionary entries ASSERTION_WITH_X_ARGUMENTS to 
  realise this, and they are, by default, set to Peano's assertion macros. You can redefine 
  them.
  
    """
    super(ReconstructPatchAndApplyFunctor,self).__init__(descend_invocation_order=1,parallel=False)
    
    assert isinstance(reconstructed_array_memory_location,ReconstructedArrayMemoryLocation), "type of argument is {}".format( type(reconstructed_array_memory_location) )
    
    if patch_overlap.dim[0] % 2 != 0:
      print( "Error: Patch associated to face has to have even number of cells. Otherwise, it is not a symmetric overlap." )
      assert( patch_overlap.dim[0] % 2 == 0 )
    if patch.dim[0] != patch.dim[1]:
      print( "Error: Can only handle square patches." )
      assert( patch.dim[0] == patch.dim[1] )
    if patch_overlap.dim[1] != patch.dim[0]:
      print( "Error: Patch of overlap and patch of cell have to match" )
      assert( patch_overlap.dim[1] == patch.dim[0] )

    self.guard          = guard
    self.no_of_unknowns = int(patch.no_of_unknowns)
    self.dofs_per_axis  = int(patch.dim[0])
    self.total_overlap  = int(patch_overlap.dim[0])
    self.overlap_name   = patch_overlap.name
    self.patch_name     = patch.name
    
    self.functor_implementation              = functor_implementation
    self.add_assertions_to_halo_exchange     = add_assertions_to_halo_exchange
    self.reconstructed_array_memory_location = reconstructed_array_memory_location


  def _add_action_set_entries_to_dictionary(self,d):
    """!
    
    Befill dictionary
    
    If you inherit from this class and specialise this routine, you can add
    more entries to the dictionary. This is important if you change the 
    actual compute kernels (templates) and use futher unknowns in there. If 
    you specialise the function, i.e. overwrite it in a subclass, please
    ensure that you use the superclass still.
    
    Please note that we do not replace entries in the functor implementation:
    The overall code created by this action set consists of different parts
    realised as different templates. The templates are subject to a text 
    replacement system and then are concatenated. One of the things that is
    replaced within the templates is the actual functor implementation 
    (CELL_FUNCTOR_IMPLEMENTATION). If this one contains again free parameters,
    it is your responsibility to first resolve those parameters, and then to 
    reset the functor in this routine.
    
    If you want to replace the functor, you can hence either first call this
    routine and then change d["CELL_FUNCTOR_IMPLEMENTATION"] or you set
    self.functor_implementation and then you call this routine, i.e. the 
    superclass.
        
    """
    d[ "GUARD" ]              = self.guard
    
    assert isinstance(self.no_of_unknowns,int)
    assert isinstance(self.dofs_per_axis,int)
    assert isinstance(self.total_overlap,int)

    d[ "UNKNOWNS" ]           = str(self.no_of_unknowns)
    d[ "DOFS_PER_AXIS" ]      = str(self.dofs_per_axis)
    d[ "OVERLAP" ]            = str(int(self.total_overlap/2))
    d[ "NUMBER_OF_DOUBLE_VALUES_IN_ORIGINAL_PATCH_2D" ]      = str(self.no_of_unknowns * self.dofs_per_axis * self.dofs_per_axis)
    d[ "NUMBER_OF_DOUBLE_VALUES_IN_ORIGINAL_PATCH_3D" ]      = str(self.no_of_unknowns * self.dofs_per_axis * self.dofs_per_axis * self.dofs_per_axis)
    d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_2D" ] = str(self.no_of_unknowns * (self.total_overlap + self.dofs_per_axis) * (self.total_overlap + self.dofs_per_axis))
    d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_3D" ] = str(self.no_of_unknowns * (self.total_overlap + self.dofs_per_axis) * (self.total_overlap + self.dofs_per_axis) * (self.total_overlap + self.dofs_per_axis))
    d[ "FACES_ACCESSOR" ]     = "fineGridFaces"  + self.overlap_name
    d[ "CELL_ACCESSOR" ]      = "fineGridCell" + self.patch_name

    d[ "ASSERTION_WITH_1_ARGUMENTS" ] = "assertion1"
    d[ "ASSERTION_WITH_2_ARGUMENTS" ] = "assertion2"
    d[ "ASSERTION_WITH_3_ARGUMENTS" ] = "assertion3"
    d[ "ASSERTION_WITH_4_ARGUMENTS" ] = "assertion4"
    d[ "ASSERTION_WITH_5_ARGUMENTS" ] = "assertion5"
    d[ "ASSERTION_WITH_6_ARGUMENTS" ] = "assertion6"
    d[ "ASSERTION_WITH_7_ARGUMENTS" ] = "assertion7"

    d[ "CELL_FUNCTOR_IMPLEMENTATION" ] = self.functor_implementation

    if self.add_assertions_to_halo_exchange:
      d[ "ASSERTION_PREFIX_FOR_HALO" ] = "false"
    else:
      d[ "ASSERTION_PREFIX_FOR_HALO" ] = "true"

    if self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.Heap or self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.HeapWithoutDelete:
      d[ "CREATE_RECONSTRUCTED_PATCH" ] = """
    #if Dimensions==2
      double* oldQWithHalo = new double[""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_2D" ] + """];
    #elif Dimensions==3
      double* oldQWithHalo = new double[""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_3D" ] + """];
    #endif
"""
    elif self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.CallStack:
      d[ "CREATE_RECONSTRUCTED_PATCH" ] = """
    #if Dimensions==2
      double oldQWithHalo[""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_2D" ] + """];
    #elif Dimensions==3
      double oldQWithHalo[""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_3D" ] + """];
    #endif
"""
    elif self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.HeapThroughTarch or self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.HeapThroughTarchWithoutDelete:
      d[ "CREATE_RECONSTRUCTED_PATCH" ] = """
    #if Dimensions==2
      double* oldQWithHalo = ::tarch::allocateMemory<double>(""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_2D" ] + """, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory);
    #elif Dimensions==3
      double* oldQWithHalo = ::tarch::allocateMemory<double>(""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_3D" ] + """, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory);
    #endif
"""
    elif self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.ManagedSharedAcceleratorDeviceMemoryThroughTarch or self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.ManagedSharedAcceleratorDeviceMemoryThroughTarchWithoutDelete:
      d[ "CREATE_RECONSTRUCTED_PATCH" ] = """
    #if Dimensions==2
      double* oldQWithHalo = ::tarch::allocateMemory<double>(""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_2D" ] + """, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory);
    #elif Dimensions==3
      double* oldQWithHalo = ::tarch::allocateMemory<double>(""" + d[ "NUMBER_OF_DOUBLE_VALUES_IN_RECONSTRUCTED_PATCH_3D" ] + """, ::tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory);
    #endif
"""
    else:
      raise Exception( "Error: memory allocation mode {} for patch reconstruction not known".format(self.reconstructed_array_memory_location) )

    if self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.Heap:
      d[ "DESTROY_RECONSTRUCTED_PATCH" ] = """
    delete[] oldQWithHalo;
"""
    elif self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.HeapThroughTarch:
      d[ "DESTROY_RECONSTRUCTED_PATCH" ] = """
    ::tarch::freeMemory(oldQWithHalo, tarch::MemoryLocation::Heap );
"""
    elif self.reconstructed_array_memory_location==ReconstructedArrayMemoryLocation.ManagedSharedAcceleratorDeviceMemoryThroughTarch:
      d[ "DESTROY_RECONSTRUCTED_PATCH" ] = """
    ::tarch::freeMemory(oldQWithHalo, tarch::MemoryLocation::ManagedSharedAcceleratorDeviceMemory );
""" 
    else:
      d[ "DESTROY_RECONSTRUCTED_PATCH" ] = ""

  def get_constructor_body(self):
    return """  _treeNumber = treeNumber;
"""    

 
  def get_action_set_name(self):
    return __name__.replace(".py", "").replace(".", "_")


  def user_should_modify_template(self):
    return False


  _Template_TouchCellFirstTime_Preamble = """
  auto serialisePatchIndex = [](tarch::la::Vector<Dimensions,int> overlapCell, int normal) {{
    int base   = 1;
    int result = 0;
    for (int d=0; d<Dimensions; d++) {{
      result += overlapCell(d) * base;
      if (d==normal) {{
        base *= {OVERLAP}*2;
      }}
      else {{
        base *= {DOFS_PER_AXIS};
      }}
    }}
    return result;
  }};
  
  if ({GUARD}) {{
    logTraceInWith1Argument( "touchCellFirstTime(...)", marker.toString() );

    {CREATE_RECONSTRUCTED_PATCH}
"""


  _Template_TouchCellFirstTime_Fill_Patch = """
    //
    // Loop over original patch (k) and copy stuff over.
    //
    dfor(sourceCell,{DOFS_PER_AXIS}) {{
      tarch::la::Vector<Dimensions,int> destinationCell = sourceCell + tarch::la::Vector<Dimensions,int>({OVERLAP});
      int sourceCellSerialised       = peano4::utils::dLinearised(sourceCell,{DOFS_PER_AXIS});
      int destinationCellSerialised  = peano4::utils::dLinearised(destinationCell,{DOFS_PER_AXIS} + 2*{OVERLAP});
      for (int j=0; j<{UNKNOWNS}; j++) {{
        oldQWithHalo[destinationCellSerialised*{UNKNOWNS}+j] = {CELL_ACCESSOR}.value[ sourceCellSerialised*{UNKNOWNS}+j ];
        {ASSERTION_WITH_4_ARGUMENTS}( oldQWithHalo[destinationCellSerialised*{UNKNOWNS}+j]==oldQWithHalo[destinationCellSerialised*{UNKNOWNS}+j], sourceCell, j, _treeNumber, marker.toString() );
      }}
    }}
"""

  
  _Template_TouchCellFirstTime_Fill_Halos = """
    //
    // Bring in the auxiliary patches, i.e. befill halo
    //
    for(int d=0; d<Dimensions; d++) {{
      logTraceInWith1Argument( "touchCellFirstTime(...)::loopOverFace", d );
      //
      // d-loop over all dimensions except d. The vector k's entry d is set
      // to 0. We start with the left/bottom face, i.e. the one closer to the 
      // coordinate system's origin.
      //
      dfore(k,{DOFS_PER_AXIS},d,0) {{
        for (int i=0; i<{OVERLAP}; i++) {{
          tarch::la::Vector<Dimensions,int> destinationCell = k + tarch::la::Vector<Dimensions,int>({OVERLAP});
          tarch::la::Vector<Dimensions,int> sourceCell      = k;
          destinationCell(d) = i;
          sourceCell(d)      = i;
          
          int destinationCellSerialised   = peano4::utils::dLinearised(destinationCell,{DOFS_PER_AXIS} + 2*{OVERLAP});
          int sourceCellSerialised        = serialisePatchIndex(sourceCell,d);

          for (int j=0; j<{UNKNOWNS}; j++) {{
            oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ] = {FACES_ACCESSOR}(d).value[ sourceCellSerialised*{UNKNOWNS}+j ];
            {ASSERTION_WITH_7_ARGUMENTS}( 
              {ASSERTION_PREFIX_FOR_HALO} or 
              oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ]==oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ], 
              sourceCell, destinationCell, j, d, marker.toString(), _treeNumber, marker.toString()
            );
          }}

          destinationCell(d) = i+{DOFS_PER_AXIS}+{OVERLAP};
          sourceCell(d)      = i+{OVERLAP};

          destinationCellSerialised   = peano4::utils::dLinearised(destinationCell,{DOFS_PER_AXIS} + 2*{OVERLAP});
          sourceCellSerialised        = serialisePatchIndex(sourceCell,d);
          for (int j=0; j<{UNKNOWNS}; j++) {{
            oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ] = {FACES_ACCESSOR}(d+Dimensions).value[ sourceCellSerialised*{UNKNOWNS}+j ];
            {ASSERTION_WITH_7_ARGUMENTS}( 
              {ASSERTION_PREFIX_FOR_HALO} or 
              oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ]==oldQWithHalo[ destinationCellSerialised*{UNKNOWNS}+j ], 
              sourceCell, destinationCell, j, d, marker.toString(), _treeNumber, marker.toString()
            );
          }}
        }}
      }}
      logTraceOut( "touchCellFirstTime(...)::loopOverFace" );
    }}
"""


  _Template_TouchCellFirstTime_Core = """
    double* newQ = {CELL_ACCESSOR}.value;
    {CELL_FUNCTOR_IMPLEMENTATION}
    
    {DESTROY_RECONSTRUCTED_PATCH}
"""


  _Template_TouchCellFirstTime_Epilogue = """
    logTraceOut( "touchCellFirstTime(...)" );
  }}
"""


  def get_body_of_operation(self,operation_name):
    result = "\n"
    if operation_name==ActionSet.OPERATION_TOUCH_CELL_FIRST_TIME:
      d = {}
      self._add_action_set_entries_to_dictionary(d)
      result  = self._Template_TouchCellFirstTime_Preamble.format(**d)
      result += self._Template_TouchCellFirstTime_Fill_Patch.format(**d)
      result += self._Template_TouchCellFirstTime_Fill_Halos.format(**d)
      result += self._Template_TouchCellFirstTime_Core.format(**d)
      result += self._Template_TouchCellFirstTime_Epilogue.format(**d)
      pass 
    return result


  def get_attributes(self):
    return """int  _treeNumber;
"""   


  def get_includes(self):
    return """
#include <functional>
#include "peano4/utils/Loop.h"
#include "tarch/multicore/Core.h"
"""
