# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from peano4.solversteps.ActionSet import ActionSet


import jinja2


class PlotGridInPeanoBlockFormat(ActionSet):
  NoMetaFile     = "no-meta-file"
  CountTimeSteps = "count-time-steps"
  
  def __init__(self, filename, cell_unknown=None, time_stamp_evaluation=NoMetaFile, guard_predicate="true", additional_includes=""):
    """
      Plot only the grid structure
      
      :Attibutes:
      
      filename: String
         Name of the output file.
         
      cell_unknown: None or cell data   
         If you use cell unknowns, pass any unknown in. As we do not dump
         any semantic information about unknowns, it does not matter which 
         one you choose. If you don't have cell unknowns at all, pass in 
         None.
         
      time_stamp_evaluation: String
         C++ expression returning a double. You can also hand in NoMetaFile
         or CountTimeSteps.         
          
    """
    super( PlotGridInPeanoBlockFormat, self ).__init__()
    
    self.d = {}
    self.d[ "FILENAME" ]     = filename
    self.d[ "CELL_WRAPPER" ] = "Cell"
    if cell_unknown!=None:
      self.d[ "CELL_WRAPPER" ] += cell_unknown.name
    self.d[ "GUARD_PREDICATE" ]    = guard_predicate
    self.additional_includes       = additional_includes
    self.d[ "TIMESTAMP" ]          = time_stamp_evaluation
    

  __Template_Constructor = """
  _writer      = nullptr;
  _dataWriter  = nullptr;
  _treeNumber  = treeNumber;

  // An MPI lock (critical section) would be important!
    
  logDebug( "PlotGrid2PlotGridInPeanoBlockFormat1()", "created tree instance for " << treeNumber );
"""


  def get_constructor_body(self):
    return self.__Template_Constructor.format(**self.d)


  __Template_EndTraversal = jinja2.Template("""
  assertion1(_dataWriter!=nullptr,_treeNumber);

  _dataWriter->close();
  _writer->writeToFile();
  
  delete _dataWriter;
  delete _writer;

  _dataWriter = nullptr;
  _writer     = nullptr;
""")

    
  def get_destructor_body(self):
    return ""


  def get_body_of_getGridControlEvents(self):
    return "  return std::vector< peano4::grid::GridControlEvent >();\n" 


  def get_action_set_name(self):
    return __name__.replace(".py", "").replace(".", "_")


  def user_should_modify_template(self):
    return False


  __Template_TouchCellFirstTime = jinja2.Template(""" 
  if ( {{GUARD_PREDICATE}} ) {
  int vertexIndices[TwoPowerD];

  int indices = _writer->plotPatch(
    marker.x() - marker.h() * 0.5,
    marker.h()
  );

  assertion( _dataWriter!=nullptr );
  
  double markerData[] = {
    marker.hasBeenRefined() ? 1.0 : 0.0,
    marker.willBeRefined() ? 1.0 : 0.0,
    marker.isLocal() ? 1.0 : 0.0,
    marker.hasBeenEnclaveCell() ? 1.0 : 0.0,
    marker.willBeEnclaveCell() ? 1.0 : 0.0
  };
 _dataWriter->plotCell(indices,markerData);
 }
""")


  __Template_BeginTraversal = jinja2.Template("""
  tarch::mpi::Lock lock( _semaphore );
  
  static int counter = -1;
  counter++;

  std::ostringstream snapshotFileName;
  snapshotFileName << "{{FILENAME}}-" << counter;

  if (tarch::mpi::Rank::getInstance().getNumberOfRanks()>1 ) {
    snapshotFileName << "-rank-" << tarch::mpi::Rank::getInstance().getRank();
  }

  {% if TIMESTAMP==\"""" + NoMetaFile + """\" %}
  _writer = new tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter(
    Dimensions, snapshotFileName.str(), "{{FILENAME}}",
    tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter::IndexFileMode::NoIndexFile,
    0.0
  );
  {% elif TIMESTAMP==\"""" + CountTimeSteps + """\" %}
  static int timeStep = -1;
  timeStep++;
  _writer = new tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter(
    Dimensions, snapshotFileName.str(), "{{FILENAME}}",
    tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter::IndexFileMode::AppendNewData,
    timeStep
  );
  {% else %}
  _writer = new tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter(
    Dimensions, snapshotFileName.str(), "{{FILENAME}}",
    tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter::IndexFileMode::AppendNewData,
    {{TIMESTAMP}}
  );    
  {% endif %}
      
  _dataWriter = _writer->createCellDataWriter( "cell-marker", 1, 5, "is-refined,has-been-refined,local,enclave" );
""")


  def get_body_of_operation(self,operation_name):
    result = "\n"
    if operation_name==ActionSet.OPERATION_TOUCH_CELL_FIRST_TIME:
      result = self.__Template_TouchCellFirstTime.render(**self.d) 
    if operation_name==ActionSet.OPERATION_BEGIN_TRAVERSAL:
      result = self.__Template_BeginTraversal.render(**self.d)             
    if operation_name==ActionSet.OPERATION_END_TRAVERSAL:
      result = self.__Template_EndTraversal.render(**self.d)             
    return result


  def get_attributes(self):
    return """
    static tarch::mpi::BooleanSemaphore                                              _semaphore;
    
    int                _treeNumber;

    tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter*                  _writer;
    tarch::plotter::griddata::blockstructured::PeanoTextPatchFileWriter::CellDataWriter*  _dataWriter;
"""


  def get_includes(self):
    return """
#include "tarch/plotter/griddata/blockstructured/PeanoTextPatchFileWriter.h"
#include "tarch/multicore/Lock.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/mpi/Lock.h"
#include "tarch/mpi/BooleanSemaphore.h"
#include "peano4/parallel/SpacetreeSet.h"
""" + self.additional_includes


  def get_static_initialisations(self,full_qualified_classname):
    return """
tarch::mpi::BooleanSemaphore  """ + full_qualified_classname + """::_semaphore;
"""
    


