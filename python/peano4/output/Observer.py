# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import os
import re
import jinja2


from .Helper import write_file
from .Helper import using_cuda_backend

from .Overwrite import Overwrite

from .Jinja2TemplatedHeaderFile import Jinja2TemplatedHeaderFile


class Observer(object):
    default_overwrite = True

    def __init__(
        self,
        classname,
        namespace,
        subdirectory,
        included_actions,
        vertices,
        faces,
        cells,
    ):
        """
        Included actions is a list of qualified actions which are used
        """
        self.classname = classname
        self.namespace = namespace
        self.subdirectory = subdirectory

        self.included_actions = included_actions
        self.vertices = vertices
        self.faces = faces
        self.cells = cells

        self.d = {}
        self.d["OPEN_NAMESPACE"] = ""
        self.d["CLOSE_NAMESPACE"] = ""
        self.d["FULL_QUALIFIED_CLASSNAME"] = ""
        for i in namespace:
            self.d["OPEN_NAMESPACE"] += "namespace " + i + "{\n"
            self.d["CLOSE_NAMESPACE"] += "}\n"
            self.d["FULL_QUALIFIED_CLASSNAME"] += i + "::"
        self.d["CLASSNAME"] = classname
        self.d["FULL_QUALIFIED_CLASSNAME"] += classname
        self.d["INCLUDE_GUARD"] = (
            "_" + self.d["FULL_QUALIFIED_CLASSNAME"].replace("::", "_").upper() + "_H_"
        )
        self.d["INCLUDES"] = ""
        self.d["ATTRIBUTES"] = ""

        for action in self.included_actions:
            self.d["INCLUDES"] += '#include "'
            self.d["INCLUDES"] += action.replace("::", "/")
            self.d["INCLUDES"] += '.h"\n'
            self.d["ATTRIBUTES"] += "    "
            self.d["ATTRIBUTES"] += action
            self.d["ATTRIBUTES"] += "   _actionSet"
            self.d["ATTRIBUTES"] += str(self.included_actions.index(action))
            self.d["ATTRIBUTES"] += "; \n"

        self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] = ""
        for cell in cells:
            if cells.index(cell) != 0:
                self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS"] += ","
                self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS"] += ","
                self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] += ","
                self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] += ","
            self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS"] += (
                "repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(0)"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS"] += (
                "repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(1)"
            )
            self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] += (
                "repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(0)"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_PICK_ENTRY"] += (
                "repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(1)"
            )

        self.d["MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_CELL_EVENT"] = self.d[
            "MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS"
        ]
        self.d["MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_CELL_EVENT"] = self.d[
            "MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS"
        ]

        self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_PICK_ENTRY"] = ""
        for face in faces:
            if faces.index(face) != 0:
                self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS"] += ","
                self.d["MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS"] += ","
                self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY"] += ","
                self.d[
                    "MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_PICK_ENTRY"
                ] += ","
            self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS"] += (
                "peano4::datamanagement::FaceEnumerator<"
                + face.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoTimesD-1) )"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS"] += (
                "peano4::datamanagement::FaceEnumerator<"
                + face.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoTimesD*2-1) )"
            )
            self.d["MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY"] += (
                "peano4::datamanagement::FaceEnumerator<"
                + face.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoTimesD-1) )(pick)"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_PICK_ENTRY"] += (
                "peano4::datamanagement::FaceEnumerator<"
                + face.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoTimesD*2-1) )(pick)"
            )

        self.d["MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS"] = ""
        self.d["MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"] = ""
        self.d["MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"] = ""
        for vertex in vertices:
            if vertices.index(vertex) != 0:
                self.d["MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS"] += ","
                self.d["MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS"] += ","
                self.d[
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"
                ] += ","
                self.d[
                    "MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"
                ] += ","
            self.d["MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS"] += (
                "peano4::datamanagement::VertexEnumerator<"
                + vertex.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoPowerD-1) )"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS"] += (
                "peano4::datamanagement::VertexEnumerator<"
                + vertex.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoPowerD*2-1) )"
            )
            self.d["MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"] += (
                "peano4::datamanagement::VertexEnumerator<"
                + vertex.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoPowerD-1) )(pick)"
            )
            self.d["MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY"] += (
                "peano4::datamanagement::VertexEnumerator<"
                + vertex.get_full_qualified_type()
                + ">( &repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top(TwoPowerD*2-1) )(pick)"
            )

    def __generate_header(self, overwrite, directory):
        headerfile_template = (
            os.path.realpath(__file__)
            .replace(".pyc", ".template.h")
            .replace(".py", ".template.h")
        )
        md = self.mkSubDict(["ATTRIBUTES", "INCLUDES"])
        header = Jinja2TemplatedHeaderFile(
            headerfile_template,
            self.classname,
            self.namespace,
            self.subdirectory,
            md,
            self.default_overwrite,
        )
        header.generate(overwrite, directory)
        del header

    TemplateConstructor = """

{{FULL_QUALIFIED_CLASSNAME}}::{{CLASSNAME}}(int spacetreeId):
  _spacetreeId( spacetreeId ) {{MAPPING_INITIALISATION_LIST}}
{}


  """

    def __generate_constructor(self, output_file):
        #
        # Constructor
        #
        self.d["MAPPING_INITIALISATION_LIST"] = ""
        for actions in self.included_actions:
            self.d["MAPPING_INITIALISATION_LIST"] += ", _actionSet"
            self.d["MAPPING_INITIALISATION_LIST"] += str(
                self.included_actions.index(actions)
            )
            self.d["MAPPING_INITIALISATION_LIST"] += "(spacetreeId)"
        output_file.write(jinja2.Template(self.TemplateConstructor).render(**self.d))

    TemplateClone = """

peano4::grid::TraversalObserver* {{FULL_QUALIFIED_CLASSNAME}}::clone(int spacetreeId) {
  return new {{CLASSNAME}}(spacetreeId);
}

  """

    TemplateBeginTraversal = """

void {{FULL_QUALIFIED_CLASSNAME}}::beginTraversal( const tarch::la::Vector<Dimensions,double>&  x, const tarch::la::Vector<Dimensions,double>&  h ) {
  logTraceInWith2Arguments( "beginTraversal(...)", x, h );
  //
  // Invoke beginTraversal() on the actions
  //
{{MAPPING_BEGIN_TRAVERSAL_CALLS}}

  //
  // Fill call stacks with dummies which represent level 0 such that we can
  // call standard action routines on level 1 with parents. Without these
  // statements, a top(1) call would raise an assertion.
  //
{{INITIAL_PUSH_TO_OUTPUT_STREAMS}}
  logTraceOutWith2Arguments( "beginTraversal(...)", x, h );
}

  """

    def __generate_beginTraversal(self, output_file):
        self.d["MAPPING_BEGIN_TRAVERSAL_CALLS"] = ""
        for action in self.included_actions:
            self.d["MAPPING_BEGIN_TRAVERSAL_CALLS"] += "  _actionSet"
            self.d["MAPPING_BEGIN_TRAVERSAL_CALLS"] += str(
                self.included_actions.index(action)
            )
            self.d["MAPPING_BEGIN_TRAVERSAL_CALLS"] += ".beginTraversal();\n"

        self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] = ""
        for cell in self.cells:
            self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] += (
                "  repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->push( "
                + cell.get_full_qualified_type()
                + "() );\n"
            )
            pass

        for face in self.faces:
            self.d[
                "INITIAL_PUSH_TO_OUTPUT_STREAMS"
            ] += "  for (int i=0; i<TwoTimesD; i++) {\n"
            self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] += (
                "    repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->push( "
                + face.get_full_qualified_type()
                + "() );\n"
            )
            self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] += "  }\n"
            pass

        for vertex in self.vertices:
            self.d[
                "INITIAL_PUSH_TO_OUTPUT_STREAMS"
            ] += "  for (int i=0; i<TwoPowerD; i++) {\n"
            self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] += (
                "    repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->push( "
                + vertex.get_full_qualified_type()
                + "() );\n"
            )
            self.d["INITIAL_PUSH_TO_OUTPUT_STREAMS"] += "  }\n"
            pass

        output_file.write(jinja2.Template(self.TemplateBeginTraversal).render(**self.d))

    TemplateEndTraversal = """

void {{FULL_QUALIFIED_CLASSNAME}}::endTraversal( const tarch::la::Vector<Dimensions,double>&  x, const tarch::la::Vector<Dimensions,double>&  h ) {
  logTraceInWith2Arguments( "endTraversal(...)", x, h );
  {{MAPPING_END_TRAVERSAL_CALLS}}
  {{FINAL_POP_FROM_INPUT_STREAMS}}
  logTraceOutWith2Arguments( "endTraversal(...)", x, h );
}

  """

    def __generate_endTraversal(self, output_file):
        self.d["MAPPING_END_TRAVERSAL_CALLS"] = ""
        for action in self.included_actions:
            self.d["MAPPING_END_TRAVERSAL_CALLS"] += "  _actionSet"
            self.d["MAPPING_END_TRAVERSAL_CALLS"] += str(
                self.included_actions.index(action)
            )
            self.d["MAPPING_END_TRAVERSAL_CALLS"] += ".endTraversal();\n"

        self.d["FINAL_POP_FROM_INPUT_STREAMS"] = ""
        for cell in self.cells:
            self.d["FINAL_POP_FROM_INPUT_STREAMS"] += (
                "  repositories::DataRepository::_"
                + cell.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pop();\n"
            )
            pass

        for face in self.faces:
            self.d[
                "FINAL_POP_FROM_INPUT_STREAMS"
            ] += "  for (int i=0; i<TwoTimesD; i++) {\n"
            self.d["FINAL_POP_FROM_INPUT_STREAMS"] += (
                "    repositories::DataRepository::_"
                + face.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pop();\n"
            )
            self.d["FINAL_POP_FROM_INPUT_STREAMS"] += "  }\n"
            pass

        for vertex in self.vertices:
            self.d[
                "FINAL_POP_FROM_INPUT_STREAMS"
            ] += "  for (int i=0; i<TwoPowerD; i++) {\n"
            self.d["FINAL_POP_FROM_INPUT_STREAMS"] += (
                "    repositories::DataRepository::_"
                + vertex.get_logical_type_name()
                + "Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pop();\n"
            )
            self.d["FINAL_POP_FROM_INPUT_STREAMS"] += "  }\n"
            pass

        output_file.write(jinja2.Template(self.TemplateEndTraversal).render(**self.d))

    TemplatePrepareTraversal = """

void {{FULL_QUALIFIED_CLASSNAME}}::prepareTraversal() {
  logTraceIn( "prepareTraversal(...)" );
  {{MAPPING_PREPARE_TRAVERSAL_CALLS}}
  logTraceOut( "prepareTraversal(...)" );
}

  """

    def __generate_prepareTraversal(self, output_file):
        self.d["MAPPING_PREPARE_TRAVERSAL_CALLS"] = ""
        for action in self.included_actions:
            self.d["MAPPING_PREPARE_TRAVERSAL_CALLS"] += "  "
            self.d["MAPPING_PREPARE_TRAVERSAL_CALLS"] += action
            self.d["MAPPING_PREPARE_TRAVERSAL_CALLS"] += "::prepareTraversal();\n"
        output_file.write(
            jinja2.Template(self.TemplatePrepareTraversal).render(**self.d)
        )

    TemplateUnprepareTraversal = """

void {{FULL_QUALIFIED_CLASSNAME}}::unprepareTraversal() {
  logTraceIn( "unprepareTraversal(...)" );
  {{MAPPING_UNPREPARE_TRAVERSAL_CALLS}}
  logTraceOut( "unprepareTraversal(...)" );
}

  """

    def __generate_unprepareTraversal(self, output_file):
        self.d["MAPPING_UNPREPARE_TRAVERSAL_CALLS"] = ""
        for action in self.included_actions:
            self.d["MAPPING_UNPREPARE_TRAVERSAL_CALLS"] += "  "
            self.d["MAPPING_UNPREPARE_TRAVERSAL_CALLS"] += action
            self.d["MAPPING_UNPREPARE_TRAVERSAL_CALLS"] += "::unprepareTraversal();\n"
        output_file.write(
            jinja2.Template(self.TemplateUnprepareTraversal).render(**self.d)
        )

    def __format_template_per_action(
        self, output_file, template, reverse_order=False, manual_dict=None
    ):
        """

        Takes the specified template file, iterates over actions and pastes
        the template into the output file once per action. Per action, the dictionary's
        entries are updated. Otherwise, the dictionary remains unchanged.


        output_file:
          Handle on output file

        """
        local_actions = [x for x in self.included_actions]
        if reverse_order:
            local_actions.reverse()

        for action in local_actions:
            if manual_dict is None:
                self.d["ACTIVE_ACTION_SET"] = "_actionSet" + str(
                    self.included_actions.index(action)
                )
                self.d["ACTIVE_ACTION_SET_FULL_QUALIFIED_NAME"] = action
                if output_file is not None:
                    output_file.write(jinja2.Template(template).render(**self.d))
                else:
                    return jinja2.Template(template).render(**self.d)
            else:
                manual_dict["ACTIVE_ACTION_SET"] = "_actionSet" + str(
                    self.included_actions.index(action)
                )
                manual_dict["ACTIVE_ACTION_SET_FULL_QUALIFIED_NAME"] = action
                if output_file is not None:
                    output_file.write(jinja2.Template(template).render(**manual_dict))
                else:
                    return jinja2.Template(template).render(**manual_dict)

    TemplateGetGridControlEvents_Prologue = """

std::vector< peano4::grid::GridControlEvent > {{FULL_QUALIFIED_CLASSNAME}}::getGridControlEvents() const {
  std::vector< peano4::grid::GridControlEvent > result;
"""

    TemplateGetGridControlEvents_MappingCall = """
  {
    const std::vector< peano4::grid::GridControlEvent > actionResult = {{ACTIVE_ACTION_SET}}.getGridControlEvents();
    result.insert(result.begin(),actionResult.begin(),actionResult.end());
  }
"""

    TemplateGetGridControlEvents_Epilogue = """
  return result;
}


"""

    def __generate_getGridControlEvents(self, output_file):
        output_file.write(
            jinja2.Template(self.TemplateGetGridControlEvents_Prologue).render(
                **{"FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"]}
            )
        )
        self.__format_template_per_action(
            output_file, self.TemplateGetGridControlEvents_MappingCall, manual_dict={}
        )
        output_file.write(
            jinja2.Template(self.TemplateGetGridControlEvents_Epilogue).render(**{})
        )

    def __generate_clone(self, output_file):
        output_file.write(jinja2.Template(self.TemplateClone).render(**self.d))

    TemplateEnterCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::enterCell( const peano4::grid::GridTraversalEvent&  event ) {
  logTraceInWith2Arguments( "enterCell(...)", _spacetreeId, event.toString() );
"""

    TemplateLoadCell_VertexLoad = """
  // Load vertex {{logical_type_name}}
  std::function<bool ()> loadVertex{{logical_type_name}} = 
    [&]()->bool {
    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pushBlock( TwoPowerD );
    for (int i=0; i<TwoPowerD; i++) {
      int inVertexStack          = event.getVertexDataFrom(i);
      int outVertexStackPosition = event.getVertexDataTo(i);
      logDebug("loadCell(...)", "vertex stack " << inVertexStack << "->pos-" << outVertexStackPosition );

      peano4::datamanagement::VertexMarker  marker(event,outVertexStackPosition);

      bool dataResidesOnIntermediateStack = 
        not peano4::grid::PeanoCurve::isInOutStack(inVertexStack)
        and
        inVertexStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        inVertexStack!=peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        inVertexStack!=peano4::grid::TraversalObserver::NoData;
      bool dataResidesOnInputStack        = 
        peano4::grid::PeanoCurve::isInOutStack(inVertexStack)
        and
        ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ));
      bool initialiseStackContentForNewData = 
        ::peano4::grid::computeOnData( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ))
        or
        inVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        or
        inVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        ;

      if ( dataResidesOnIntermediateStack or dataResidesOnInputStack ) {
        logDebug( 
          "loadCell(...)", 
          "load data of vertex {{name}} with " << marker.x(outVertexStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        assertion4( not repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inVertexStack))->empty(), event.toString(), peano4::datamanagement::VertexMarker(event).toString(), _spacetreeId, inVertexStack);
        {{full_qualified_type}} data = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inVertexStack))->pop() );
        #if PeanoDebug>0
        if ( peano4::grid::PeanoCurve::isInOutStack(inVertexStack) ) {
          assertionVectorNumericalEquals7( data.getDebugX(), marker.x(outVertexStackPosition), event.toString(), data.getDebugX(), marker.toString(), inVertexStack, i, outVertexStackPosition, _spacetreeId );
          assertionVectorNumericalEquals6( data.getDebugH(), marker.h(),                       event.toString(), data.getDebugX(), marker.toString(), inVertexStack, i, _spacetreeId );
        }
        #endif
        view.set(outVertexStackPosition,data);
      }
      else if ( initialiseStackContentForNewData ) {
        {{full_qualified_type}} data;
        #if PeanoDebug>0
        logDebug( 
          "loadCell(...)", 
          "initialise meta data of new vertex {{name}} with " << marker.x(outVertexStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        data.setDebugX( marker.x(outVertexStackPosition) );
        data.setDebugH( marker.h() );
        #endif
        view.set(outVertexStackPosition,data);
      }
      else {
        #if PeanoDebug>0
        logDebug( 
          "loadCell(...)", 
          "initialise meta data of unused vertex {{name}} with " << marker.x(outVertexStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        view.get(outVertexStackPosition).setDebugX( marker.x(outVertexStackPosition) );
        view.get(outVertexStackPosition).setDebugH( marker.h() );
        #endif
      }
    }
    return false;
  };
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    loadVertex{{logical_type_name}}
  ));
  
"""

    TemplateEnterCell_VertexMappingCall = """
  // Handle vertex {{logical_type_name}}
  {
    peano4::datamanagement::VertexMarker  marker(event);

    for (int i=0; i<TwoPowerD; i++) {
      int inVertexStack  = event.getVertexDataFrom(i);
      int pick          = event.getVertexDataTo(i);   // the vertex position

      marker.select(pick);

      logDebug( "enterCell(...)", inVertexStack << "->" << pick << " (is-local=" << marker.isLocal() << ")" );

      if (
        inVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        marker.isLocal()
      ) {
        // Take care about the coarse grid accesses: Faces and cells are not yet loaded.
        // Therefore we don't use the usual shift of @f$ 2 \cdot 2d @f$ or @f$ 2 \cdot 2^d @f$
        // but only half of it.
        {{ACTIVE_ACTION_SET}}.createPersistentVertex(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
        {{ACTIVE_ACTION_SET}}.touchVertexFirstTime(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        inVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.createHangingVertex(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        peano4::grid::PeanoCurve::isInOutStack(inVertexStack)
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchVertexFirstTime(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
    }
  }
"""

    TemplateLoadCell_FaceLoad = """
  // Load face {{logical_type_name}}
  std::function<bool ()> loadFace{{logical_type_name}} = 
    [&]()->bool {
    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pushBlock( TwoTimesD );
    for (int i=0; i<TwoTimesD; i++) {
      int inFaceStack          = event.getFaceDataFrom(i);
      int outFaceStackPosition = event.getFaceDataTo(i);
      logDebug("loadCell(...)", "face stack " << inFaceStack << "->pos-" << outFaceStackPosition );

      peano4::datamanagement::FaceMarker  marker(event,outFaceStackPosition);

      bool dataResidesOnIntermediateStack = 
        not peano4::grid::PeanoCurve::isInOutStack(inFaceStack)
        and
        inFaceStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        inFaceStack!=peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        inFaceStack!=peano4::grid::TraversalObserver::NoData;
      bool dataResidesOnInputStack        = 
        peano4::grid::PeanoCurve::isInOutStack(inFaceStack)
        and
        ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ));
      bool initialiseStackContentForNewData = 
        ::peano4::grid::computeOnData( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ))
        or
        inFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        or
        inFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        ;

      if ( dataResidesOnIntermediateStack or dataResidesOnInputStack ) {
        logDebug( 
          "loadCell(...)", 
          "load data of face {{name}} with " << marker.x(outFaceStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        assertion5( 
          not repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inFaceStack))->empty(), 
          event.toString(), 
          peano4::datamanagement::FaceMarker(event).toString(), 
          _spacetreeId, inFaceStack,
          ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
                        marker
                        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
          ))
        );
        {{full_qualified_type}} data = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inFaceStack))->pop() );
        #if PeanoDebug>0
        if ( peano4::grid::PeanoCurve::isInOutStack(inFaceStack) ) {
          assertionVectorNumericalEquals5( data.getDebugX(), marker.x(outFaceStackPosition), data.getDebugX(), data.getDebugH(), marker.toString(), outFaceStackPosition, _spacetreeId );
          assertionVectorNumericalEquals5( data.getDebugH(), marker.h(),                     data.getDebugX(), data.getDebugH(), marker.toString(), outFaceStackPosition, _spacetreeId );
        }
        #endif
        view.set(outFaceStackPosition,data);
      }
      else if ( initialiseStackContentForNewData ) {
        {{full_qualified_type}} data;
        #if PeanoDebug>0
        logDebug( 
          "loadCell(...)", 
          "initialise meta data of new face {{name}} with " << marker.x(outFaceStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        data.setDebugX( marker.x(outFaceStackPosition) );
        data.setDebugH( marker.h() );
        #endif
        view.set(outFaceStackPosition,data);
      }
      else {
        #if PeanoDebug>0
        logDebug( 
          "loadCell(...)", 
          "initialise meta data of unused face {{name}} with " << marker.x(outFaceStackPosition) << " x " << marker.h() << 
          ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
              marker
              {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
            ))
        );
        view.get(outFaceStackPosition).setDebugX( marker.x(outFaceStackPosition) );
        view.get(outFaceStackPosition).setDebugH( marker.h() );
        #endif
      }
    }
    return false;
  };
  
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    loadFace{{logical_type_name}}
  ));
  
"""

    TemplateEnterCell_FaceMappingCall = """
  // Handle face {{logical_type_name}}
  {
    peano4::datamanagement::FaceMarker marker( event );
    for (int i=0; i<TwoTimesD; i++) {
      int inFaceStack = event.getFaceDataFrom(i);
      int pick        = event.getFaceDataTo(i);

      marker.select(pick);

      assertion4( 
        marker.isLocal() or not event.getIsCellLocal(), 
        marker.toString(), 
        event.toString(), 
        i, 
        _spacetreeId 
      );

      if (
        inFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.createPersistentFace(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
        {{ACTIVE_ACTION_SET}}.touchFaceFirstTime(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        inFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.createHangingFace(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        peano4::grid::PeanoCurve::isInOutStack(inFaceStack)
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchFaceFirstTime(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
    }
  }
"""

    TemplateLoadCell_CellLoad = """
  // Load cell {{logical_type_name}}
  std::function<bool ()> loadCell{{logical_type_name}} =  
    [&]()->bool {
    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->pushBlock( 1 );

    peano4::datamanagement::CellMarker  marker(event);

    const int inCellStack  = event.getCellData();
    const int outCellStack = peano4::grid::PeanoCurve::CallStack;

    bool dataResidesOnInputStack        = 
      inCellStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
      and
      inCellStack!=peano4::grid::TraversalObserver::NoData
      and
      ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ));
    bool initialiseStackContentForNewData = 
      ::peano4::grid::computeOnData( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      or
      inCellStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
      ;
 
    if ( dataResidesOnInputStack ) {
      logDebug( 
        "loadCell(...)", 
        "load data of cell {{name}} with " << marker.toString() << 
        ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
           marker
           {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ))
      );
      assertion3( not repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inCellStack))->empty(), event.toString(), _spacetreeId, inCellStack);
      {{full_qualified_type}} data = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,inCellStack))->pop() );

      #if PeanoDebug>0
      assertionVectorNumericalEquals4( data.getDebugX(), marker.x(), data.getDebugX(), data.getDebugH(), marker.toString(), _spacetreeId );
      assertionVectorNumericalEquals4( data.getDebugH(), marker.h(), data.getDebugX(), data.getDebugH(), marker.toString(), _spacetreeId );
      #endif

      view.set(0,data);
    }
    else if ( initialiseStackContentForNewData ) {
      {{full_qualified_type}} data;

      #if PeanoDebug>0
      logDebug( 
        "loadCell(...)", 
        "initialise meta data of new cell {{name}} with " << marker.toString() << 
        ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
            marker
            {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
          ))
      );
      data.setDebugX( marker.x() );
      data.setDebugH( marker.h() );
      #endif
      view.set(0,data);
    }
    else {
      #if PeanoDebug>0
      logDebug( 
        "loadCell(...)", 
        "initialise meta data of unused cell {{name}} with " << marker.toString() << 
        ": " << ::peano4::grid::toString( {{full_qualified_type}}::loadStoreComputeFlag(
            marker
            {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
          ))
      );
      view.get(0).setDebugX( marker.x() );
      view.get(0).setDebugH( marker.h() );
      #endif
    }
    return false;
  };
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    loadCell{{logical_type_name}}
  ));
  
"""

    TemplateEnterCell_CellMappingCall = """
  // Invoke creational events on cell {{logical_type_name}}
  {
    peano4::datamanagement::CellMarker marker( event );
    if (
      event.getCellData()==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
      and
      marker.isLocal()
    ) {
      {{ACTIVE_ACTION_SET}}.createCell(
        marker,
        {{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_COMMA_}}
        {{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_COMMA_}}
        {{MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
      );
    }
  }
"""

    TemplateEnterCell_MappingCall = """
  {
    peano4::datamanagement::CellMarker marker( event );
    if (
      marker.isLocal()
    ) {
      {{ACTIVE_ACTION_SET}}.touchCellFirstTime(
         marker
        {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_CELL_EVENT}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_CELL_EVENT}}
      );
    }
  }
"""

    TemplateEnterCell_Epilogue = """
  logTraceOut( "enterCell(...)" );
}


"""

    def generateDictEntry(self, key):
        """
        Some logic to produce new dictionary entries from existing ones.
        This deals mainly with having to avoid C++ code with consequtive commata
        as some dict entries can be an empty string.
        """

        fields = key.split(",")
        nonempty = [self.d[f] for f in fields if not f == "" and not self.d[f] == ""]
        if len(nonempty) > 0:
            s = ",".join(nonempty)
            if key.startswith(",") and not s.startswith(","):
                s = "," + s
            if key.endswith(",") and not s.endswith(","):
                s += ","
        else:
            s = ""

        return s

    def mkSubDict(self, keys):
        """

        Create the particular subdictionary that's used for an expression.

        """
        temp = {}
        for k in keys:
            if not k in self.d:
                temp[k] = self.generateDictEntry(k)
            else:
                temp[k] = self.d[k]
        result = {}
        for k in temp:
            result[k.replace(",", "_COMMA_").strip()] = temp[k]
        return result

    def __generate_loadCell(self, output_file):
        output_file.write(
            jinja2.Template(self.TemplateLoadCell_Prologue).render(**self.d)
        )

        for vertex in self.vertices:
            temp = {
                "name": vertex.name,
                "enumeration_type": vertex.get_enumeration_type(),
                "logical_type_name": vertex.get_logical_type_name(),
                "full_qualified_type": vertex.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": vertex.additional_load_and_store_arguments,
            }
            self.d["name"] = vertex.name
            output_file.write(
                jinja2.Template(self.TemplateLoadCell_VertexLoad).render(**temp)
            )

        for face in self.faces:
            temp = {
                "name": face.name,
                "enumeration_type": face.get_enumeration_type(),
                "logical_type_name": face.get_logical_type_name(),
                "full_qualified_type": face.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": face.additional_load_and_store_arguments,
            }
            self.d["name"] = face.name
            output_file.write(
                jinja2.Template(self.TemplateLoadCell_FaceLoad).render(**temp)
            )

        for cell in self.cells:
            temp = {
                "name": cell.name,
                "logical_type_name": cell.get_logical_type_name(),
                "full_qualified_type": cell.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": cell.additional_load_and_store_arguments,
            }
            self.d["name"] = cell.name
            output_file.write(
                jinja2.Template(self.TemplateLoadCell_CellLoad).render(**temp)
            )

        output_file.write(jinja2.Template(self.TemplateLoadCell_Epilogue).render({}))

    def __generate_enterCell(self, output_file):
        """
        Generates enter cell
        """
        output_file.write(
            jinja2.Template(self.TemplateEnterCell_Prologue).render(**self.d)
        )

        if len(self.vertices) > 0:
            md = self.mkSubDict(
                [
                    "name",
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateEnterCell_VertexMappingCall,
                False,
                manual_dict=md,
            )

        if len(self.faces) > 0:
            md = self.mkSubDict(
                [
                    "name",
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY",
                    "MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS",
                    "MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                    ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateEnterCell_FaceMappingCall,
                False,
                manual_dict=md,
            )

        if len(self.cells) > 0:
            md = self.mkSubDict(
                [
                    "name",
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS,",
                    "MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS,",
                    "MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS",
                    ",MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS",
                    ",MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateEnterCell_CellMappingCall,
                False,
                manual_dict=md,
            )

        md = self.mkSubDict(
            [
                ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_CELL_EVENT",
                ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_CELL_EVENT",
            ]
        )

        self.__format_template_per_action(
            output_file, self.TemplateEnterCell_MappingCall, False, manual_dict=md
        )
        output_file.write(jinja2.Template(self.TemplateEnterCell_Epilogue).render({}))

    TemplateLoadCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::loadCell( const peano4::grid::GridTraversalEvent&  event ) {
  logTraceInWith2Arguments( "loadCell(...)", _spacetreeId, event.toString() );
  
  std::vector< tarch::multicore::Task* >  tasks;
"""

    TemplateStoreCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::storeCell( const peano4::grid::GridTraversalEvent&  event ) {
  logTraceInWith2Arguments( "storeCell(...)", _spacetreeId, event.toString() );
  std::vector< tarch::multicore::Task* >  tasks;
"""

    TemplateLeaveCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::leaveCell( const peano4::grid::GridTraversalEvent&  event ) {
  logTraceInWith2Arguments( "leaveCell(...)", _spacetreeId, event.toString() );
"""

    TemplateLeaveCell_MappingCall = """
  {
    peano4::datamanagement::CellMarker marker( event );
    if (
      marker.isLocal()
    ) {
      {{ACTIVE_ACTION_SET}}.touchCellLastTime(
         marker
        {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_CELL_EVENT}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_CELL_EVENT}}
      );
    }
  }
"""

    TemplateLeaveCell_CellStore_MappingCall = """
  {
    peano4::datamanagement::CellMarker marker( event );
    if (
      event.getCellData()==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
      and
      marker.isLocal()
    ) {
    {{ACTIVE_ACTION_SET}}.destroyCell(
      marker
        {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS}}
        ,{{MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY}}
        {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
      );
    }
  }
"""

    TemplateStoreCell_CellStore = """
  // Handle cell {{logical_type_name}}
  std::function<bool ()> storeCell{{logical_type_name}} = 
    [&]()->bool {
    const int inCellStack   = peano4::grid::PeanoCurve::CallStack;
    const int outCellStack  = event.getCellData();
    logDebug("storeCell(...)", "cell stack " << inCellStack << "->pos-" << outCellStack );

    peano4::datamanagement::CellMarker  marker(event);

    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->popBlock( 1 );

    bool dataShallBePushedOntoIntermediateStack = 
        not peano4::grid::PeanoCurve::isInOutStack(outCellStack)
        and
        outCellStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        outCellStack!=peano4::grid::TraversalObserver::NoData;
    bool dataShallBePushedOntoOutputStack        = 
        peano4::grid::PeanoCurve::isInOutStack(outCellStack)
        and
        ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ));

    if ( dataShallBePushedOntoIntermediateStack or dataShallBePushedOntoOutputStack ) {
      repositories::DataRepository::_{{logical_type_name}}Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,outCellStack))->push( view.get(0) );
    }
    else {
      logDebug(
        "storeCell(...)",
        "do not store cell {{name}} with " << marker.x() << " x " << marker.h() << ":" <<
        "  destroy=" << (outCellStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity) <<
        ", no-data=" << (outCellStack==peano4::grid::TraversalObserver::NoData) <<
        ", is-in/out=" << peano4::grid::PeanoCurve::isInOutStack(outCellStack)
      );
    }
    return false;
  };
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    storeCell{{logical_type_name}}
  ));
  
"""

    TemplateLeaveCell_FaceStore_MappingCall = """
  // Handle face {{logical_type_name}}
  {
    peano4::datamanagement::FaceMarker  marker(event);

    for (int i=0; i<TwoTimesD; i++) {
      int outFaceStack      = event.getFaceDataTo(i);
      int pick              = event.getFaceDataFrom(i);

      marker.select(pick);

      if (
        outFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchFaceLastTime(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
        {{ACTIVE_ACTION_SET}}.destroyPersistentFace(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        outFaceStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.destroyHangingFace(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        peano4::grid::PeanoCurve::isInOutStack(outFaceStack)
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchFaceLastTime(
           marker
          {{_COMMA_MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS}}
          ,{{MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
    }
  }
"""

    TemplateStoreCell_FaceStore = """
  // Store face {{logical_type_name}}
  std::function<bool ()> storeFace{{logical_type_name}} = 
    [&]()->bool {
    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->popBlock( TwoTimesD );
    for (int i=0; i<TwoTimesD; i++) {
      int inFaceStackPosition  = event.getFaceDataFrom(i);
      int outFaceStack         = event.getFaceDataTo(i);
      logDebug("storeCell(...)", "pos-" << inFaceStackPosition << "->face stack " << outFaceStack );

      peano4::datamanagement::FaceMarker  marker(event,inFaceStackPosition);

      {{full_qualified_type}}& data = view.get(inFaceStackPosition);

      bool dataShallBePushedOntoIntermediateStack = 
        not peano4::grid::PeanoCurve::isInOutStack(outFaceStack)
        and
        outFaceStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        outFaceStack!=peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        outFaceStack!=peano4::grid::TraversalObserver::NoData;
      bool dataShallBePushedOntoOutputStack        = 
        peano4::grid::PeanoCurve::isInOutStack(outFaceStack)
        and
        ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ));
      
      if ( dataShallBePushedOntoIntermediateStack or dataShallBePushedOntoOutputStack ) {
        logDebug( "storeCell(...)", "store face {{name}} with " << marker.x(inFaceStackPosition) << " x " << marker.h() );
        repositories::DataRepository::_{{logical_type_name}}Stack.getForPush( repositories::DataRepository::DataKey(_spacetreeId,outFaceStack))->push(data);
      }
      else {
        logDebug( "storeCell(...)", "do not store face {{name}} with " << marker.x(inFaceStackPosition) << " x " << marker.h() );
      }
    }
    return false;
  };
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    storeFace{{logical_type_name}}
  ));
  
"""

    TemplateLeaveCell_VertexStore_MappingCall = """
  // Handle vertex {{logical_type_name}}
  {
    peano4::datamanagement::VertexMarker  marker(event);

    for (int i=0; i<TwoPowerD; i++) {
      int outVertexStack        = event.getVertexDataTo(i);
      int pick                  = event.getVertexDataFrom(i);

      marker.select(pick);

      if (
        outVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchVertexLastTime(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
        {{ACTIVE_ACTION_SET}}.destroyPersistentVertex(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        outVertexStack==peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.destroyHangingVertex(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
      else if (
        peano4::grid::PeanoCurve::isInOutStack(outVertexStack)
        and
        marker.isLocal()
      ) {
        {{ACTIVE_ACTION_SET}}.touchVertexLastTime(
           marker
          ,{{MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY}}
          {{_COMMA_MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS_COMMA_MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS}}
        );
      }
    }
  }
"""

    TemplateStoreCell_VertexStore = """
  // Store vertex {{logical_type_name}}
  std::function<bool ()> storeVertex{{logical_type_name}} = 
    [&]()->bool {
    auto view = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->popBlock( TwoPowerD );
    for (int i=0; i<TwoPowerD; i++) {
      int inVertexStackPosition  = event.getVertexDataFrom(i);
      int outVertexStack         = event.getVertexDataTo(i);
      logDebug("storeCell(...)", "pos-" << inVertexStackPosition << "->vertex stack " << outVertexStack);

      peano4::datamanagement::VertexMarker  marker(event,inVertexStackPosition);

      {{full_qualified_type}}& data = view.get(inVertexStackPosition);


      bool dataShallBePushedOntoIntermediateStack = 
        not peano4::grid::PeanoCurve::isInOutStack(outVertexStack)
        and
        outVertexStack!=peano4::grid::TraversalObserver::CreateOrDestroyPersistentGridEntity
        and
        outVertexStack!=peano4::grid::TraversalObserver::CreateOrDestroyHangingGridEntity
        and
        outVertexStack!=peano4::grid::TraversalObserver::NoData;
      bool dataShallBePushedOntoOutputStack        = 
        peano4::grid::PeanoCurve::isInOutStack(outVertexStack)
        and
        ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        ));

      if ( dataShallBePushedOntoIntermediateStack or dataShallBePushedOntoOutputStack ) {
        logDebug( "storeCell(...)", "store vertex {{name}} with " << marker.x(inVertexStackPosition) << " x " << marker.h() );
        repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(repositories::DataRepository::DataKey(_spacetreeId,outVertexStack))->push(data);
      }
      else {
        logDebug( "storeCell(...)", "do not store vertex {{name}} with " << marker.x(inVertexStackPosition) << " x " << marker.h() );
      }
    }
    return false;
  };
  
  tasks.push_back( new tarch::multicore::TaskWithoutCopyOfFunctor(
    tarch::multicore::Task::DontFuse,
    tarch::multicore::Task::DefaultPriority,
    storeVertex{{logical_type_name}}
  ));
  
"""

    TemplateLeaveCell_Epilogue = """
  logTraceOutWith1Argument( "leaveCell(...)", _spacetreeId );
}


"""

    TemplateLoadCell_Epilogue = """
  tarch::multicore::spawnAndWait(tasks);
  logTraceOutWith1Argument( "loadCell(...)", _spacetreeId );
}


"""

    TemplateStoreCell_Epilogue = """
  tarch::multicore::spawnAndWait( tasks );
  logTraceOutWith1Argument( "storeCell(...)", _spacetreeId );
}


"""

    def __generate_leaveCell(self, output_file):
        """
        Generates enter cell
        """
        output_file.write(
            jinja2.Template(self.TemplateLeaveCell_Prologue).render(
                **self.mkSubDict(["FULL_QUALIFIED_CLASSNAME"])
            )
        )

        md = self.mkSubDict(
            [
                ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_CELL_EVENT",
                ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS_CELL_EVENT",
            ]
        )
        self.__format_template_per_action(
            output_file, self.TemplateLeaveCell_MappingCall, True, manual_dict=md
        )

        if len(self.cells) > 0:
            md = self.mkSubDict(
                [
                    ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateLeaveCell_CellStore_MappingCall,
                True,
                manual_dict=md,
            )

        if len(self.faces) > 0:
            md = self.mkSubDict(
                [
                    "name",
                    ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                    ",MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateLeaveCell_FaceStore_MappingCall,
                True,
                manual_dict=md,
            )

        if len(self.vertices) > 0:
            md = self.mkSubDict(
                [
                    "name",
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_FINE_GRID_CELL_ARGUMENTS",
                    ",MAPPING_SIGNATURE_COARSE_GRID_VERTICES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_FACES_ARGUMENTS,MAPPING_SIGNATURE_COARSE_GRID_CELL_ARGUMENTS",
                    "MAPPING_SIGNATURE_FINE_GRID_VERTICES_ARGUMENTS_PICK_ENTRY",
                ]
            )
            self.__format_template_per_action(
                output_file,
                self.TemplateLeaveCell_VertexStore_MappingCall,
                True,
                manual_dict=md,
            )

        output_file.write(jinja2.Template(self.TemplateLeaveCell_Epilogue).render({}))

    def __generate_storeCell(self, output_file):
        output_file.write(
            jinja2.Template(self.TemplateStoreCell_Prologue).render(**self.d)
        )

        cells_inverted = [x for x in self.cells]
        faces_inverted = [x for x in self.faces]
        vertices_inverted = [x for x in self.vertices]

        cells_inverted.reverse()
        faces_inverted.reverse()
        vertices_inverted.reverse()

        for cell in cells_inverted:
            temp = {
                "name": cell.name,
                "logical_type_name": cell.get_logical_type_name(),
                "full_qualified_type": cell.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": cell.additional_load_and_store_arguments,
            }
            self.d["name"] = cell.name
            output_file.write(
                jinja2.Template(self.TemplateStoreCell_CellStore).render(**temp)
            )

        for face in faces_inverted:
            temp = {
                "name": face.name,
                "enumeration_type": face.get_enumeration_type(),
                "logical_type_name": face.get_logical_type_name(),
                "full_qualified_type": face.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": face.additional_load_and_store_arguments,
            }
            self.d["name"] = face.name
            output_file.write(
                jinja2.Template(self.TemplateStoreCell_FaceStore).render(**temp)
            )

        for vertex in vertices_inverted:
            temp = {
                "name": vertex.name,
                "enumeration_type": vertex.get_enumeration_type(),
                "logical_type_name": vertex.get_logical_type_name(),
                "full_qualified_type": vertex.get_full_qualified_type(),
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": vertex.additional_load_and_store_arguments,
            }
            self.d["name"] = vertex.name
            output_file.write(
                jinja2.Template(self.TemplateStoreCell_VertexStore).render(**temp)
            )

        output_file.write(jinja2.Template(self.TemplateStoreCell_Epilogue).render({}))

    TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::exchangeAllVerticalDataExchangeStacks( int masterId ) {
  logTraceInWith2Arguments( "exchangeAllVerticalDataExchangeStacks(...)", masterId, _spacetreeId  );
"""

    TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Exchange = """
  peano4::parallel::SpacetreeSet::exchangeAllVerticalDataExchangeStacks(
    {{DATASET}},
    _spacetreeId,
    masterId
  );
"""

    TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Epilogue = """
  logTraceOut( "exchangeAllVerticalDataExchangeStacks(...)" );
}


"""

    TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::exchangeAllHorizontalDataExchangeStacks( bool symmetricDataCardinality ) {
  logTraceInWith2Arguments( "exchangeAllHorizontalDataExchangeStacks(...)", symmetricDataCardinality, _spacetreeId  );
"""

    TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Exchange = """
  peano4::parallel::SpacetreeSet::exchangeAllHorizontalDataExchangeStacks(
    {{DATASET}},
    _spacetreeId,
    symmetricDataCardinality
  );
"""

    TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Epilogue = """
  logTraceOut( "exchangeAllHorizontalDataExchangeStacks(...)" );
}


"""

    TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::exchangeAllPeriodicBoundaryDataStacks() {
  logTraceInWith1Argument( "exchangeAllPeriodicBoundaryDataStacks()", _spacetreeId  );
"""

    TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Exchange = """
  peano4::parallel::SpacetreeSet::exchangeAllPeriodicBoundaryDataStacks(
    {{DATASET}},
    _spacetreeId
  );
"""

    TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Epilogue = """
  logTraceOut( "exchangeAllPeriodicBoundaryDataStacks()" );
}


"""

    TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::streamDataFromSplittingTreeToNewTree(int newWorker) {
  logTraceInWith2Arguments( "streamDataFromSplittingTreeToNewTree(int)", _spacetreeId, newWorker );
"""

    TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Exchange = """
  peano4::parallel::SpacetreeSet::streamDataFromSplittingTreeToNewTree(
    {{DATASET}},
    _spacetreeId,
    newWorker
  );
"""

    TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Epilogue = """
  logTraceOut( "streamDataFromSplittingTreeToNewTree(int)");
}


"""

    TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::streamDataFromJoiningTreeToMasterTree(int master) {
  logTraceInWith2Arguments( "streamDataFromJoiningTreeToMasterTree(int)", _spacetreeId, master );
"""

    TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Exchange = """
  peano4::parallel::SpacetreeSet::streamDataFromJoiningTreeToMasterTree(
    {{DATASET}},
    _spacetreeId,
    master
  );
"""

    TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Epilogue = """
  logTraceOut( "streamDataFromJoiningTreeToMasterTree(int)");
}


"""

    TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::finishAllOutstandingSendsAndReceives() {
  logTraceInWith1Argument( "finishAllOutstandingSendsAndReceives()", _spacetreeId );
"""

    TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Exchange = """
  peano4::parallel::SpacetreeSet::finishAllOutstandingSendsAndReceives(
    {{DATASET}},
    _spacetreeId
  );
"""

    TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Epilogue = """
  logTraceOut( "finishAllOutstandingSendsAndReceives()");
}


"""

    TemplateSendVertex_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::sendVertex(int position, int toStack, ::peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent&  event) {
  logTraceInWith4Arguments( "sendVertex(int,int,int)", position, toStack, event.toString(), _spacetreeId );

"""

    TemplateSendVertex_SendCall = """
  {
    peano4::datamanagement::VertexMarker  marker(event,position);

    const {{full_qualified_type}}& data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(
      repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack)
    )->top(TwoPowerD-1-position);
    if (
      ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      (
        data.send(marker)
        or
        context == TraversalObserver::SendReceiveContext::ForkDomain
        or
        context == TraversalObserver::SendReceiveContext::JoinDomain
      )
    ) {
      logDebug( "sendVertex(...)", "send out " << data.toString() << " to stack " << toStack << " on tree " << _spacetreeId << " for marker " << marker.toString() );

      repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
        _spacetreeId, toStack
      ) -> push(data);
    }
  }
"""

    TemplateSendVertex_Epilogue = """
  logTraceOut( "sendVertex(int,int,int)");
}


"""

    TemplateSendFace_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::sendFace(int position, int toStack, ::peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent&  event) {
  logTraceInWith4Arguments( "sendFace(int,int,int)", position, toStack, event.toString(), _spacetreeId );

"""

    TemplateSendFace_SendCall = """
  {
    peano4::datamanagement::FaceMarker  marker(event,position);
    const {{full_qualified_type}}& data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(
      repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack)
    )->top(TwoTimesD-1-position);
    if (
      ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      (data.send(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        )
        or
        context == TraversalObserver::SendReceiveContext::ForkDomain
        or
        context == TraversalObserver::SendReceiveContext::JoinDomain
      )
    ) {
      logDebug( "sendFace(...)", "send out " << data.toString() << " to stack " << toStack << " on tree " << _spacetreeId << " for marker " << marker.toString() );

      repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
        _spacetreeId, toStack
      ) -> push(data);
    }
  }
"""

    TemplateSendFace_Epilogue = """
  logTraceOut( "sendFace(int,int,int)");
}


"""

    TemplateSendCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::sendCell(int toStack, ::peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent&  event) {
  logTraceInWith3Arguments( "sendCell(int,int,int)", toStack, event.toString(), _spacetreeId );

"""

    TemplateSendCell_SendCall = """
  {
    peano4::datamanagement::CellMarker  marker(event);
    const {{full_qualified_type}}&  data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPop( repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack))->top();
    if (
      ::peano4::grid::storePersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      (
        data.send(
          marker
          {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
        )
        or
        context == TraversalObserver::SendReceiveContext::ForkDomain
        or
        context == TraversalObserver::SendReceiveContext::JoinDomain
      )
    ) {
      logDebug( "sendCell(...)", "send out " << data.toString() << " to stack " << toStack << " on tree " << _spacetreeId );

      repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
        _spacetreeId, toStack
      ) -> push(data);
    }
  }
"""

    TemplateSendCell_Epilogue = """
  logTraceOut( "sendCell(int,int,int)");
}


"""

    TemplateReceiveAndMergeCell_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::receiveAndMergeCell(int fromStack, peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent& event ) {
  logTraceInWith3Arguments( "receiveAndMergeCell(...)", fromStack, event.toString(), _spacetreeId );
"""

    TemplateReceiveAndMergeCell_Epilogue = """
  logTraceOut( "receiveAndMergeCell(...)");
}


"""

    TemplateReceiveAndMergeVertex_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::receiveAndMergeVertex(int position, int fromStack, peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent& event) {
  logTraceInWith4Arguments( "receiveAndMergeVertex(...)", position, fromStack, event.toString(), _spacetreeId );
"""

    TemplateReceiveAndMergeVertex_Epilogue = """
  logTraceOut( "receiveAndMergeVertex(...)");
}


"""

    TemplateReceiveAndMergeFace_Epilogue = """
  logTraceOut( "receiveAndMergeFace(...)");
}


"""

    TemplateReceiveAndMergeFace_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::receiveAndMergeFace(int position, int fromStack, peano4::grid::TraversalObserver::SendReceiveContext context, const peano4::grid::GridTraversalEvent& event) {
  logTraceInWith4Arguments( "receiveAndMergeFace(...)", position, fromStack, event.toString(), _spacetreeId );
"""

    TemplateReceiveAndMergeVertex_ReceiveAndMergeCalls = """
  {
    peano4::datamanagement::VertexMarker  marker(event,position);
    {{full_qualified_type}}& data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
      repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack)
    )->top(TwoPowerD-1-position);

    if (
      ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      data.receiveAndMerge(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      )
    ) {
      assertion5(
        not repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(_spacetreeId, fromStack)->empty(),
        _spacetreeId, fromStack, event.toString(), position, marker.toString()
      );

      auto   incomingData = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(
        _spacetreeId, fromStack
      )->pop() );

      logDebug( "receiveAndMergeVertex(...)", "merge " << incomingData.toString() << " into " << data.toString() );

      if (context==::peano4::grid::TraversalObserver::SendReceiveContext::PeriodicBoundaryDataSwap) {
        // @todo Different to faces. As we can have diagonal exchange, too,
        //       we only know that it has to be unequal
        assertion8(
          data.getDebugX()!=incomingData.getDebugX(),
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(),
          fromStack, event.toString(), marker.toString(), _spacetreeId );
        assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }
      else {
        assertionVectorNumericalEquals6(
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(), fromStack, event.toString(), marker.toString(), _spacetreeId );
          assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }

      data.merge(context, incomingData, marker, _spacetreeId);
    }
  }
"""

    TemplateReceiveAndMergeFace_ReceiveAndMergeCalls = """
  {
    peano4::datamanagement::FaceMarker  marker(event,position);
    {{full_qualified_type}}& data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
      repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack)
    )->top(TwoTimesD-1-position);

    if (
      ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      data.receiveAndMerge(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      )
    ) {
      auto   incomingData = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(
        _spacetreeId, fromStack
      )->pop() );

      logDebug( "receiveAndMergeFace(...)", "merge " << incomingData.toString() << " into " << data.toString() << " within marker " << marker.toString() );

      if (context==::peano4::grid::TraversalObserver::SendReceiveContext::PeriodicBoundaryDataSwap) {
        assertion8(
          tarch::la::countEqualEntries(data.getDebugX(), incomingData.getDebugX())==Dimensions-1,
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(),
          fromStack, event.toString(), marker.toString(), _spacetreeId );
        assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }
      else {
        assertionVectorNumericalEquals6(
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(), fromStack, event.toString(), marker.toString(), _spacetreeId );
          assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }

      data.merge(context,incomingData, marker, _spacetreeId);
    }
  }
"""

    TemplateReceiveAndMergeCell_ReceiveAndMergeCalls = """
  {
    peano4::datamanagement::CellMarker  marker(event);

    {{full_qualified_type}}& data = repositories::DataRepository::_{{logical_type_name}}Stack.getForPush(
      repositories::DataRepository::DataKey(_spacetreeId,peano4::grid::PeanoCurve::CallStack)
    )->top();

    if (
      ::peano4::grid::loadPersistently( {{full_qualified_type}}::loadStoreComputeFlag(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      ))
      and
      data.receiveAndMerge(
        marker
        {% for arg in ADDITIONAL_LOAD_STORE_ARGUMENTS %}, {{arg[0]}} {% endfor %}
      )
    ) {
      auto   incomingData = std::move( repositories::DataRepository::_{{logical_type_name}}Stack.getForPop(
        _spacetreeId, fromStack
      )->pop() );

      logDebug( "receiveAndMergeCell(...)", "merge " << incomingData.toString() << " into " << data.toString() );

      if (context==::peano4::grid::TraversalObserver::SendReceiveContext::PeriodicBoundaryDataSwap) {
        assertion8(
          tarch::la::countEqualEntries(data.getDebugX(), incomingData.getDebugX())==Dimensions-1,
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(),
          fromStack, event.toString(), marker.toString(), _spacetreeId );
        assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }
      else {
        assertionVectorNumericalEquals6(
          data.getDebugX(), incomingData.getDebugX(),
          data.getDebugH(), incomingData.getDebugH(), fromStack, event.toString(), marker.toString(), _spacetreeId );
          assertionVectorNumericalEquals6(
          data.getDebugH(), incomingData.getDebugH(),
          data.getDebugX(), incomingData.getDebugX(), fromStack, event.toString(), marker.toString(), _spacetreeId );
      }

      data.merge(context,incomingData, marker, _spacetreeId);
    }
  }
"""

    TemplateExchangeRoutines_deleteAllStacks_Prologue = """
void {{FULL_QUALIFIED_CLASSNAME}}::deleteAllStacks() {
  logTraceInWith1Argument( "deleteAllStacks()", _spacetreeId );
"""

    TemplateExchangeRoutines_deleteAllStacks_Exchange = """
  peano4::parallel::SpacetreeSet::deleteAllStacks(
    {{DATASET}},
    _spacetreeId
  );
"""

    TemplateExchangeRoutines_deleteAllStacks_Epilogue = """
  logTraceOut( "deleteAllStacks()");
}


"""

    def __generate_exchange_routines(self, output_file):
        s = ""

        generic_dict_for_prologue_and_epilogue = {
            "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"]
        }

        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllVerticalDataExchangeStacks_Epilogue
        ).render({})

        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllHorizontalDataExchangeStacks_Epilogue
        ).render({})

        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_exchangeAllPeriodicBoundaryDataStacks_Epilogue
        ).render({})

        s += jinja2.Template(
            self.TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_streamDataFromSplittingTreeToNewTree_Epilogue
        ).render({})

        s += jinja2.Template(
            self.TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_streamDataFromJoiningTreeToMasterTree_Epilogue
        ).render({})

        s += jinja2.Template(
            self.TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_finishAllOutstandingSendsAndReceives_Epilogue
        ).render({})

        s += jinja2.Template(self.TemplateSendVertex_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for vertex in self.vertices:
            temp = {
                "name": vertex.name,
                "enumeration_type": vertex.get_enumeration_type(),
                "logical_type_name": vertex.get_logical_type_name(),
                "full_qualified_type": vertex.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": vertex.additional_load_and_store_arguments,
            }
            s += jinja2.Template(self.TemplateSendVertex_SendCall).render(**temp)
        s += jinja2.Template(self.TemplateSendVertex_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(self.TemplateSendFace_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for face in self.faces:
            temp = {
                "name": face.name,
                "enumeration_type": face.get_enumeration_type(),
                "logical_type_name": face.get_logical_type_name(),
                "full_qualified_type": face.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": face.additional_load_and_store_arguments,
            }
            s += jinja2.Template(self.TemplateSendFace_SendCall).render(**temp)
        s += jinja2.Template(self.TemplateSendFace_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(self.TemplateSendCell_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for cell in self.cells:
            temp = {
                "name": cell.name,
                "logical_type_name": cell.get_logical_type_name(),
                "full_qualified_type": cell.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": cell.additional_load_and_store_arguments,
            }
            s += jinja2.Template(self.TemplateSendCell_SendCall).render(**temp)
        s += jinja2.Template(self.TemplateSendCell_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(self.TemplateReceiveAndMergeCell_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for cell in self.cells:
            temp = {
                "name": cell.name,
                "logical_type_name": cell.get_logical_type_name(),
                "full_qualified_type": cell.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": cell.additional_load_and_store_arguments,
            }
            s += jinja2.Template(
                self.TemplateReceiveAndMergeCell_ReceiveAndMergeCalls
            ).render(**temp)
        s += jinja2.Template(self.TemplateReceiveAndMergeCell_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(self.TemplateReceiveAndMergeFace_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for face in self.faces:
            temp = {
                "name": face.name,
                "enumeration_type": face.get_enumeration_type(),
                "logical_type_name": face.get_logical_type_name(),
                "full_qualified_type": face.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": face.additional_load_and_store_arguments,
            }
            s += jinja2.Template(
                self.TemplateReceiveAndMergeFace_ReceiveAndMergeCalls
            ).render(**temp)
        s += jinja2.Template(self.TemplateReceiveAndMergeFace_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(self.TemplateReceiveAndMergeVertex_Prologue).render(
            **generic_dict_for_prologue_and_epilogue
        )
        for vertex in self.vertices:
            temp = {
                "name": vertex.name,
                "enumeration_type": vertex.get_enumeration_type(),
                "logical_type_name": vertex.get_logical_type_name(),
                "full_qualified_type": vertex.get_full_qualified_type(),
                "FULL_QUALIFIED_CLASSNAME": self.d["FULL_QUALIFIED_CLASSNAME"],
                "ADDITIONAL_LOAD_STORE_ARGUMENTS": vertex.additional_load_and_store_arguments,
            }
            s += jinja2.Template(
                self.TemplateReceiveAndMergeVertex_ReceiveAndMergeCalls
            ).render(**temp)
        s += jinja2.Template(self.TemplateReceiveAndMergeVertex_Epilogue).render(
            **generic_dict_for_prologue_and_epilogue
        )

        s += jinja2.Template(
            self.TemplateExchangeRoutines_deleteAllStacks_Prologue
        ).render(**generic_dict_for_prologue_and_epilogue)
        for cell in self.cells:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_deleteAllStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + cell.get_logical_type_name()
                    + "Stack"
                }
            )
        for face in self.faces:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_deleteAllStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + face.get_logical_type_name()
                    + "Stack"
                }
            )
        for vertex in self.vertices:
            s += jinja2.Template(
                self.TemplateExchangeRoutines_deleteAllStacks_Exchange
            ).render(
                **{
                    "DATASET": "repositories::DataRepository::_"
                    + vertex.get_logical_type_name()
                    + "Stack"
                }
            )
        s += jinja2.Template(
            self.TemplateExchangeRoutines_deleteAllStacks_Epilogue
        ).render(**generic_dict_for_prologue_and_epilogue)

        output_file.write(s)

    TemplateImplementationFilePrologue = """
#include "{{CLASSNAME}}.h"
#include "repositories/DataRepository.h"

#include "peano4/grid/PeanoCurve.h"

#include "peano4/datamanagement/VertexEnumerator.h"
#include "peano4/datamanagement/VertexMarker.h"
#include "peano4/datamanagement/FaceEnumerator.h"
#include "peano4/datamanagement/FaceMarker.h"
#include "peano4/datamanagement/CellMarker.h"

#include "peano4/parallel/SpacetreeSet.h"


tarch::logging::Log {{FULL_QUALIFIED_CLASSNAME}}::_log( "{{FULL_QUALIFIED_CLASSNAME}}" );

"""

    def __generate_implementation(self, overwrite, full_qualified_filename):
        if write_file(overwrite, self.default_overwrite, full_qualified_filename):
            import inspect, os

            cuda = using_cuda_backend()
            if cuda:
                full_qualified_filename = full_qualified_filename.replace(".cpp", ".cu")

            print(
                "{} written by {}".format(
                    full_qualified_filename,
                    os.path.basename(inspect.getfile(self.__class__)),
                )
            )

            output_file = open(full_qualified_filename, "w")
            output_file.write(
                jinja2.Template(self.TemplateImplementationFilePrologue).render(
                    **self.d
                )
            )

            self.__generate_constructor(output_file)
            self.__generate_clone(output_file)
            self.__generate_getGridControlEvents(output_file)
            self.__generate_beginTraversal(output_file)
            self.__generate_endTraversal(output_file)
            self.__generate_prepareTraversal(output_file)
            self.__generate_unprepareTraversal(output_file)
            self.__generate_loadCell(output_file)
            self.__generate_enterCell(output_file)
            self.__generate_leaveCell(output_file)
            self.__generate_storeCell(output_file)
            self.__generate_exchange_routines(output_file)

    #
    # @todo Sollte man mit Jinja 2 nicht mehr brauchen
    #
    def get_cpp_file_name(self):
        return self.subdirectory + "/" + self.classname + ".cpp"

    def generate(self, overwrite, directory):
        if not os.path.exists(directory + "/" + self.subdirectory):
            os.mkdir(directory + "/" + self.subdirectory)

        cpp_filename = directory + "/" + self.get_cpp_file_name()

        self.__generate_header(overwrite, directory)
        self.__generate_implementation(overwrite, cpp_filename)
