# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import os
import re
import jinja2

from copy import deepcopy

from .Helper import write_file
from .Helper import using_cuda_backend

from .Overwrite import Overwrite

class Jinja2TemplatedHeaderImplementationFilePair(object):
  """!
Represents a pair of C++ header and implementation files

The content of these files is represented by jinja2 templates. While
sorting out these templates, the class defines a couple of default
dictionary entries besides the ones manually handed over by the user.
  """

  def __init__(self,
               headfile_template,
               cppfile_template,
               classname,
               namespace,
               subdirectory,
               dictionary,
               default_overwrite=True,
               apply_iteratively=False):
    """
      The template files should be fully qualified
      classname is a string
      namespace is a (possibly empty) list of strings

      cppfile_template can be None
    """
    self.headerfile_template  = headfile_template
    self.cppfile_template     = cppfile_template
    self.classname            = classname
    self.namespace            = namespace
    self.subdirectory         = subdirectory
    self.default_overwrite    = default_overwrite
    self.apply_iteratively    = apply_iteratively

    self.d = deepcopy(dictionary)
    self.d[ "CLASSNAME" ] = classname
    self.d[ "NAMESPACE" ] = []
    self.d[ "FULL_QUALIFIED_NAMESPACE" ] = ""

    for i in namespace:
      self.d["NAMESPACE"].append(i)
      self.d[ "FULL_QUALIFIED_NAMESPACE" ] += "::" + i

  def __generate_file(self,overwrite,full_qualified_filename,template_file):
    if template_file!=None and write_file(overwrite,self.default_overwrite,full_qualified_filename):
      import inspect
      cuda = using_cuda_backend()
      if cuda:
        full_qualified_filename = full_qualified_filename.replace(".cpp", ".cu")

      print( "{} written by {} (from template {})".format(full_qualified_filename, os.path.basename(inspect.getfile(self.__class__)),  template_file))

      template_loader = jinja2.FileSystemLoader(searchpath=os.path.split(template_file)[0])
      templateEnv = jinja2.Environment(loader=template_loader, undefined=jinja2.DebugUndefined)
      template = templateEnv.get_template( os.path.split(template_file)[1] )

      if self.apply_iteratively:
        old_size                = 0
        rendered_text           = template.render(self.d)
        while "{{" in rendered_text and len(rendered_text)!=old_size:
          old_size      = len(rendered_text)
          template      = jinja2.Template(rendered_text, undefined=jinja2.DebugUndefined)
          rendered_text = template.render( self.d )

      with open( full_qualified_filename, "w" ) as output:
        output.write( template.render(self.d) )

  def generate(self,overwrite,directory):
    if not os.path.exists( directory + "/" + self.subdirectory ):
      os.mkdir(directory + "/" + self.subdirectory)

    header_filename = directory + "/" + self.subdirectory + "/" + self.classname + ".h"

    cuda = using_cuda_backend()
    if not cuda:
      cpp_filename    = directory + "/" + self.subdirectory + "/" + self.classname + ".cpp"
    else:
      cpp_filename    = directory + "/" + self.subdirectory + "/" + self.classname + ".cu"

    self.__generate_file(overwrite,header_filename,self.headerfile_template)
    self.__generate_file(overwrite,cpp_filename   ,self.cppfile_template)
