# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .DaStGenToLegacyTool import DaStGenToLegacyTool
from .DoF                 import DoF



class DaStGen(DoF):
  def __init__(self, name):
    """
      dim should be a two-tuple or a three-tuple of integers, so a construction 
      could look like
      
        cell_data = peano4.datamodel.Patch( (6,6,6), "Fluid" )

      attributes_double  Series of strings which encode doubles.
      attributes_enum    Series of tuples with a string plus an array of strings.
    """
    super(DaStGen, self).__init__(name)
    self.generator = DaStGenToLegacyTool(self)
    self.attributes_double  = []
    self.attributes_integer = []
    self.attributes_enum    = []


  def add_double_scalar(self, name):
    self.attributes_double.append( name )


  def add_integer_scalar(self, name):
    self.attributes_integer.append( name )

    
  def add_enum(self, enum_name, enum_variants ):
    """
      enum  Is a set or series of strings.
    """
    self.attributes_enum.append( (enum_name, enum_variants) )
    #  attributes_enum    Series of tuples with a string plus an array of strings.
    
