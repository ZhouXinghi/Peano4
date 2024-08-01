# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .DoF                import DoF
from .DynamicArrayOverPrimitivesToStdVector import DynamicArrayOverPrimitivesToStdVector


class DynamicArrayOverPrimitives(DoF):
  def __init__(self, name, primitive_type="double"):
    """
     type should refer to a C++ primitive (which is also supported by MPI).

    """
    super(DynamicArrayOverPrimitives, self).__init__(name)
    self.primitive_type = primitive_type
    self.generator      = DynamicArrayOverPrimitivesToStdVector(self)




