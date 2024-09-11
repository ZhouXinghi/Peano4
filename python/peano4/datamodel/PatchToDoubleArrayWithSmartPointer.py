# This file is part of the Peano project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import peano4.output.Jinja2TemplatedHeaderImplementationFilePair
import os

from .DoF import DoFAssociation
from .PatchToDoubleArray import FloatTypes
from .PatchToDoubleArray import PatchToDoubleArray


class PatchToDoubleArrayWithSmartPointer(PatchToDoubleArray):
    """!

    Realise patch via smart pointers

    Converter which maps the patch 1:1 onto a double array administered
    through a smart pointer. This is an alternative to the patch realisation
    on the call stack, which can become really costly and lead to call stack
    overflows.

    """

    def __init__(self, patch, float_type="double"):
        """ """
        super(PatchToDoubleArrayWithSmartPointer, self).__init__(patch, float_type)

    def get_stack_container(self):
        return (
            "peano4::stacks::STDVectorStackOverSmartPointers< "
            + self.data.get_full_qualified_type()
            + " >"
        )

    def construct_output(self, output):
        d = self._get_dictionary_for_output()
        """
      Pass in a version of output
    """
        output.makefile.add_cpp_file(
            self.data.namespace[-1] + "/" + self.data.name + ".cpp", generated=True
        )
        templatefile_prefix = (
            os.path.realpath(__file__).replace(".pyc", "").replace(".py", "")
        )
        generated_files = peano4.output.Jinja2TemplatedHeaderImplementationFilePair(
            templatefile_prefix + ".template.h",
            templatefile_prefix + ".template.cpp",
            self.data.name,
            self.data.namespace,
            self.data.namespace[-1],
            d,
            True,
        )
        output.add(generated_files)

    def get_header_file_include(self):
        return (
            """
#include "peano4/stacks/STDVectorStackOverSmartPointers.h"
#include \""""
            + self.data.namespace[-1]
            + """/"""
            + self.data.name
            + """.h"
"""
        )
