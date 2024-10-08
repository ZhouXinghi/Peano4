##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# Generate the Kernels.h header
#


from .abstractModelBaseClass import AbstractModelBaseClass


class KernelsHeaderModel(AbstractModelBaseClass):

    def generateCode(self):
        
        print(self.context["solverName"])

        self.context["solverNamespace"] = self.context["solverName"].split("::")[0]
        self.context["solverNamespaceList"] = self.context["solverName"].split("::")
        self.context["solverNamespaceList"].pop(-1)
        self.context["solverClass"] = self.context["solverName"].split("::")[-1]

        self.render((self.context["kernelType"], "Kernels_h.template"), "Kernels.h")
