# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

from SWE import SWE

if __name__ == "__main__":
    swe = SWE()

    my_parser = swe.setup_parser()
    my_parser.set_defaults(
        solver="RusanovRiemannGlobalAdaptiveFV",
        end_time=250.0,
        width="[10000, 10000]",  # [m]
        min_depth=2,
        amr_levels=6,
        patch_size=4,
        peano_dir="../../../../",
    )
    my_args = my_parser.parse_args()

    my_solver = swe.setup_solver(my_args)
    my_solver.add_user_solver_includes(
        """
#include "../SWE.h"
#include "../FWave.h"
#include "../HLLEM.h"
#include "../Rusanov.h"
#include "../TopologyParser.h"
"""
    )
    my_solver.add_solver_constants(
        """
TopologyParser topologyParser = TopologyParser(
  \"artificial_bath_1000.nc\", \
  \"artificial_displ_1000.nc\");
"""
    )

    my_project = swe.setup_project(my_args, my_solver)

    my_project.constants.export_constexpr_with_type(
        "DRY_TOLERANCE", str(1e-2), "double"
    )

    my_project.output.makefile.add_cpp_file("SWE.cpp")
    my_project.output.makefile.add_cpp_file("../SWE.cpp")
    my_project.output.makefile.add_cpp_file("../FWave.cpp")
    my_project.output.makefile.add_cpp_file("../HLLEM.cpp")
    my_project.output.makefile.add_cpp_file("../Rusanov.cpp")
    my_project.output.makefile.add_cpp_file("../NetCDFReader.cpp")
    my_project.output.makefile.add_cpp_file("../TopologyParser.cpp")

    swe.build(my_args, my_project)

    print(my_args)
    print(my_solver)
    print(my_project)
