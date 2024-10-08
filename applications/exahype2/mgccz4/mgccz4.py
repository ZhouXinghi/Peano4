import os
import argparse

import peano4
import exahype2


modes = { 
  "release": peano4.output.CompileMode.Release,
  "trace":   peano4.output.CompileMode.Trace,
  "assert":  peano4.output.CompileMode.Asserts, "stats":  peano4.output.CompileMode.Stats,
  "debug":   peano4.output.CompileMode.Debug,
}

floatparams = {
        "GLMc0":1.5, "GLMc":1.2, "GLMd":2.0, "GLMepsA":1.0, "GLMepsP":1.0,
        "GLMepsD":1.0, "itau":1.0, "k1":0.0, "k2":0.0, "k3":0.0, "eta":0.0,
        "f":0.0, "g":0.0, "xi":0.0, "e":1.0, "c":1.0, "mu":0.2, "ds":1.0,
        "sk":0.0, "bs":0.0}
intparams = {"LapseType":0}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ExaHyPE 2 - MGCCZ4 script')
    parser.add_argument("-cs",   "--cell-size",       dest="max_h",           type=float, default=0.4,  help="Mesh size" )
    parser.add_argument("-ps",   "--patch-size",      dest="patch_size",      type=int, default=6,    help="Patch size, i.e. number of volumes per cell per direction" )
    parser.add_argument("-plt",  "--plot-step-size",  dest="plot_step_size",  type=float, default=0.04, help="Plot step size (0 to switch it off)" )
    parser.add_argument("-m",    "--mode",            dest="mode",            default="release",  help="|".join(modes.keys()) )
    parser.add_argument("-ext",  "--extension",       dest="extension",       choices=["none", "gradient", "adm"],   default="none",  help="Pick extension, i.e. what should be plotted on top. Default is none" )
    parser.add_argument("-impl", "--implementation",  dest="implementation",  choices=["ader-fixed", "fv-fixed", "fv-fixed-enclave", "fv-adaptive" ,"fv-adaptive-enclave", "fv-optimistic-enclave", "fv-fixed-gpu"], required="True",  help="Pick solver type" )
    parser.add_argument("-no-pbc",  "--no-periodic-boundary-conditions",      dest="periodic_bc", action="store_false", default="True",  help="switch on or off the periodic BC" )
    parser.add_argument("-et",   "--end-time",        dest="end_time",        type=float, default=1.0, help="End (terminal) time" )


    for k, v in floatparams.items(): parser.add_argument("--{}".format(k), dest="MGCCZ4{}".format(k), type=float, default=v, help="default: %(default)s")
    for k, v in intparams.items():   parser.add_argument("--{}".format(k), dest="MGCCZ4{}".format(k), type=int, default=v, help="default: %(default)s")

    args = parser.parse_args()

    SuperClass = None
    
    if args.implementation=="fv-fixed":
       SuperClass = exahype2.solvers.fv.GenericRusanovFixedTimeStepSize
    if args.implementation=="fv-fixed-enclave":
       SuperClass = exahype2.solvers.fv.GenericRusanovFixedTimeStepSizeWithEnclaves
    if args.implementation=="fv-adaptive":
       SuperClass = exahype2.solvers.fv.GenericRusanovAdaptiveTimeStepSize
    if args.implementation=="fv-adaptive-enclave":
       SuperClass = exahype2.solvers.fv.GenericRusanovAdaptiveTimeStepSizeWithEnclaves
    if args.implementation=="fv-optimistic-enclave":
       SuperClass = exahype2.solvers.fv.GenericRusanovOptimisticTimeStepSizeWithEnclaves
    if args.implementation=="ader-fixed":
       SuperClass = exahype2.solvers.aderdg.NonFusedGenericRusanovFixedTimeStepSize
    if args.implementation=="fv-fixed-gpu":
       SuperClass = exahype2.solvers.fv.GenericRusanovFixedTimeStepSizeWithAccelerator

    class MGCCZ4Solver( SuperClass ):
      def __init__(self, name, patch_size, min_h, max_h ):
        unknowns = {
          "G":6,
          "K":6,
          "theta":1,
          "Z":3,
          "lapse":1,
          "shift":3,
          "b":3,
          "dLapse":3,
          "dxShift":3,
          "dyShift":3,
          "dzShift":3,
          "dxG":6,
          "dyG":6,
          "dzG":6,
          "traceK":1,
          "phi":1,
          "P":3,
          "K0":1,
          "MGphi":1,
          "Kphi":1,
          "Pi":3
        }

        number_of_unknowns = sum(unknowns.values())

        if SuperClass==exahype2.solvers.fv.GenericRusanovFixedTimeStepSize or \
           SuperClass==exahype2.solvers.fv.GenericRusanovFixedTimeStepSizeWithAccelerator or \
           SuperClass==exahype2.solvers.fv.GenericRusanovFixedTimeStepSizeWithEnclaves:
          SuperClass.__init__(
            self,
            name=name, patch_size=patch_size,
            unknowns=number_of_unknowns,
            auxiliary_variables=0,
            min_h=min_h, max_h=max_h,
            time_step_size=1e-2
          )
        else:
          SuperClass.__init__(
            self,
            name=name, patch_size=patch_size,
            unknowns=number_of_unknowns,
            auxiliary_variables=0,
            min_h=min_h, max_h=max_h,
            time_step_relaxation=0.1
          )

        self._solver_template_file_class_name = SuperClass.__name__

        #pde = exahype2.sympy.PDE(unknowns=self._unknowns,auxiliary_variables=self._auxiliary_variables,dimensions=3)

        self.set_implementation(
          boundary_conditions=exahype2.solvers.PDETerms.User_Defined_Implementation,
          flux=exahype2.solvers.PDETerms.None_Implementation,
          ncp=exahype2.solvers.PDETerms.User_Defined_Implementation,
          source_term=exahype2.solvers.PDETerms.User_Defined_Implementation
        )

      def get_user_action_set_includes(self):
        """
         We take this routine to add a few additional include statements.
        """
        return SuperClass.get_user_action_set_includes(self) + """
    #include "../MGCCZ4Kernels.h"
    #include "exahype2/PatchUtils.h"
    """


      def add_constraint_verification(self):
        """

         Add the constraint verification code

         We introduce new auxiliary variables. Prior to each time step, I
         compute the Laplacian and store it in the auxiliary variable. This
         is kind of a material parameter F(Q) which does not feed back into
         the solution.

         Changing the number of unknowns a posteriori is a delicate update
         to the solver, so we invoke the constructor again to be on the safe
         side.

        """
        self._auxiliary_variables = 6

        self.set_preprocess_reconstructed_patch_kernel( """
        const int patchSize = """ + str( self._patch.dim[0] ) + """;
        double volumeH = ::exahype2::getVolumeLength(marker.h(),patchSize);
        const int n_a_v=6;
        dfor(cell,patchSize) {
          tarch::la::Vector<Dimensions,int> currentCell = cell + tarch::la::Vector<Dimensions,int>(1);

          // This constraint routine will evaluate both the solution per voxel
          // plus the derivative. The latter is something that we don't have, i.e.
          // we have to reconstruct it manually.

          // See the docu in PDE.h
          double gradQ[3*64]={ 0 };

          // Lets look left vs right and compute the gradient. Then, lets
          // loop up and down. So we look three times for the respective
          // directional gradients
          for (int d=0; d<3; d++) {
            tarch::la::Vector<Dimensions,int> leftCell  = currentCell;
            tarch::la::Vector<Dimensions,int> rightCell = currentCell;
            leftCell(d)  -= 1;
            rightCell(d) += 1;
            const int leftCellSerialised  = peano4::utils::dLinearised(leftCell, patchSize + 2*1);
            const int rightCellSerialised = peano4::utils::dLinearised(rightCell,patchSize + 2*1);
            for(int i=0; i<64; i++) {
              gradQ[d*64+i] = ( oldQWithHalo[rightCellSerialised*(64+n_a_v)+i] - oldQWithHalo[leftCellSerialised*(64+n_a_v)+i] ) / 2.0 / volumeH;
            }
          }

          // We will use a Fortran routine to compute the constraints per
          // Finite Volume
          double constraints[n_a_v]={ 0 };

          // Central cell
          const int cellSerialised  = peano4::utils::dLinearised(currentCell, patchSize + 2*1);

          admconstraints(constraints,oldQWithHalo+cellSerialised*(64+n_a_v),gradQ);

          for(int i=0;i<n_a_v;i++){
            oldQWithHalo[cellSerialised*(64+n_a_v)+64+i] = constraints[i];
          }
        }
    """)

        self.create_data_structures()
        self.create_action_sets()


      def add_derivative_calculation(self):
        """

         Add the constraint verification code

         We introduce new auxiliary variables. Prior to each time step, I
         compute the Laplacian and store it in the auxiliary variable. This
         is kind of a material parameter F(Q) which does not feed back into
         the solution.

         Changing the number of unknowns a posteriori is a delicate update
         to the solver, so we invoke the constructor again to be on the safe
         side.

        """
        self._auxiliary_variables = 64*3

        self.set_preprocess_reconstructed_patch_kernel( """
        const int patchSize = """ + str( self._patch.dim[0] ) + """;
        double volumeH = ::exahype2::getVolumeLength(marker.h(),patchSize);
        dfor(cell,patchSize) {
          tarch::la::Vector<Dimensions,int> currentCell = cell + tarch::la::Vector<Dimensions,int>(1);
          const int cellSerialised  = peano4::utils::dLinearised(currentCell, patchSize + 2*1);

          // Lets look left vs right and compute the gradient. Then, lets
          // loop up and down. So we look three times for the respective
          // directional gradients
          for (int d=0; d<3; d++) {
            tarch::la::Vector<Dimensions,int> leftCell  = currentCell;
            tarch::la::Vector<Dimensions,int> rightCell = currentCell;
            leftCell(d)  -= 1;
            rightCell(d) += 1;
            const int leftCellSerialised  = peano4::utils::dLinearised(leftCell, patchSize + 2*1);
            const int rightCellSerialised = peano4::utils::dLinearised(rightCell,patchSize + 2*1);
            for (int i=0; i<64; i++) {
              oldQWithHalo[cellSerialised*(64*4)+64+i*3+d] =
                ( oldQWithHalo[rightCellSerialised*(64*4)+i] - oldQWithHalo[leftCellSerialised*(64*4)+i] ) / 2.0 / volumeH;
            }
          }
        }
    """)

        self.create_data_structures()
        self.create_action_sets()



    userwarnings = []

    project = exahype2.Project( ["examples", "exahype2", "mgccz4"], "mgccz4" )

    is_aderdg = False
    solver_name = "MGCCZ4"
    try:
      if SuperClass==exahype2.solvers.aderdg.NonFusedGenericRusanovFixedTimeStepSize:
        is_aderdg = True
        order = 3
        unknowns = 64
        time_step_size = 0.001
    except Exception as e:
        msg = "Warning: ADER-DG no supported on this machine"
        print(msg)
        userwarnings.append((msg,e))

    if is_aderdg:
      solver_name    = "ADERDG" + solver_name
    else:
      solver_name    = "FiniteVolume" + solver_name

    if SuperClass == exahype2.solvers.fv.GenericRusanovFixedTimeStepSizeWithAccelerator:
      solver_name += "OnGPU"

    if is_aderdg:
      my_solver = exahype2.solvers.aderdg.NonFusedGenericRusanovFixedTimeStepSize(
          solver_name, order, unknowns, 0, #auxiliary_variables
          exahype2.solvers.aderdg.Polynomials.Gauss_Legendre,
          args.max_h, args.max_h, time_step_size,
          flux = None,
          ncp  = exahype2.solvers.PDETerms.User_Defined_Implementation,
          sources = exahype2.solvers.PDETerms.User_Defined_Implementation
      )
    else:
      my_solver = MGCCZ4Solver(solver_name, args.patch_size, args.max_h, args.max_h)

      if args.extension=="gradient":
        my_solver.add_derivative_calculation()
      if args.extension=="adm":
        my_solver.add_constraint_verification()

    myscenario = 2 # 0...gaugewave-c++  1...linearwave 2...forcedflat

    solverconstants=""
    for k, v in floatparams.items(): solverconstants+= "static constexpr double {} = {};\n".format("MGCCZ4{}".format(k), eval('args.MGCCZ4{}'.format(k)))
    for k, v in intparams.items():   solverconstants+= "static constexpr int {} = {};\n".format("MGCCZ4{}".format(k), eval('args.MGCCZ4{}'.format(k)))
    solverconstants+= "static constexpr int Scenario = {};\n".format(myscenario)

    my_solver.setSolverConstants(solverconstants)

    project.add_solver(my_solver)


    build_mode = modes[args.mode]
    
    dimensions = 3

    if args.periodic_bc:
      print( "Periodic BC set")
      periodic_boundary_conditions = [True,True,True]          # Periodic BC
    else:
      msg = "WARNING: Periodic BC deactivated"
      print(msg)
      periodic_boundary_conditions = [False,False,False]
      userwarnings.append((msg,None))

    project.set_global_simulation_parameters(
      dimensions,               # dimensions
      [-0.5, -0.5, -0.5],  [1.0, 1.0, 1.0],
      args.end_time,                 # end time
      0.0, args.plot_step_size,   # snapshots
      periodic_boundary_conditions
    )

    project.set_Peano4_installation("../../..", build_mode)

    #project.set_output_path( "/cosma6/data/dp004/dc-zhan3/tem" )
    #probe_point = [0,0,0]
    #project.add_plot_filter( probe_point,[0.0,0.0,0.0],1 )

    project.set_load_balancing("toolbox::loadbalancing::strategies::RecursiveSubdivision")

    peano4_project = project.generate_Peano4_project(verbose=True)

    peano4_project.output.makefile.add_CXX_flag( "-DMGCCZ4EINSTEIN" )
    peano4_project.output.makefile.add_cpp_file( "InitialValues.cpp" )
    peano4_project.output.makefile.add_cpp_file( "MGCCZ4Kernels.cpp" )

    # NOTE these lines are required to build with the fortran routines --- this will also require to uncomment some
    # includes etc
    #
    peano4_project.output.makefile.add_Fortran_flag( "-DMGCCZ4EINSTEIN -DDim3" )
    peano4_project.output.makefile.add_Fortran_flag( "-lstdc++ -fdefault-real-8 -fdefault-double-8 -cpp -std=legacy -ffree-line-length-512 -fPIC" )
    # peano4_project.output.makefile.add_CXX_flag( "-fPIE -DMGCCZ4EINSTEIN" )
    peano4_project.output.makefile.add_linker_flag( "-lstdc++ -fPIC -lgfortran" )
    # peano4_project.output.makefile.add_linker_flag( "-lstdc++ -fPIC -L/usr/lib/x86_64-linux-gnu -lgfortran" )

    # This might work for Intel (not tested)
    #peano4_project.output.makefile.add_Fortran_flag( "-r8 -cpp -auto -qopenmp-simd -O2" )
    #peano4_project.output.makefile.add_linker_flag( "-lstdc++ -fPIC" )
    # you might need -lifcore

    peano4_project.output.makefile.add_Fortran_module( "MainVariables.f90" )

    peano4_project.output.makefile.add_Fortran_files(
       ["PDE.f90 "]
     )


    # peano4_project.constants.export_string( "Scenario", "gaugewave-c++" )
    # peano4_project.constants.add_target_begin()
    # for k, v in floatparams.items(): peano4_project.constants.export_constexpr_with_type("MGCCZ4{}".format(k), eval('args.MGCCZ4{}'.format(k)), "double")
    # peano4_project.constants.add_target_end()
    # peano4_project.constants.add_target_begin()
    # for k, v in intparams.items():   peano4_project.constants.export_constexpr_with_type("MGCCZ4{}".format(k), eval('args.MGCCZ4{}'.format(k)), "int")
    # peano4_project.constants.add_target_end()
    # project.export_const_with_type("NumberOfUnknowns", 64, "int")
    #peano4_project.constants.export_string( "Scenario", "linearwave-c++" )

    peano4_project.generate( throw_away_data_after_generation=False )

    peano4_project.build( make_clean_first = True )

    # Remind the user of warnings!
    if len(userwarnings) >0:
        print("Please note that these warning occured before the build:")
        for msg, exception in userwarnings:
            if exception is None:
                print(msg)
            else: print(msg, "Exception: {}".format(exception))
