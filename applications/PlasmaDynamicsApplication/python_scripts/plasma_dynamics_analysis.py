from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#TODO: test DEM bounding box

import os
import sys
import math
import time as timer
import weakref

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidTransportApplication import *

from analysis_stage import AnalysisStage

""" import CFD_DEM_coupling
import swimming_DEM_procedures as SDP
import swimming_DEM_gid_output """
import CFD_DEM_for_plasma_dynamics_coupling
import plasma_dynamics_procedures
#import variables_management

def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()


# Import MPI modules if needed. This way to do this is only valid when using OpenMPI.
# For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    # import DEM_material_test_script_mpi as DEM_material_test_script
    Say('Running under MPI...........\n')
else:
    # DEM Application
    import DEM_procedures
    # import DEM_material_test_script
    Say('Running under OpenMP........\n')


class PlasmaDynamicsLogger(object):
    def __init__(self, do_print_file=False):
        self.terminal = sys.stdout
        self.console_output_file_name = 'console_output.txt'
        self.path_to_console_out_file = os.getcwd()
        self.path_to_console_out_file += '/' + self.console_output_file_name
        self.do_print_file = do_print_file
        if self.do_print_file:
            self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        if self.do_print_file:
            self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

    def getvalue(self):
        return self.terminal.getvalue()


class PlasmaDynamicsAnalysis(AnalysisStage):
    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __init__(self, model, parameters = Parameters("{}")):
        sys.stdout = PlasmaDynamicsLogger()
        self.StartTimer()
        self.model = model
        self.main_path = os.getcwd()
        self.project_parameters = parameters
        self.vars_man = variables_management.VariablesManager(self.project_parameters)

        # storing some frequently used variables
        self.time_step = self.project_parameters["MaxTimeStep"].GetDouble()
        self.end_time   = self.project_parameters["FinalTime"].GetDouble()
        self.do_print_results = self.project_parameters["do_print_results_option"].GetBool()

    # To-do: for the moment, provided for compatibility
    def _CreateSolver(self):
        import plasma_dynamics_solver
        return plasma_dynamics_solver.PlasmaDynamicsSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self.fluid_solution._GetSolver(),
                                                     self.disperse_phase_solution._GetSolver(),
                                                     self.vars_man)

