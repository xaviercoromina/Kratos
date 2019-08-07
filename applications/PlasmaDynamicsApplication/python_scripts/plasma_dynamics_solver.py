from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from python_solver import PythonSolver

# Import applications
import KratosMultiphysics.PlasmaDynamicsApplication as PlasmaDynamicsApplication
import plasma_dynamics_procedures 
import CFD_DEM_for_plasma_dynamics_coupling
import parameters_tools_for_plasma_dynamics as PT
import sys
lib_path = os.path.abspath(os.path.join(__file__, '..', 'SwimmingDEMApplication', 'python_scripts','derivative_recovery'))
sys.path.append(lib_path)

import derivative_recovery.derivative_recovery_strategy as derivative_recoverer
from mpmath import *

import math

def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()

class PlasmaDynamicsSolver(PythonSolver):
    def _ValidateSettings(self, project_parameters):

        default_processes_settings = Parameters("""{
                "python_module" : "calculate_nodal_area_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "CalculateNodalAreaProcess",
                "Parameters"    : {
                    "model_part_name" : "FluidModelPart",
                    "domain_size" : 3,
                    "fixed_mesh": false
                }
            }

        """)

        if not project_parameters["processes"].Has('non_optional_solver_processes'):
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

        else: # reconstruct non_optional_solver_processes list making sure calculate_nodal_area_process is not added twice
            non_optional_processes_list = list(project_parameters["processes"]["non_optional_solver_processes"])
            project_parameters["processes"].Remove("non_optional_solver_processes")
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

            for process in non_optional_processes_list:
                if process["python_module"].GetString() != 'calculate_nodal_area_process':
                    project_parameters["processes"]["non_optional_solver_processes"].Append(process)

        non_optional_solver_processes = project_parameters["processes"]["non_optional_solver_processes"]
        non_optional_solver_processes.Append(default_processes_settings)
        nodal_area_process_parameters = non_optional_solver_processes[non_optional_solver_processes.size() -1]["Parameters"]
        nodal_area_process_parameters["model_part_name"].SetString(self.fluid_solver.main_model_part.Name)
        nodal_area_process_parameters["domain_size"].SetInt(self.fluid_domain_dimension)

        if self.fluid_solver.settings.Has('move_mesh_flag'):
            the_mesh_moves = self.fluid_solver.settings["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)

        return project_parameters

    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        # Validate settings
        self.field_utility = field_utility
        self.vars_man = variables_manager
        self.fluid_domain_dimension = project_parameters["fluid_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.fluid_solver = fluid_solver
        self.dem_solver = dem_solver
        self.project_parameters = self._ValidateSettings(project_parameters)
        #self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        self.coupling_level_type = project_parameters["coupling"]["coupling_level_type"].GetInt()
        self.interaction_start_time = project_parameters["coupling"]["interaction_start_time"].GetDouble()
        self.integration_scheme = project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        self.fluid_dt = fluid_solver.settings["time_step"].GetDouble()
        self.do_solve_dem = project_parameters["custom_dem"]["do_solve_dem"].GetBool()
        self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool()

        self.step = 0
        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        self.first_DEM_iteration = True
        #self.ConstructStationarityTool()
        self.ConstructDerivativeRecoverer()
        self.ConstructHistoryForceUtility()
        # Call the base Python solver constructor
        super(PlasmaDynamicsSolver, self).__init__(model, project_parameters)


    def AddVariables(self):
        super(PlasmaDynamicsSolver,self).AddVariables()

        self.fluid_solver.main_model_part.ProcessInfo.SetValue(DENSITY, 0.0)  #should always be 0 when solving the Poisson equation
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(ABSORPTION_COEFFICIENT, 0.0) #should always be 0 when solving the Poisson equation
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(CONDUCTIVITY, 1.0) #should always be 1.0 when solving the Poisson equation
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(FLUID_FRACTION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_OLD)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_RATE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(ELECTRIC_POTENTIAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(ELECTRIC_FIELD)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(FLUID_ION_DENSITY)
        #self.dem_solver.spheres_model_part.AddNodalSolutionStepVariable(MACROPARTICLE_ION_DENSITY)
        #self.dem_solver.spheres_model_part.AddNodalSolutionStepVariable(PARTICLE_ION_VELOCITY)

        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(FLUID_FRACTION, 0, 1.0)
            node.SetSolutionStepValue(FLUID_FRACTION_OLD, 0, 1.0)


    def ConstructStationarityTool(self):
        """         self.stationarity = False
        self.stationarity_counter = self.GetStationarityCounter()
        self.stationarity_tool = plasma_dynamics_procedures.StationarityAssessmentTool(
            self.project_parameters["max_pressure_variation_rate_tol"].GetDouble(),
            plasma_dynamics_procedures.FunctionsCalculator()
            ) """
        pass

    def _ConstructProjectionModule(self):
        ###############
        #2D parameters:
        ###############
        self.h_min = 0.01 #TODO: 2D parameter only, this must be set from interface and the method must be checked for 2D
        n_balls = 1
        fluid_volume = 10
        # the variable n_particles_in_depth is only relevant in 2D problems
        self.project_parameters.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))

        #########################################################
        # creating a projection module for the fluid-DEM coupling
        #########################################################
        projection_module = CFD_DEM_for_plasma_dynamics_coupling.ProjectionModule(
        self.fluid_solver.main_model_part,
        self.dem_solver.spheres_model_part,
        self.dem_solver.all_model_parts.Get("RigidFacePart"),
        self.project_parameters,
        self.vars_man.coupling_dem_vars,
        self.vars_man.coupling_fluid_vars,
        self.vars_man.time_filtered_vars,
        flow_field=self.field_utility,
        domain_size=self.fluid_domain_dimension
        )

        projection_module.UpdateDatabase(self.h_min)

        return projection_module
        

    def ConstructDerivativeRecoverer(self):
        self.derivative_recovery_counter = self.GetRecoveryCounter()
        self.using_hinsberg_method = bool(self.project_parameters["basset_force_type"].GetInt() >= 3 or
                                          self.project_parameters["basset_force_type"].GetInt() == 1)
        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_solver.main_model_part,
            plasma_dynamics_procedures.FunctionsCalculator(self.fluid_domain_dimension)) 

    def ConstructHistoryForceUtility(self):
        """         self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        self.basset_force_tool = PlasmaDynamicsApplication.BassetForceTools() """
        pass

    def GetStationarityCounter(self):
        return plasma_dynamics_procedures.Counter(
            steps_in_cycle=self.project_parameters["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationary_problem_option"].GetBool())

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt())
        return plasma_dynamics_procedures.Counter(1, 1, there_is_something_to_recover)


    def GetHistoryForceQuadratureCounter(self):
        """         for prop in self.project_parameters["properties"].values():
            if prop["plasma_dynamics_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["plasma_dynamics_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()

                    return plasma_dynamics_procedures.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1) """

        return plasma_dynamics_procedures.Counter(is_dead=True)


    def AdvanceInTime(self, time):
        self.time = self.dem_solver.AdvanceInTime(time)
        # self.calculating_fluid_in_current_step = bool(time >= self.next_time_to_solve_fluid)
        # if self.calculating_fluid_in_current_step:
        #     self.next_time_to_solve_fluid = self.fluid_solver.AdvanceInTime(time)
        #     self.fluid_step += 1 
        self.step += 1 

        return self.time

    def UpdateALEMeshMovement(self, time): # TODO: move to derived solver
        pass

    def CalculateMinElementSize(self):
        return self.h_min

    def AssessStationarity(self):
        """         Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.Assess(self.fluid_solver.main_model_part)
        self.stationarity_counter.Deactivate(self.stationarity) """
        pass

    # Compute nodal quantities to be printed that are not generated as part of the
    # solution algorithm. For instance, the pressure gradient, which is not used for
    # the coupling but can be of interest.
    def ComputePostProcessResults(self):
        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            self._GetProjectionModule().ComputePostProcessResults(self.dem_solver.spheres_model_part.ProcessInfo)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step


    def Predict(self):
        # if self.CannotIgnoreFluidNow():
        #     self.fluid_solver.Predict() 
        pass


    def ApplyForwardCoupling(self, alpha='None'):
        self._GetProjectionModule().ApplyForwardCoupling(alpha)


    def _GetProjectionModule(self):
        if not hasattr(self, 'projection_module'):
            self.projection_module = self._ConstructProjectionModule()
        return self.projection_module


    def InitializeSolutionStep(self):
        pass

    def ProjectFromParticles(self):
        self._GetProjectionModule().ProjectFromParticles()

    def SolveSolutionStep(self):
        # update possible movements of the fluid mesh
        self.UpdateALEMeshMovement(self.time)

        # Solving the fluid part
        #Say('Solving Fluid... (', self.fluid_solver.main_model_part.NumberOfElements(0), 'elements )\n')
        #self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool() 

        # if self.CannotIgnoreFluidNow():
        #     self.SolveFluidSolutionStep()
        # else:
        #     Say("Skipping solving system for the fluid phase...\n")

        # Check for stationarity: this is useful for steady-state problems, so that
        # the calculation stops after reaching the solution.
        """         if self.stationarity_counter.Tick():
            self.AssessStationarity() """

        self.derivative_recovery_counter.Activate(self.time > self.interaction_start_time and self.calculating_fluid_in_current_step)

        if self.derivative_recovery_counter.Tick():
            self.recovery.Recover()

        # Solving the disperse-phase component
        Say('Solving DEM... (', self.dem_solver.spheres_model_part.NumberOfElements(0), 'elements )')
        self.SolveDEM()

    def SolveFluidSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        

    def SolveDEMSolutionStep(self):
        self.dem_solver.SolveSolutionStep()

    def SolveDEM(self):
        #self.PerformEmbeddedOperations() TO-DO: it's crashing

        #Forward Coupling
        # it_is_time_to_forward_couple = (self.time >= self.interaction_start_time
        #                                 and self.coupling_level_type)

        # alpha = 1.0 - (self.next_time_to_solve_fluid - self.time) / self.fluid_dt
        

        # if it_is_time_to_forward_couple or self.first_DEM_iteration:
        #     self.ApplyForwardCoupling(alpha)

        """         if self.quadrature_counter.Tick():
            self.AppendValuesForTheHistoryForce() """

        # Performing the time integration of the DEM part

        if self.do_solve_dem:
            self.SolveDEMSolutionStep()

        self.first_DEM_iteration = False

    def FinalizeSolutionStep(self):
        # mp.dps=30 #digit precision for calculation of exp

        # #Backward Coupling
        # for node in self.fluid_solver.main_model_part.Nodes:

        #     #Getting the electric potential Phi using the variable TEMPERATURE in the FluidTransportApplication
        #     electric_potential = node.GetSolutionStepValue(TEMPERATURE)  #TODO: inverse the way
        #     node.SetSolutionStepValue(ELECTRIC_POTENTIAL, electric_potential)

        #     #Calculate electron density on each fluid model node TODO: put it in C++
        #     n_0 = 2.0*1e11 # electron density when Phi = 0, n_electron = n_0 * 10^(n_1) TODO: put it in the json
        #     #n_1 = 0.0 
        #     # e = 1.60*10**(-19) C  Coulomb charge
        #     # k = 1.38*10**(-23) J/K  Boltzmann constant
        #     # T_e = 58025 K  (= 5eV) Electron temperature TODO: put it in the json
        #     # n_e = n_0 * exp(e * Phi / (k * T_e))
        #     electric_constant = 0.2  # e / (k * T_e)
        #     fluid_electron_density = n_0 * exp(electric_constant * electric_potential)

        #     node.SetSolutionStepValue(FLUID_ELECTRON_DENSITY, fluid_electron_density)


        #     fluid_ion_density = node.GetSolutionStepValue(FLUID_ION_DENSITY)  
        #     #RHS = 0.0
        #     RHS = 1.81*(10**(-8))*(fluid_electron_density-fluid_ion_density)
        #     print("RHS is equal to:")
        #     print(RHS)
        #     node.SetSolutionStepValue(HEAT_FLUX, RHS)

        for node in self.dem_solver.spheres_model_part.Nodes:
            particle_ion_velocity = node.GetSolutionStepValue(VELOCITY)
            node.SetSolutionStepValue(PARTICLE_ION_VELOCITY, particle_ion_velocity)


    def AppendValuesForTheHistoryForce(self):
        if self.using_hinsberg_method:
            self.basset_force_tool.AppendIntegrandsWindow(self.dem_solver.spheres_model_part)
        elif self.project_parameters["basset_force_type"].GetInt() == 2:
            self.basset_force_tool.AppendIntegrands(self.dem_solver.spheres_model_part)

    def ImportModelPart(self): # TODO: implement this
        pass

    def GetComputingModelPart(self):
        return self.dem_solver.spheres_model_part