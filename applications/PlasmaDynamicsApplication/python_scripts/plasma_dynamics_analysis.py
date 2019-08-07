from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#TODO: test DEM bounding box

import os
import sys
import math
import time as timer
import weakref
from mpmath import *

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *

from analysis_stage import AnalysisStage

import CFD_DEM_for_plasma_dynamics_coupling
import plasma_dynamics_procedures
import variables_management_for_plasma_dynamics
import plasma_dynamics_gid_output

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

        self.SetProjectParameters(parameters)

        self.vars_man = variables_management_for_plasma_dynamics.VariablesManager(self.project_parameters)

        self._GetDEMAnalysis().coupling_analysis = weakref.proxy(self)

        self._GetFluidAnalysis().coupling_analysis = weakref.proxy(self)

        self.procedures = weakref.proxy(self._GetDEMAnalysis().procedures)
        self.report = DEM_procedures.Report()

        self._GetDEMAnalysis().SetAnalyticFaceWatcher()

        # defining member variables for the model_parts (for convenience)
        self.fluid_model_part = self._GetFluidAnalysis().fluid_model_part
        self.spheres_model_part = self._GetDEMAnalysis().spheres_model_part
        self.cluster_model_part = self._GetDEMAnalysis().cluster_model_part
        self.rigid_face_model_part = self._GetDEMAnalysis().rigid_face_model_part
        self.dem_inlet_model_part = self._GetDEMAnalysis().dem_inlet_model_part
        self.vars_man.ConstructListsOfVariables(self.project_parameters)
        super(PlasmaDynamicsAnalysis, self).__init__(model, self.project_parameters) # TODO: The DEM jason is now interpreted as the coupling json. This must be changed

    def StartTimer(self):
        self.timer = timer
        self.simulation_start_time = timer.time()


    def SetProjectParameters(self, parameters):
        self.project_parameters = parameters
        self.time_step = self.project_parameters["time_stepping"]["time_step"].GetDouble()
        self.end_time   = self.project_parameters["problem_data"]["end_time"].GetDouble()
        self.do_print_results = self.project_parameters["do_print_results_option"].GetBool()
        self.fluid_parameters = self.project_parameters['fluid_parameters']
        self.time_to_calculate_poisson_step = self.project_parameters["time_stepping"]["time_to_calculate_poisson_step"].GetDouble()
        self.full_PIC_option = self.project_parameters["full_PIC_option"].GetBool()

        self.is_the_first_time = True

        if self.full_PIC_option:
            self.electron_time_step = self.project_parameters["time_stepping"]["electron_time_step"].GetDouble()

        # First, read the parameters generated from the interface
        import plasma_dynamics_default_input_parameters as only_plasma_dynamics_defaults
        self.project_parameters.ValidateAndAssignDefaults(only_plasma_dynamics_defaults.GetDefaultInputParameters())

        # Second, set the default 'beta' parameters (candidates to be moved to the interface)
        self.SetBetaParameters()

        # Third, make sure the parameters passed to the different (sub-)analyses is coherent
        # with the general parameters
        self.ModifyInputParametersForCoherence()



    # Set input parameters that have not yet been transferred to the interface
    # import the configuration data as read from the GiD
    def SetBetaParameters(self):
        Add = self.project_parameters.AddEmptyValue
        if self.project_parameters["custom_dem"]["type_of_dem_inlet"].GetString() == 'ForceImposed':
            Add("inlet_force_vector").SetVector(Vector([0., 0., 1.])) # TODO: generalize

        # Setting body_force_per_unit_mass_variable_name
        Add("body_force_per_unit_mass_variable_name").SetString('BODY_FORCE')


    # This step is added to allow modifications to the possibly incompatibilities
    # between the individual parameters coming from each sub-application
    # (i.e., fluid and dem apps)
    def ModifyInputParametersForCoherence(self):
        # Making all time steps exactly commensurable
        output_time = self.project_parameters["output_interval"].GetDouble()
        self.output_time = int(output_time / self.time_step) * self.time_step
        self.project_parameters["output_interval"].SetDouble(self.output_time)

        if self.full_PIC_option:
            if self.time_step != 0.0 and self.electron_time_step != 0.0:
                self.electron_time_step = self.time_step / (int(self.time_step / self.electron_time_step))

        self.fluid_time_step = self.fluid_parameters["solver_settings"]["time_step"].GetDouble()
        self.fluid_end_time = self.fluid_parameters["problem_data"]["end_time"].GetDouble()

        # if self.fluid_time_step < self.time_step:
        #     error_message = ('The fluid time step (' + str(self.fluid_time_step)
        #                      + ') must be larger or equal than the overall time step (' + str(self.time_step)
        #                      + ')!')
        #     raise Exception(error_message)


        # self.fluid_time_step = int(self.fluid_time_step / self.time_step) * self.time_step
        # self.fluid_parameters["solver_settings"]["time_step"].SetDouble(self.fluid_time_step)
        self.project_parameters["dem_parameters"]["MaxTimeStep"].SetDouble(self.time_step)
        translational_scheme_name = self.project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        self.project_parameters["dem_parameters"]["TranslationalIntegrationScheme"].SetString(translational_scheme_name)

        time_steps_per_stationarity_step = self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt()
        self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].SetInt(max(1, int(time_steps_per_stationarity_step)))

        if self.project_parameters["coupling"]["coupling_level_type"].GetInt() > 1:
            self.project_parameters["stationarity"]["stationary_problem_option"].SetBool(False)

        self.SetDoSolveDEMVariable()


    def SetDoSolveDEMVariable(self):
        self.do_solve_dem = self.project_parameters["custom_dem"]["do_solve_dem"].GetBool()


    def ReadDispersePhaseModelParts(self,
                                    starting_node_Id=0,
                                    starting_elem_Id=0,
                                    starting_cond_Id=0):
        creator_destructor = self._GetDEMAnalysis().creator_destructor
        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(self.fluid_model_part)  #self.fluid_model_part
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(self.fluid_model_part)
        max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(self.fluid_model_part)
        self._GetDEMAnalysis().BaseReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def Run(self):
        super(PlasmaDynamicsAnalysis, self).Run()

        return self.GetReturnValue()

    def GetReturnValue(self):
        return 0.0

    def Initialize(self):
        Say('Initializing simulation...\n')
        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        if self.do_print_results:
            [self.post_path, data_and_results, self.graphs_path, MPI_results] = \
            self.procedures.CreateDirectories(str(self.main_path),
                                            str(self.project_parameters["problem_data"]["problem_name"].GetString()),
                                            self.run_code)
            plasma_dynamics_procedures.CopyInputFilesIntoFolder(self.main_path, self.post_path)
            self.MPI_results = MPI_results

        self.FluidInitialize()

        self.DispersePhaseInitialize()

        if self.full_PIC_option:
            self.FillParticlesSubModelParts()

        self.SetAllModelParts()

        if self.project_parameters.Has('plasma_dynamics_output_processes') and self.do_print_results:
            gid_output_options = self.project_parameters["plasma_dynamics_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            write_conditions_option = result_file_configuration["gidpost_flags"]["WriteConditionsFlag"].GetString() == "WriteConditions"
            deformed_mesh_option = result_file_configuration["gidpost_flags"]["WriteDeformedMeshFlag"].GetString() == "WriteDeformed"
            old_gid_output_post_options_dict = {'GiD_PostAscii':'Ascii','GiD_PostBinary':'Binary','GiD_PostAsciiZipped':'AsciiZipped'}
            old_gid_output_multiple_file_option_dict = {'SingleFile':'Single','MultipleFiles':'Multiples'}
            post_mode_key = result_file_configuration["gidpost_flags"]["GiDPostMode"].GetString()
            multiple_files_option_key = result_file_configuration["gidpost_flags"]["MultiFileFlag"].GetString()

            self.plasma_dynamics_gid_io = \
            plasma_dynamics_gid_output.PlasmaDynamicsGiDOutput(
                file_name = self.project_parameters["problem_data"]["problem_name"].GetString(),
                vol_output = result_file_configuration["body_output"].GetBool(),
                post_mode = old_gid_output_post_options_dict[post_mode_key],
                multifile = old_gid_output_multiple_file_option_dict[multiple_files_option_key],
                deformed_mesh = deformed_mesh_option,
                write_conditions = write_conditions_option
                )

            self.plasma_dynamics_gid_io.initialize_plasma_dynamics_results(
                self.spheres_model_part,
                self.cluster_model_part,
                self.rigid_face_model_part,
                self.mixed_model_part)


        super(PlasmaDynamicsAnalysis, self).Initialize()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        """         plasma_dynamics_procedures.ApplySimilarityTransformations(
            self.fluid_model_part,
            self.project_parameters["similarity_transformation_type"].GetInt(),
            self.project_parameters["model_over_real_diameter_factor"].GetDouble()
            ) """

        if self.do_print_results:
            self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = plasma_dynamics_procedures.IOTools(self.project_parameters)

        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(
            self.spheres_model_part)

        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            default_meso_scale_length_needed = (
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].GetDouble() <= 0.0 and
                self.spheres_model_part.NumberOfElements(0) > 0)

            if default_meso_scale_length_needed:
                biggest_size = (2 * dem_physics_calculator.CalculateMaxNodalVariable(self.spheres_model_part, Kratos.RADIUS))
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].SetDouble(20 * biggest_size)

            elif self.spheres_model_part.NumberOfElements(0) == 0:
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].SetDouble(1.0)

        # creating a custom functions calculator for the implementation of
        # additional custom functions
        fluid_domain_dimension = self.project_parameters["fluid_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.custom_functions_tool = plasma_dynamics_procedures.FunctionsCalculator(fluid_domain_dimension)

        # creating a stationarity assessment tool
        """         self.stationarity_tool = plasma_dynamics_procedures.StationarityAssessmentTool(
            self.project_parameters["max_pressure_variation_rate_tol"].GetDouble(),
            self.custom_functions_tool
            ) """

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()


        Say('Initialization Complete\n')

        ##################################################

        #    I N I T I A L I Z I N G    T I M E    L O O P

        ##################################################
        self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
        self.interaction_start_time = self.project_parameters["coupling"]["interaction_start_time"].GetDouble()

        #fluid_time is independent from the ion/electron time line, it is the own time parameter of the Poisson equation solver
        #time_to_calculate_poisson is the time parameter belonging to the ion/electron time line which says when to solve the space charge equation (Poisson)
        self.fluid_time = self.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self.time_to_calculate_poisson = self.project_parameters["time_stepping"]["start_time_to_calculate_poisson"].GetDouble()

        self.time_to_calculate_Boltzmann_electron_density = self.project_parameters["time_stepping"]["start_time_Boltzmann_electron_density"].GetDouble()

        if self.full_PIC_option:
            self.electron_time = self.project_parameters["time_stepping"]["electron_start_time"].GetDouble()

        self.fluid_time_step = self._GetFluidAnalysis()._GetSolver().ComputeDeltaTime()
        self.time_step = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.time_step
        self.cluster_model_part.ProcessInfo[DELTA_TIME] = self.time_step

        if self.full_PIC_option:
            self.electron_model_part.ProcessInfo[DELTA_TIME] = self.time_step
            self.ion_model_part.ProcessInfo[DELTA_TIME] = self.time_step
        
        
        self.stationarity = False
    
        self.upstream_potential = self.project_parameters["upstream_potential"].GetDouble()
        self.downstream_potential = self.project_parameters["downstream_potential"].GetDouble()
        self.upstream_electron_density = self.project_parameters["upstream_electron_density"].GetDouble()
        self.downstream_electron_density = self.project_parameters["downstream_electron_density"].GetDouble()    

        self.coupling_level_type = self.project_parameters["coupling"]["coupling_level_type"].GetInt()

        # setting up loop counters:
        self.DEM_to_fluid_counter = self.GetBackwardCouplingCounter()
        #self.stationarity_counter = self.GetStationarityCounter()
        self.print_counter = self.GetPrintCounter()
        self.debug_info_counter = self.GetDebugInfo()
        self.particles_results_counter = self.GetParticlesResultsCounter()
        self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        # Phantom
        self._GetDEMAnalysis().analytic_data_counter = self.ProcessAnalyticDataCounter()
        self.mat_deriv_averager           = plasma_dynamics_procedures.Averager(1, 3)
        self.laplacian_averager           = plasma_dynamics_procedures.Averager(1, 3)

        self.report.total_steps_expected = int(self.end_time / self.time_step)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DEM_procedures.PostUtils(self.project_parameters['dem_parameters'], self.spheres_model_part)


        plasma_dynamics_procedures.InitializeVariablesWithNonZeroValues(self.project_parameters, 
                                                                        self.fluid_model_part, 
                                                                        self.spheres_model_part)

        # ANALYTICS BEGIN
        self.project_parameters.AddEmptyValue("perform_analytics_option").SetBool(False)

        """         if self.project_parameters["perform_analytics_option"].GetBool():
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(
                self.fluid_model_part,
                self.fluid_time_step,
                self.end_time,
                variables_to_measure,
                steps_between_measurements
                )
            point_coors = [0.0, 0.0, 0.01]
            target_node = plasma_dynamics_procedures.FindClosestNode(self.fluid_model_part, point_coors)
            target_id = target_node.Id
            Say(target_node.X, target_node.Y, target_node.Z)
            Say(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            Say(gauge.variables) """
        # ANALYTICS END

        """         import derivative_recovery.derivative_recovery_strategy as derivative_recoverer

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_model_part,
            self.custom_functions_tool) """

        self.FillHistoryForcePrecalculatedVectors()


        if self.do_print_results:
            self._Print()


    def GetRunCode(self):
        return ""


    def FluidInitialize(self):
        self.fluid_model_part = self._GetFluidAnalysis().fluid_model_part
        self._GetFluidAnalysis().vars_man = self.vars_man
        self._GetFluidAnalysis().Initialize()

        self.AddExtraProcessInfoVariablesToFluid()

        plasma_dynamics_procedures.AddExtraDofs(self.fluid_model_part,
                         self.spheres_model_part,
                         self.cluster_model_part,
                         self.dem_inlet_model_part,
                         self.vars_man)

    def AddExtraProcessInfoVariablesToFluid(self):
        self.vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.project_parameters, self.fluid_model_part)

    def DispersePhaseInitialize(self):
        self.vars_man.__class__.AddNodalVariables(self.spheres_model_part, self.vars_man.dem_vars)
        self.vars_man.__class__.AddNodalVariables(self.rigid_face_model_part, self.vars_man.rigid_faces_vars)
        self.vars_man.__class__.AddNodalVariables(self.dem_inlet_model_part, self.vars_man.inlet_vars)
        self.vars_man.AddExtraProcessInfoVariablesToDispersePhaseModelPart(self.project_parameters,
                                                                           self._GetDEMAnalysis().spheres_model_part)

        self._GetDEMAnalysis().Initialize()

    def FillParticlesSubModelParts(self):
        info_particle = 0.0
        if not self.spheres_model_part.HasSubModelPart("ElectronParticlePart"):
            self.spheres_model_part.CreateSubModelPart('ElectronParticlePart')
        self.electron_model_part = self.spheres_model_part.GetSubModelPart('ElectronParticlePart')

        electron_particle_ids = [elem.Id for elem in self.spheres_model_part.Elements if (elem.Calculate(ELECTRIC_FIELD_PROJECTED_TO_PARTICLE, info_particle, self.spheres_model_part.ProcessInfo)==2.0)]
        self.electron_model_part.AddElements(electron_particle_ids)

        if not self.spheres_model_part.HasSubModelPart("IonParticlePart"):
            self.spheres_model_part.CreateSubModelPart('IonParticlePart')
        self.ion_model_part = self.spheres_model_part.GetSubModelPart('IonParticlePart')
        ion_particle_ids = [elem.Id for elem in self.spheres_model_part.Elements if (elem.Calculate(ELECTRIC_FIELD_PROJECTED_TO_PARTICLE, info_particle, self.spheres_model_part.ProcessInfo)==1.0)]
        self.ion_model_part.AddElements(ion_particle_ids)

    def SetAllModelParts(self):
        self.all_model_parts = weakref.proxy(self._GetDEMAnalysis().all_model_parts)

        # defining a fluid model
        self.all_model_parts.Add(self.fluid_model_part) 

        # defining a model part for the mixed part
        self.all_model_parts.Add(self.model.CreateModelPart("MixedPart"))

        self.mixed_model_part = self.all_model_parts.Get('MixedPart')


    def SetPostUtils(self):
          # creating a Post Utils object that executes several post-related tasks
        self.post_utils = plasma_dynamics_procedures.PostUtils(self.plasma_dynamics_gid_io,
                                        self.project_parameters,
                                        self.vars_man,
                                        self.fluid_model_part,
                                        self.spheres_model_part,
                                        self.cluster_model_part,
                                        self.rigid_face_model_part,
                                        self.mixed_model_part)



    def GetBackwardCouplingCounter(self):
        return plasma_dynamics_procedures.Counter(1, 1, self.project_parameters["coupling"]["coupling_level_type"].GetInt() > 1)

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt()) 
        return plasma_dynamics_procedures.Counter(1, 1, there_is_something_to_recover)

    def GetStationarityCounter(self):
        return plasma_dynamics_procedures.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationarity"]["stationary_problem_option"].GetBool())

    def GetPrintCounter(self):
        counter = plasma_dynamics_procedures.Counter(steps_in_cycle=int(self.output_time / self.time_step + 0.5),
                              beginning_step=int(self.output_time / self.time_step),
                              is_dead = not self.do_print_results)
        return counter

    def GetDebugInfo(self):
        return plasma_dynamics_procedures.Counter(
            self.project_parameters["debug_tool_cycle"].GetInt(),
            1,
            self.project_parameters["print_debug_info_option"].GetBool())

    def GetParticlesResultsCounter(self):
        return plasma_dynamics_procedures.Counter(
            self.project_parameters["print_particles_results_cycle"].GetInt(),
            1,
            self.project_parameters["print_particles_results_option"].GetBool())

    def GetHistoryForceQuadratureCounter(self):
        """         for prop in self.project_parameters["properties"].values():
            if prop["plasma_dynamics_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["plasma_dynamics_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()

                    return plasma_dynamics_procedures.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1) """

        return plasma_dynamics_procedures.Counter(is_dead=True)

    def ProcessAnalyticDataCounter(self):
        return plasma_dynamics_procedures.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_analytic_processing_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["do_process_analytic_data"].GetBool())

    def GetVolumeDebugTool(self):
        return plasma_dynamics_procedures.ProjectionDebugUtils(
            self.project_parameters["fluid_domain_volume"].GetDouble(),
            self.fluid_model_part,
            self.spheres_model_part,
            self.custom_functions_tool) 



    def FillHistoryForcePrecalculatedVectors(self): # TODO: more robust implementation
        # Warning: this estimation is based on a constant time step for DEM.
        # This is usually the case, but could not be so.
        """         for prop in self.project_parameters["properties"].values():
            if prop["plasma_dynamics_law_parameters"].Has("history_force_parameters"):

                if prop["plasma_dynamics_law_parameters"]["history_force_parameters"]["name"].GetString() != 'default':
                    total_number_of_steps = int(self.end_time / self.project_parameters["MaxTimeStep"].GetDouble()) + 20
                    history_force_parameters = prop["plasma_dynamics_law_parameters"]["history_force_parameters"]
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()
                    self._GetSolver().basset_force_tool.FillDaitcheVectors(
                        total_number_of_steps,
                        history_force_parameters["quadrature_order"].GetInt(),
                        time_steps_per_quadrature_step)

                    if history_force_parameters.Has("mae_parameters"):
                        mae_parameters = history_force_parameters["mae_parameters"]
                        time_window = mae_parameters["window_time_interval"].GetDouble()
                        quadrature_dt = time_steps_per_quadrature_step * self.time_step
                        number_of_quadrature_steps_in_window = int(time_window / quadrature_dt)
                        if mae_parameters["do_use_mae"].GetBool():
                            self._GetSolver().basset_force_tool.FillHinsbergVectors(
                            self.spheres_model_part,
                            mae_parameters["m"].GetInt(),
                            number_of_quadrature_steps_in_window)
                            break """
        pass


    def _Print(self):
        os.chdir(self.post_path)

        if self.particles_results_counter.Tick():
            self.io_tools.PrintParticlesResults(
                self.vars_man.variables_to_print_in_file,
                self.time,
                self.spheres_model_part)

        self.post_utils.Writeresults(self.time)
        os.chdir(self.main_path)


    def InitializeSolutionStep(self):
        self.TellTime()
        if self.time >= self.time_to_calculate_poisson:
            if self.project_parameters["dem_parameters"]["dem_inlet_option"].GetBool(): 
                plasma_dynamics_procedures.InitializeVariablesWithNonZeroValues(self.project_parameters, 
                                                                        self.fluid_model_part, 
                                                                        self.spheres_model_part)
            
            # applying DEM-to-fluid coupling (backward coupling) for ion density

            if self.is_the_first_time:
                self.is_the_first_time = False
            else:
                if self.DEM_to_fluid_counter.Tick() and self.time >= self.interaction_start_time:
                    #self._GetSolver()._GetProjectionModule().ProjectFromParticles()
                    self._GetSolver().ProjectFromParticles()  
                if self.time >= self.time_to_calculate_Boltzmann_electron_density: 
                    self.DoBackwardCoupling()
            self.SolveFluid()
            self.SetElectricPotential()
            self.time_to_calculate_poisson += self.time_to_calculate_poisson_step
        
        self._GetDEMAnalysis().InitializeSolutionStep()
        # if self._GetSolver().CannotIgnoreFluidNow():
        #     self._GetFluidAnalysis().InitializeSolutionStep()
        super(PlasmaDynamicsAnalysis, self).InitializeSolutionStep()
        
        #Forward Coupling (before DEM solving)
        it_is_time_to_forward_couple = (self.time >= self.interaction_start_time
                                        and self.coupling_level_type)

        alpha = 1.0 - (self.time_to_calculate_poisson - self.time) / self.time_to_calculate_poisson_step
        #alpha = 1.0
        

        if it_is_time_to_forward_couple:
            self.ApplyForwardCoupling(alpha)


    def SetElectricPotential(self):
        for node in self.fluid_model_part.Nodes:

            #Getting the electric potential Phi using the variable TEMPERATURE in the FluidTransportApplication
            electric_potential = node.GetSolutionStepValue(TEMPERATURE)  #TODO: inverse the way
            node.SetSolutionStepValue(ELECTRIC_POTENTIAL, electric_potential)

    def DoBackwardCoupling(self):
        mp.dps=30 #digit precision for calculation of exp, otherwise gives an overflow error

        #Backward Coupling
        for node in self.fluid_model_part.Nodes:

            #Calculate electron density on each fluid model node TODO: put it in C++
            #n_0 = 1.0*1e4 # electron density when Phi = 0, n_electron = n_0 * 10^(n_1) TODO: put it in the json
            #n_1 = 0.0 
            # e = 1.60*10**(-19) C  Coulomb charge
            # k = 1.38*10**(-23) J/K  Boltzmann constant
            # T_e = 58025 K  (= 5eV) Electron temperature TODO: put it in the json
            # n_e = n_0 * exp(e * Phi / (k * T_e))
            electric_constant = 0.2  # e / (k * T_e)

            fluid_ion_density = node.GetSolutionStepValue(FLUID_ION_DENSITY)

            electric_potential = node.GetSolutionStepValue(ELECTRIC_POTENTIAL)
            
            #Downstream Region
            if node.Z > 0.01495:
                if self.downstream_potential >= electric_potential:
                    fluid_electron_density = self.downstream_electron_density * exp(electric_constant * (electric_potential-self.downstream_potential))
                else:
                    fluid_electron_density = self.downstream_electron_density * (1 + electric_constant * (electric_potential-self.downstream_potential))
                RHS = -1.81*(10**(-8))*(fluid_electron_density-fluid_ion_density)
                node.SetSolutionStepValue(HEAT_FLUX, RHS)
            #Grid Region
            elif node.Z <= 0.01495 and node.Z >= 0.0072:
                fluid_electron_density = 0
                RHS = 0.0
                node.SetSolutionStepValue(HEAT_FLUX, RHS)
            #Upstream Region
            else:
                if self.upstream_potential >= electric_potential:
                    fluid_electron_density = self.upstream_electron_density * exp(electric_constant * (electric_potential-self.upstream_potential))
                else:
                    fluid_electron_density = self.upstream_electron_density * (1 + electric_constant * (electric_potential-self.upstream_potential))
                RHS = -1.81*(10**(-8))*(fluid_electron_density-fluid_ion_density)
                node.SetSolutionStepValue(HEAT_FLUX, RHS)
            #fluid_electron_density = n_0 * exp(electric_constant * electric_potential)

            node.SetSolutionStepValue(FLUID_ELECTRON_DENSITY, fluid_electron_density)

    def SolveFluid(self):
        self.SolveFluidCounter = 0
        Say('Solving Fluid... (', self._GetFluidAnalysis()._GetSolver().main_model_part.NumberOfElements(0), 'elements )\n')
        while self.KeepAdvancingSolutionLoopForFluid() and self.SolveFluidCounter<110:
            self.fluid_time = self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo.GetValue(TIME)
            self.fluid_time = self._GetFluidAnalysis()._GetSolver().AdvanceInTime(self.fluid_time)
            self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo.SetValue(TIME,self.fluid_time)
            fluid_dt = self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo[DELTA_TIME]
            self._GetSolver().fluid_step += 1 
            self.TellFluidTime()
            self._GetFluidAnalysis().InitializeSolutionStep()
            self._GetFluidAnalysis()._GetSolver().Predict()
            is_converged = self._GetFluidAnalysis()._GetSolver().SolveSolutionStep()
            self._GetFluidAnalysis().CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self._GetFluidAnalysis().FinalizeSolutionStep()
            self.SolveFluidCounter += 1
        #self._GetFluidAnalysis().OutputSolutionStep()
        #self._GetSolver().fluid_step = 0 
        self.fluid_time = self.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo.SetValue(TIME,self.fluid_time)
        self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo.SetValue(DELTA_TIME,
                                                  self.fluid_parameters["solver_settings"]["time_step"].GetDouble())
        self._GetFluidAnalysis()._GetSolver().main_model_part.ProcessInfo.SetValue(STEP, 0)

    def KeepAdvancingSolutionLoopForFluid(self):
        return self.fluid_time < self.fluid_end_time

    def TellTime(self):
        Say('DEM time: ', str(self.time) + ', step: ', self._GetSolver().step)
        Say('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time, '\n')

    def TellFluidTime(self):
        Say('fluid time: ', str(self.fluid_time) + ', step: ', self._GetSolver().fluid_step)
        Say('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time, '\n')


    def ApplyForwardCoupling(self, alpha='None'):
        #self._GetSolver()._GetProjectionModule().ApplyForwardCoupling(alpha)
        self._GetSolver().ApplyForwardCoupling(alpha)


    def FinalizeSolutionStep(self):
        # printing if required
        # if self._GetSolver().CannotIgnoreFluidNow():
        #     self._GetFluidAnalysis().FinalizeSolutionStep() 

        self._GetDEMAnalysis().FinalizeSolutionStep()

        # applying DEM-to-fluid coupling (backward coupling)

        # if self.DEM_to_fluid_counter.Tick() and self.time >= self.project_parameters["coupling"]["interaction_start_time"].GetDouble():
        #     self._GetSolver().projection_module.ProjectFromParticles() 

        # coupling checks (debugging)
        if self.debug_info_counter.Tick():
            self.dem_volume_tool.UpdateDataAndPrint(
                self.project_parameters["fluid_domain_volume"].GetDouble()) 

        super(PlasmaDynamicsAnalysis, self).FinalizeSolutionStep()

    def OutputSolutionStep(self):
        # printing if required

        if self.print_counter.Tick():
            self.ComputePostProcessResults()
            self._Print()

        super(PlasmaDynamicsAnalysis, self).OutputSolutionStep()

    def ComputePostProcessResults(self):
        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            #self._GetSolver()._GetProjectionModule().ComputePostProcessResults(self.spheres_model_part.ProcessInfo)
            self._GetSolver().ComputePostProcessResults()
    def SetInletWatcher(self):
        self.watcher_analyser.SetInlet(self.DEM_inlet)

    def SetInlet(self):
        if self.project_parameters["dem_inlet_option"].GetBool():
            # Constructing the inlet and initializing it
            # (must be done AFTER the self.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible.
            # This should be generalized.
            if self.project_parameters["type_of_dem_inlet"].GetString() == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.dem_inlet_model_part)
            elif self.project_parameters["type_of_dem_inlet"].GetString() == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.dem_inlet_model_part, self.project_parameters["inlet_force_vector"].GetVector())

            self._GetDEMAnalysis().DEM_inlet = self.DEM_inlet
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self._GetDEMAnalysis().creator_destructor)


    def SetAnalyticParticleWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(
            analytic_particle_watcher=self.particle_watcher,
            path=self.main_path)


    def Finalize(self):
        Say('Finalizing simulation...\n')
        if self.do_print_results:
            self.plasma_dynamics_gid_io.finalize_results()

        self.PerformFinalOperations(self.time)

        self._GetFluidAnalysis().Finalize()

        self.TellFinalSummary(self.time, self._GetSolver().step, self._GetSolver().fluid_step)



    def PerformFinalOperations(self, time=None):
        os.chdir(self.main_path)

        if self.do_print_results:
            del self.post_utils
            


    def TellFinalSummary(self, time, dem_step, fluid_step, message_n_char_width=60):
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time

        if simulation_elapsed_time and dem_step and fluid_step: 
            elapsed_time_per_unit_dem_step = simulation_elapsed_time / dem_step
            elapsed_time_per_unit_fluid_step = simulation_elapsed_time / fluid_step

        else:
            elapsed_time_per_unit_dem_step = 0.0
            elapsed_time_per_unit_fluid_step = 0.0

        final_message = ('\n\n'
                         + '*' * message_n_char_width + '\n'
                         + 'CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.' + '\n'
                         + 'Total number of DEM steps run: ' + str(dem_step) + '\n'
                         + 'Total number of fluid steps run: ' + str(fluid_step) + '\n'
                         + 'Elapsed time: ' + '%.5f'%(simulation_elapsed_time) + ' s ' + '\n'
                         + ',, per fluid time step: ' + '%.5f'%(elapsed_time_per_unit_fluid_step) + ' s ' + '\n'
                         + ',, per DEM time step: ' + '%.5f'%(elapsed_time_per_unit_dem_step) + ' s ' + '\n'
                         + '*' * message_n_char_width + '\n')

        Say(final_message)


    # To-do: for the moment, provided for compatibility
    def _CreateSolver(self):
        import plasma_dynamics_solver
        return plasma_dynamics_solver.PlasmaDynamicsSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self._GetFluidAnalysis()._GetSolver(),
                                                     self._GetDEMAnalysis()._GetSolver(),
                                                     self.vars_man)

    def _GetDEMAnalysis(self):
        if not hasattr(self, '_disperse_phase_analysis'):
            import fluid_coupled_DEM_plasma_dynamics_analysis as DEM_analysis
            self._disperse_phase_analysis = DEM_analysis.FluidCoupledDEMPDAnalysisStage(self.model, self.project_parameters)

        return self._disperse_phase_analysis

    def _GetFluidAnalysis(self):
        if not hasattr(self, '_fluid_phase_analysis'):
            import DEM_coupled_fluid_plasma_dynamics_analysis as fluid_analysis
            self._fluid_phase_analysis = fluid_analysis.DEMCoupledFluidPlasmaDynamicsAnalysis(self.model, self.project_parameters, self.vars_man)
            self._fluid_phase_analysis.main_path = self.main_path 
        return self._fluid_phase_analysis
        

    def GetFieldUtility(self):
        return None
