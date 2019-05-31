# General comments to do: version of variables_management.py adapted to plasma dynamics application

from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidTransportApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *
import parameters_tools_for_plasma_dynamics as PT
def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)

class VariablesManager:
    @staticmethod
    def EliminateRepeatedValuesFromList(redundant_list):
        clean_list = []

        for var in redundant_list:

            if var in clean_list:
                redundant_list.remove(var)

            clean_list += [var]

    @staticmethod
    def AddNodalVariables(model_part, variable_list):

        for var in variable_list:
            model_part.AddNodalSolutionStepVariable(var)

    def __init__(self, parameters):
        self.project_parameters = parameters

    # constructing lists of variables to add
    # * Performing modifications to the input parameters for consistency (provisional until interface does it)
    # * Choosing the variables to be printed
    # * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
    #       instance to be filled with the other phase's info through the coupling process
    # * Listing nodal variables to be added to the model parts (memory will be allocated for them).
    #       Note that additional variables may be added as well by the fluid and/or DEM strategies.
    @staticmethod
    def AddFrameOfReferenceRelatedVariables(parameters, model_part):
        frame_of_reference_type = parameters["frame_of_reference_type"].GetInt()
        model_part.ProcessInfo.SetValue(FRAME_OF_REFERENCE_TYPE, frame_of_reference_type)

        if frame_of_reference_type == 1: # Rotating frame
            angular_velocity_of_frame = Vector(3)
            angular_velocity_of_frame[:] = [parameters["angular_velocity_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]

            model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME, angular_velocity_of_frame)

            if frame_of_reference_type >= 2: # General frame
                angular_velocity_of_frame_old = Vector(3)
                angular_velocity_of_frame_old[:] = [parameters["angular_velocity_of_frame_old" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                acceleration_of_frame_origin = Vector(3)
                acceleration_of_frame_origin[:] = [parameters["acceleration_of_frame_origin" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                angular_acceleration_of_frame = Vector(3)
                angular_acceleration_of_frame[:] = [parameters["angular_acceleration_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME_OLD, angular_velocity_of_frame_old)
                model_part.ProcessInfo.SetValue(ACCELERATION_MOVING_FRAME_ORIGIN, acceleration_of_frame_origin)
                model_part.ProcessInfo.SetValue(ANGULAR_ACCELERATION_MOVING_FRAME, angular_acceleration_of_frame)

    def AddExtraProcessInfoVariablesToFluidModelPart(self, parameters, fluid_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, fluid_model_part)

        fluid_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)
        #fluid_model_part.ProcessInfo.SetValue(ELECTRIC_POTENTIAL, 0.0)
        #fluid_model_part.ProcessInfo.SetValue(FLUID_ION_DENSITY, 1.0) TODO:remove this if not necessary

        fluid_electron_density = parameters["properties"][2]["density_parameters"]["fluid_electron_density"].GetDouble()
        fluid_neutral_density = parameters["properties"][2]["density_parameters"]["fluid_neutral_density"].GetDouble()
        fluid_model_part.ProcessInfo.SetValue(FLUID_ELECTRON_DENSITY, fluid_electron_density)
        fluid_model_part.ProcessInfo.SetValue(FLUID_NEUTRAL_DENSITY, fluid_neutral_density)

        electric_field = Vector(3)
        magnetic_field = Vector(3)
        fluid_model_part.ProcessInfo.SetValue(ELECTRIC_FIELD, electric_field)
        fluid_model_part.ProcessInfo.SetValue(MAGNETIC_FIELD, magnetic_field)        


    def AddExtraProcessInfoVariablesToDispersePhaseModelPart(self, parameters, dem_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, dem_model_part)
        dem_model_part.ProcessInfo.SetValue(COUPLING_TYPE, parameters["coupling_level_type"].GetInt())

        dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, parameters["fluid_model_type"].GetInt())


        external_electric_field = Vector(3)
        external_electric_field[0] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_electric_field_X"].GetDouble()
        external_electric_field[1] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_electric_field_Y"].GetDouble()
        external_electric_field[2] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_electric_field_Z"].GetDouble()

        external_magnetic_field = Vector(3)
        external_magnetic_field[0] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_magnetic_field_X"].GetDouble()
        external_magnetic_field[1] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_magnetic_field_Y"].GetDouble()  
        external_magnetic_field[2] = parameters["properties"][1]["electromagnetic_field_parameters"]["external_magnetic_field_Z"].GetDouble()

        dem_model_part.ProcessInfo.SetValue(EXTERNAL_ELECTRIC_FIELD, external_electric_field)
        dem_model_part.ProcessInfo.SetValue(EXTERNAL_MAGNETIC_FIELD, external_magnetic_field)

        #electric_field_projected_to_particle = Vector(3)
        #particle_ion_velocity = Vector(3)
        #dem_model_part.ProcessInfo.SetValue(ELECTRIC_FIELD_PROJECTED_TO_PARTICLE, electric_field_projected_to_particle)
        #dem_model_part.ProcessInfo.SetValue(PARTICLE_ION_VELOCITY, particle_ion_velocity)

        number_of_particles_in_a_ion_macroparticle = parameters["properties"][2]["density_parameters"]["number_of_particles_in_a_ion_macroparticle"].GetInt()
        dem_model_part.ProcessInfo.SetValue(NUMBER_OF_PARTICLES_IN_A_ION_MACROPARTICLE, number_of_particles_in_a_ion_macroparticle)





    def ConstructListsOfVariables(self, parameters):
        # PRINTING VARIABLES
        # constructing lists of variables to be printed
        self.gauss_points_results = []

        if self.project_parameters.Has('plasma_dynamics_output_processes'):
            gid_output_options = self.project_parameters["plasma_dynamics_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            gauss_point_results = result_file_configuration["gauss_point_results"]
            self.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]

        self.ConstructListsOfResultsToPrint(parameters)

        # COUPLING VARIABLES
        # listing the variables involved in the fluid-particles coupling

        if parameters["coupling"]["coupling_level_type"].GetInt():
            self.ConstructListsOfVariablesForCoupling(parameters)

        # VARIABLES TO ADD
        # listing nodal variables to be added to the model parts (memory will be allocated for them)

        # fluid variables
        self.fluid_vars = []
        self.fluid_vars += self.fluid_printing_vars
        self.fluid_vars += self.coupling_fluid_vars

        #TODO: remove if not necessary
        if (parameters["pressure_grad_recovery_type"].GetInt() > 1
            or parameters["material_acceleration_calculation_type"].GetInt() == 7
            or parameters["laplacian_calculation_type"].GetInt() > 1):
            self.fluid_vars += [Kratos.NODAL_WEIGHTS]

        # dem variables
        self.dem_vars = []
        self.dem_vars += self.dem_printing_vars
        self.dem_vars += self.coupling_dem_vars
        self.dem_vars += [VELOCITY_OLD]

        self.dem_vars += [DISPLACEMENT_OLD]
        self.dem_vars += [VELOCITY_OLD_OLD]


        if (parameters["custom_dem"]["translational_integration_scheme"].GetString()
            in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}):
            self.dem_vars += [VELOCITY_OLD]
            self.dem_vars += [ADDITIONAL_FORCE_OLD]
            self.dem_vars += [AUX_VEL]


        self.dem_vars += [PARTICLE_SPHERICITY] # TODO: add only when needed

        # clusters variables
        self.clusters_vars = []

        # rigid faces variables
        self.rigid_faces_vars = [VELOCITY,
                                 ANGULAR_VELOCITY,
                                 DISPLACEMENT,
                                 DELTA_DISPLACEMENT,
                                 DELTA_ROTATION,
                                 CONTACT_FORCES,
                                 DEM_PRESSURE,
                                 ELASTIC_FORCES,
                                 PRESSURE,
                                 TANGENTIAL_ELASTIC_FORCES,
                                 SHEAR_STRESS,
                                 NODAL_AREA,
                                 VELOCITY_OLD]

        if parameters["custom_fluid"]["embedded_option"].GetBool():
            self.rigid_faces_vars += [FORCE]
            self.rigid_faces_vars += [POSITIVE_FACE_PRESSURE]
            self.rigid_faces_vars += [NEGATIVE_FACE_PRESSURE]

        self.fluid_vars += self.rigid_faces_vars

        # inlet variables
        self.inlet_vars = self.dem_vars

    def ConstructListsOfResultsToPrint(self, parameters):

        #Construction of the list of DEM parameters to print
        dem_list = self.project_parameters["dem_nodal_results"]
        self.dem_nodal_results = [key for key in dem_list.keys() if dem_list[key].GetBool()]

        self.clusters_nodal_results = []
        self.rigid_faces_nodal_results = []

        #Construction of the list of fluid parameters to print
        fluid_list = self.project_parameters["fluid_nodal_results"]
        self.fluid_nodal_results = [key for key in fluid_list.keys() if fluid_list[key].GetBool()]


        #Adding results to print in case some parameters are used (redundancy is managed with EliminateRepeatedValuesFromList)
        if parameters["ElementType"].GetString() == "IonParticle3D":
            self.dem_nodal_results += ["EXTERNAL_APPLIED_FORCE"]
            self.dem_nodal_results += ["ELECTRIC_FIELD_PROJECTED_TO_PARTICLE"]
            self.dem_nodal_results += ["PARTICLE_ION_VELOCITY"]
            self.dem_nodal_results += ["MACROPARTICLE_ION_DENSITY"]

            self.fluid_nodal_results += ["ELECTRIC_POTENTIAL"]
            self.fluid_nodal_results += ["FLUID_ION_DENSITY"]
            self.fluid_nodal_results += ["FLUID_ELECTRON_DENSITY"]
            self.fluid_nodal_results += ["ELECTRIC_FIELD"]
            self.fluid_nodal_results += ["NODAL_PHI_GRADIENT"]

        self.ChangeListOfFluidNodalResultsToPrint(parameters)

        #Construction of the list of mixed parameters to print
        self.mixed_nodal_results = ["VELOCITY", "DISPLACEMENT"]

        #Construction of a special list for variables to print in a file
        self.variables_to_print_in_file = ["VELOCITY","ELECTRIC_FIELD_PROJECTED_TO_PARTICLE","EXTERNAL_APPLIED_FORCE"]


        #Initialization of the lists of evaluation of the nodal variables to print
        self.dem_printing_vars = []
        self.fluid_printing_vars = []

        self.clusters_printing_vars = []
        self.rigid_faces_printing_vars = []

        self.time_filtered_vars = []

        VariablesManager.EliminateRepeatedValuesFromList(self.fluid_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.dem_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.mixed_nodal_results)

        #Evaluation of the variables we want to print

        for variable in self.mixed_nodal_results:
            self.dem_printing_vars += [eval(variable)]
            self.fluid_printing_vars += [eval(variable)]

        for var in self.mixed_nodal_results:

            if var in self.fluid_nodal_results:
                self.fluid_nodal_results.remove(var)     
            if var in self.dem_nodal_results:
                self.dem_nodal_results.remove(var)                 


        for variable in self.fluid_nodal_results:
            self.fluid_printing_vars += [eval(variable)]

        for variable in self.dem_nodal_results:
            self.dem_printing_vars += [eval(variable)]

        for variable in self.clusters_nodal_results:
            self.clusters_printing_vars += [eval(variable)]

        for variable in self.rigid_faces_nodal_results:
            self.rigid_faces_printing_vars += [eval(variable)]




    def ConstructListsOfVariablesForCoupling(self, parameters):

        ## fluid coupling variables
        self.coupling_fluid_vars = []

        #for safety for the moment
        self.coupling_fluid_vars += [MATERIAL_ACCELERATION]

        self.coupling_fluid_vars += [FLUID_FRACTION]
        self.coupling_fluid_vars += [FLUID_FRACTION_OLD]

        self.coupling_fluid_vars += [DISPERSE_FRACTION]

        self.coupling_fluid_vars += [PARTICLE_VEL_FILTERED]
        self.coupling_fluid_vars += [TIME_AVERAGED_ARRAY_3]
        self.coupling_fluid_vars += [PHASE_FRACTION]

        self.coupling_fluid_vars += [FLUID_FRACTION_GRADIENT]
        self.coupling_fluid_vars += [FLUID_FRACTION_RATE]


        self.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

        #Plasma Dynamics coupling variables: forward coupling
        self.coupling_fluid_vars += [ELECTRIC_POTENTIAL]
        self.coupling_fluid_vars += [ELECTRIC_FIELD]
        self.coupling_fluid_vars += [FLUID_ELECTRON_DENSITY]
        self.coupling_fluid_vars += [FLUID_NEUTRAL_DENSITY]
        
        #Plasma Dynamics coupling variables: backward coupling
        self.coupling_fluid_vars += [FLUID_ION_DENSITY]



        ## dem coupling variables
        self.coupling_dem_vars = []

        #for safety for the moment
        self.coupling_dem_vars += [FLUID_VEL_PROJECTED]
        self.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]
        self.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]

        self.coupling_dem_vars += [MATERIAL_FLUID_ACCEL_PROJECTED]
        self.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]
        self.coupling_dem_vars += [FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED]
        self.coupling_dem_vars += [ADDITIONAL_FORCE]  # Here for safety for the moment

        self.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]

        #Plasma Dynamics coupling variables: forward coupling
        self.coupling_dem_vars += [ELECTRIC_FIELD_PROJECTED_TO_PARTICLE]

        #Plasma Dynamics coupling variables: backward coupling
        #self.coupling_dem_vars += [MACROPARTICLE_ION_DENSITY]


        if parameters["coupling"]["backward_coupling"]["apply_time_filter_to_fluid_fraction_option"].GetBool():
            self.time_filtered_vars += [FLUID_FRACTION_FILTERED]

        if parameters["coupling"]["backward_coupling"]["filter_velocity_option"].GetBool():
            self.time_filtered_vars += [PARTICLE_VEL_FILTERED]

    def ChangeListOfFluidNodalResultsToPrint(self, parameters):
        fluid_list = self.project_parameters["fluid_nodal_results"]
        self.fluid_nodal_results.extend(key for key in fluid_list.keys() if fluid_list[key].GetBool())