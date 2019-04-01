# General comments to do: version of variables_management.py adapted to plasma dynamics application

from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidTransportApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *
import ast
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

            if frame_of_reference_type >= 2: # Gemeral frame
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
        gravity = Vector(3)
        if parameters["body_force_on_fluid_option"].GetBool():
            gravity[0] = parameters["GravityX"].GetDouble()
            gravity[1] = parameters["GravityY"].GetDouble()
            gravity[2] = parameters["GravityZ"].GetDouble()
        fluid_model_part.ProcessInfo.SetValue(GRAVITY, gravity)


    def AddExtraProcessInfoVariablesToDispersePhaseModelPart(self, parameters, dem_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, dem_model_part)
        dem_model_part.ProcessInfo.SetValue(COUPLING_TYPE, parameters["coupling_level_type"].GetInt())

        dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, parameters["fluid_model_type"].GetInt())

        for prop in parameters["properties"].values():
            if prop["plasma_dynamics_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["plasma_dynamics_law_parameters"]["history_force_parameters"]

                if (prop["plasma_dynamics_law_parameters"]["name"].GetString() != "default"
                    and history_force_parameters["name"].GetString() != "default"):
                    break


    def ConstructListsOfVariables(self, parameters):
        # PRINTING VARIABLES
        # constructing lists of variables to be printed
        self.ConstructListsOfResultsToPrint(parameters)

        # COUPLING VARIABLES
        # listing the variables involved in the fluid-particles coupling

        if parameters["coupling_level_type"].GetInt():
            self.ConstructListsOfVariablesForCoupling(parameters)

        # VARIABLES TO ADD
        # listing nodal variables to be added to the model parts (memory will be allocated for them)

        # fluid variables
        self.fluid_vars = []
        self.fluid_vars += self.fluid_printing_vars
        #self.fluid_vars += self.coupling_fluid_vars



        #self.fluid_vars += [NODAL_WEIGHTS]



        # dem variables
        self.dem_vars = []
        self.dem_vars += self.dem_printing_vars
        #self.dem_vars += self.coupling_dem_vars
        self.dem_vars += [BUOYANCY]
        self.dem_vars += [VELOCITY_OLD]

        if parameters["frame_of_reference_type"].GetInt() and parameters["basset_force_type"].GetInt() > 0:
            self.dem_vars += [DISPLACEMENT_OLD]
            self.dem_vars += [VELOCITY_OLD_OLD]

        if (parameters["TranslationalIntegrationScheme"].GetString()
            in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}
            or parameters["basset_force_type"].GetInt() > 0):
            self.dem_vars += [VELOCITY_OLD]
            self.dem_vars += [ADDITIONAL_FORCE_OLD]
            self.dem_vars += [AUX_VEL]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [DRAG_FORCE]

        self.dem_vars += [PARTICLE_SPHERICITY] # TODO: add only when needed

        if (PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')
            and parameters["add_each_hydro_force_option"].GetBool()):
            self.dem_vars += [LIFT_FORCE]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [VIRTUAL_MASS_FORCE]

        will_need_basset_force_variable = False
        for prop in parameters["properties"].values():
            if prop["plasma_dynamics_law_parameters"].Has("history_force_parameters"):
                if prop["plasma_dynamics_law_parameters"]["history_force_parameters"]["name"].GetString() != 'default':
                    will_need_basset_force_variable = True
                    break

        if will_need_basset_force_variable:
            self.dem_vars += [BASSET_FORCE]

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

        if parameters["embedded_option"].GetBool():
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
            self.dem_nodal_results += ["MACROPARTICLE_ION_DENSITY"]
            self.dem_nodal_results += ["PARTICLE_ION_VELOCITY"]

            self.fluid_nodal_results += ["ELECTRIC_POTENTIAL"]
            self.fluid_nodal_results += ["FLUID_ION_DENSITY"]
            self.fluid_nodal_results += ["FLUID_ELECTRON_DENSITY"]
            self.fluid_nodal_results += ["FLUID_NEUTRAL_DENSITY"]
            self.fluid_nodal_results += ["ELECTRIC_FIELD"]
            self.fluid_nodal_results += ["MAGNETIC_FIELD"]










        # changes on the fluid variables to print for the sake of consistency
        #self.ChangeListOfFluidNodalResultsToPrint(parameters)

        #Construction of the list of mixed parameters to print
        self.mixed_nodal_results = []

        for var in self.mixed_nodal_results:

            if var in self.fluid_nodal_results:
                self.fluid_nodal_results.remove(var)

            if var in self.dem_nodal_results:
                self.dem_nodal_results.remove(var)   

        #Construction of a special list for variables to print in a file
        self.variables_to_print_in_file = ["VELOCITY"]



        #Evaluation of the variables we want to print
        self.dem_printing_vars = []
        self.fluid_printing_vars = []

        self.clusters_printing_vars = []
        self.rigid_faces_printing_vars = []

        self.time_filtered_vars = []

     

        for variable in self.fluid_nodal_results:
            self.fluid_printing_vars += [eval(variable)]

        for variable in self.dem_nodal_results:
            self.dem_printing_vars += [eval(variable)]

        for variable in self.clusters_nodal_results:
            self.clusters_printing_vars += [eval(variable)]

        for variable in self.rigid_faces_nodal_results:
            self.rigid_faces_printing_vars += [eval(variable)]

        for variable in self.mixed_nodal_results:
            self.dem_printing_vars += [eval(variable)]
            self.fluid_printing_vars += [eval(variable)]

        VariablesManager.EliminateRepeatedValuesFromList(self.fluid_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.dem_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.mixed_nodal_results)


    def ConstructListsOfVariablesForCoupling(self, parameters):

        # fluid coupling variables
        self.coupling_fluid_vars = []


        # dem coupling variables
        self.coupling_dem_vars = []

        if parameters["coupling_level_type"].GetInt() > 0:
            self.coupling_dem_vars += [ADDITIONAL_FORCE] # Here for safety for the moment


        if parameters["coupling_level_type"].GetInt() >= 1 or parameters["fluid_model_type"].GetInt() == 0:
            self.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]


        if parameters["filter_velocity_option"].GetBool():
            self.time_filtered_vars += [PARTICLE_VEL_FILTERED]


    def ChangeListOfFluidNodalResultsToPrint(self, parameters):

        if parameters["store_full_gradient_option"].GetBool() and 'VELOCITY_GRADIENT' in self.fluid_nodal_results:
            self.fluid_nodal_results += ["VELOCITY_X_GRADIENT"]
            self.fluid_nodal_results += ["VELOCITY_Y_GRADIENT"]
            self.fluid_nodal_results += ["VELOCITY_Z_GRADIENT"]

        if parameters["fluid_model_type"].GetInt() == 0 and 'AVERAGED_FLUID_VELOCITY' in self.fluid_nodal_results:
            self.fluid_nodal_results += ["AVERAGED_FLUID_VELOCITY"]

        if parameters["fluid_model_type"].GetInt() == 1 and 'FLUID_FRACTION_GRADIENT' in self.fluid_nodal_results:
            self.fluid_nodal_results += ["FLUID_FRACTION_GRADIENT"]