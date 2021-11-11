from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_static_solver import StaticMechanicalSolver
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(main_model_part, custom_settings):
    return ROMSolver(main_model_part, custom_settings)

class ROMSolver(StaticMechanicalSolver):
    """The stationary class for ROM structural mechanics solvers.

    See structural_mechanics_static_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        super(ROMSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "rom_settings": {
            "nodal_unknowns": [ "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "number_of_rom_dofs": 3
            },
            "build_petrov_galerkin": false,
            "solve_petrov_galerkin": false,
            "solve_least_squares" : false,
            "rom_residual_settings": {
            "nodal_unknowns": [ "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "number_of_rom_dofs": 0
            }
        }
        """)#ADDEDPETROV
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        super(ROMSolver, self).AddVariables() #Adding nodal area variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        rom_parameters=self.settings["rom_settings"]
        if self.settings["build_petrov_galerkin"].GetBool():
            rom_parameters.AddBool("build_petrov_galerkin",True)
        elif self.settings["solve_petrov_galerkin"].GetBool():
            rom_parameters.AddBool("solve_petrov_galerkin", True)# ADDEDPETROV
            rom_parameters.AddInt("number_of_rom_residual_dofs", self.settings["rom_residual_settings"]["number_of_rom_dofs"].GetInt())# ADDEDPETROV
            rom_parameters.AddEmptyList("nodal_residual_unknowns") #ADDEDPETROV
            rom_parameters["nodal_residual_unknowns"].SetStringArray(self.settings["rom_residual_settings"]["nodal_unknowns"].GetStringArray())#ADDEDPETROV
        elif self.settings["solve_least_squares"].GetBool():
            rom_parameters.AddBool("solve_least_squares",True)
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
