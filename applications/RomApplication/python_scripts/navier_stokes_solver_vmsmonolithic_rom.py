# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(model, custom_settings):
    return ROMSolver(model, custom_settings)

class ROMSolver(NavierStokesSolverMonolithic):

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "rom_settings": {
            "nodal_unknowns": [ "CFD_DOFS_USED_LISTED_HERE" ],
            "number_of_rom_dofs": 3
            },
            "build_petrov_galerkin": false,
            "solve_petrov_galerkin": false,
            "rom_residual_settings": {
            "nodal_unknowns": [ "RESIDUAL_X", "RESIDUAL_Y", "RESIDUAL_P"],
            "number_of_rom_dofs": 0
            }
        }
        """)#ADDEDPETROV
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultParameters())
        return default_settings

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        rom_parameters=self.settings["rom_settings"]
        if self.settings["build_petrov_galerkin"].GetBool():
            rom_parameters.AddBool("build_petrov_galerkin",True)
        elif self.settings["solve_petrov_galerkin"].GetBool():
            rom_parameters.AddBool("solve_petrov_galerkin", True)# ADDEDPETROV
            rom_parameters.AddInt("number_of_rom_residual_dofs", self.settings["rom_residual_settings"]["number_of_rom_dofs"].GetInt())# ADDEDPETROV
            rom_parameters.AddEmptyList("nodal_residual_unknowns") #ADDEDPETROV
            rom_parameters["nodal_residual_unknowns"].SetStringArray(self.settings["rom_residual_settings"]["nodal_unknowns"].GetStringArray())#ADDEDPETROV
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
