import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import json

class FluidDynamicsAnalysisROM(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters, build_petrov_galerkin=False, solve_petrov_galerkin=False,solve_least_squares=False):
        self.build_petrov_galerkin = build_petrov_galerkin
        self.solve_petrov_galerkin = solve_petrov_galerkin
        self.solve_least_squares = solve_least_squares
        super().__init__(model,project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        with open('RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
            if self.build_petrov_galerkin:
                self.project_parameters["solver_settings"].AddBool("build_petrov_galerkin", self.build_petrov_galerkin)
            if self.solve_petrov_galerkin:
                self.project_parameters["solver_settings"].AddBool("solve_petrov_galerkin", self.solve_petrov_galerkin)
                self.project_parameters["solver_settings"].AddValue("rom_residual_settings",rom_settings["Petrov_Galerkin_basis"]["rom_settings"])
            if self.solve_least_squares:
                self.project_parameters["solver_settings"].AddBool("solve_least_squares", self.solve_least_squares)
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyAfterSolverInitialize(self):
        """Here is where the ROM_BASIS is imposed to each node"""
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            if self.solve_petrov_galerkin:
                nodal_residual_modes = data["Petrov_Galerkin_basis"]["nodal_modes"]####ADDEDPETROV
                nodal_residual_dofs = len(data["Petrov_Galerkin_basis"]["rom_settings"]["nodal_unknowns"])####ADDEDPETROV
                rom_residual_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()####ADDEDPETROV
                for node in computing_model_part.Nodes:
                    aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                    aux_residual = KratosMultiphysics.Matrix(nodal_residual_dofs,rom_residual_dofs)####ADDEDPETROV
                    for j in range(nodal_dofs):
                        Counter=str(node.Id)
                        for i in range(rom_residual_dofs):
                            if (i<rom_dofs):
                                aux[j,i] = nodal_modes[Counter][j][i]
                            aux_residual[j,i] = nodal_residual_modes[Counter][j][i]####ADDEDPETROV
                    node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                    node.SetValue(romapp.ROM_BASIS_ASSEMBLED_RESIDUALS, aux_residual)####ADDEDPETROV
                    counter+=1
            else:
                for node in computing_model_part.Nodes:
                    aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                    for j in range(nodal_dofs):
                        Counter=str(node.Id)
                        for i in range(rom_dofs):
                            aux[j,i] = nodal_modes[Counter][j][i]
                    node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                    counter+=1
