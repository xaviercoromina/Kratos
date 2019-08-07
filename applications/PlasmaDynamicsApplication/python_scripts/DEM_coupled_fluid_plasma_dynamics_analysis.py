from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.FluidTransportApplication

from fluid_transport_analysis import FluidTransportAnalysis

class DEMCoupledFluidPlasmaDynamicsAnalysis(FluidTransportAnalysis):

    def __init__(self, model, parameters=None, variables_management=None):
        self.model = model
        self.plasma_dynamics_project_parameters = parameters
        self.project_parameters = self.plasma_dynamics_project_parameters['fluid_parameters']
        self.vars_man = variables_management

        super(DEMCoupledFluidPlasmaDynamicsAnalysis, self).__init__(model, self.project_parameters)
        self.fluid_model_part = self._GetSolver().main_model_part

    def Initialize(self):
        self.AddFluidVariablesByPlasmaDynamicsAlgorithm()
        super(DEMCoupledFluidPlasmaDynamicsAnalysis, self).Initialize()

    def AddFluidVariablesByPlasmaDynamicsAlgorithm(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)

    def CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning
        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep" for fluid'
                IssueDeprecationWarning("PlasmaAnalysis", warn_msg)

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Parameters(parameter_file.read())

    model = Model()
    simulation = DEMCoupledFluidPlasmaDynamicsAnalysis(model, parameters)
    simulation.Run()
