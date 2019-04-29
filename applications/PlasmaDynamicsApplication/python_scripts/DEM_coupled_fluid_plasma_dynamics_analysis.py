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

    def RunSingleTimeStep(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

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
