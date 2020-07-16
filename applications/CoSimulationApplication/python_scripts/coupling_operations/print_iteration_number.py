from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

class PrintIterationNumberOperation(CoSimulationCouplingOperation):
    # this is a dummy implementation, this be used with iterative coupled-solvers
    def Initialize(self):
        pass

    def Finalize(self):
        pass


    def InitializeSolutionStep(self):
        self.num_coupling_operations = 0

    def FinalizeSolutionStep(self):
        print("the number of coupling operations in this step was:", self.num_coupling_operations)


    def InitializeCouplingIteration(self):
        self.num_coupling_operations += 1

    def FinalizeCouplingIteration(self):
        pass


    def Execute(self):
        pass


    def PrintInfo(self):
        pass

    def Check(self):
        pass

