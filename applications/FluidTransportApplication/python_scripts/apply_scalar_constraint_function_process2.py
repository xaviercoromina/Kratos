import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import math
import matplotlib.pyplot as plt

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarConstraintFunctionProcess2(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyScalarConstraintFunctionProcess2(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []
        self.times = []
        self.suma_list = []

    def ExecuteFinalizeSolutionStep(self):
        self.suma = 0
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        for node in self.model_part.Nodes:
            self.suma += node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)/3.0

        print ("***************")
        print ("Suma = ",self.suma)
        print ("***************")
        self.times.append (time)
        self.suma_list.append (self.suma)

    def ExecuteFinalize(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        #f = plt.figure()
        k = 0.001

        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Suma nodes', color=color)
        ax1.plot(self.times, self.suma_list, color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        fig.savefig('plot_' + 'k_' + str(k) + '.pdf', bbox_inches='tight')
        plt.show()
        # for node in model_part.Nodes:
        #     radius = math.sqrt((node.X - 1.0) ** 2.0 + (node.Y - 1.0) ** 2.0)

        #     if(radius <= 0.6):

        # #     #if((node.X >= 0.1 and node.X <= 0.2 and node.Y >= 0.4 and node.Y <= 0.6) or (node.X >= 0.3 and node.X <= 0.4 and node.Y >= 0.4 and node.Y <= 0.6)):

        # #     if((node.X >= 0.1 and node.X <= 0.2)): #or (node.X >= 0.3 and node.X <= 0.4)):
        #         node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,1000.0)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,1000.0)
        #         node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,1,1000.0)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,1,1000.0)
        #         # node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,2,1.0)
        #         # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,2,1.0)
        #     else:
        #         node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,0.0)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0.0)
        #         node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,1,0.0)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,1,0.0)
        #         # node.SetSolutionStepValue(KratosFluidTransport.PHI_THETA,2,0.0)
        #         # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,2,0.0)


    # def ExecuteInitialize(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitialize()

    # def ExecuteInitializeSolutionStep(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitializeSolutionStep()