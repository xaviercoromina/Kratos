import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import matplotlib.pyplot as plt
import csv

from math import pi
from numpy import sin, cos, arctan, exp
import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarManufacturedConstraintFunctionProcess2(Model, settings["Parameters"])


class ApplyScalarManufacturedConstraintFunctionProcess2(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # self.a1 = 2
        # self.a2 = 3
        self.k = 0.01
        self.s = 10000.0
        self.rms = 0.0

        self.times = []
        self.error_list = []
        self.rms_list = []

    def source_term(self, x, y, z, t):

        f = -self.k*((1.15470053837925 - 0.75*(1.33333333333333*self.k + 0.577350269189626)*exp(-0.866025403784439*(1.0 - y)/self.k)/(self.k**2*(1.0 - exp(-0.866025403784439/self.k))))*(4.0*x*self.k + x**2.0 + (4.0*self.k + 1.0)*(-exp(-0.5*(1.0 - x)/self.k) + exp(-0.5/self.k))/(1 - exp(-0.5/self.k))) + (2.0 - 0.25*(4.0*self.k + 1.0)*exp(-0.5*(1.0 - x)/self.k)/(self.k**2*(1 - exp(-0.5/self.k))))*(1.33333333333333*y*self.k + 0.577350269189626*y**2.0 + (1.33333333333333*self.k + 0.577350269189626)*(-exp(-0.866025403784439*(1.0 - y)/self.k) + exp(-0.866025403784439/self.k))/(1.0 - exp(-0.866025403784439/self.k)))) + self.s*(4.0*x*self.k + x**2.0 + (4.0*self.k + 1.0)*(-exp(-0.5*(1.0 - x)/self.k) + exp(-0.5/self.k))/(1 - exp(-0.5/self.k)))*(1.33333333333333*y*self.k + 0.577350269189626*y**2.0 + (1.33333333333333*self.k + 0.577350269189626)*(-exp(-0.866025403784439*(1.0 - y)/self.k) + exp(-0.866025403784439/self.k))/(1.0 - exp(-0.866025403784439/self.k))) + 0.5*(2.0*x**1.0 + 4.0*self.k - 0.5*(4.0*self.k + 1.0)*exp(-0.5*(1.0 - x)/self.k)/(self.k*(1 - exp(-0.5/self.k))))*(1.33333333333333*y*self.k + 0.577350269189626*y**2.0 + (1.33333333333333*self.k + 0.577350269189626)*(-exp(-0.866025403784439*(1.0 - y)/self.k) + exp(-0.866025403784439/self.k))/(1.0 - exp(-0.866025403784439/self.k))) + 0.866025403784439*(1.15470053837925*y**1.0 + 1.33333333333333*self.k - 0.866025403784439*(1.33333333333333*self.k + 0.577350269189626)*exp(-0.866025403784439*(1.0 - y)/self.k)/(self.k*(1.0 - exp(-0.866025403784439/self.k))))*(4.0*x*self.k + x**2.0 + (4.0*self.k + 1.0)*(-exp(-0.5*(1.0 - x)/self.k) + exp(-0.5/self.k))/(1 - exp(-0.5/self.k)))

        return f

   # def velocity(self, x, y ,z ,t):
   #     return  # TODO: paste here the function

    def analytical_solution(self, x, y, z, t):

        u = (x ** 2.0 + 4.0 * self.k * x + (1.0 + 4.0 * self.k) * (exp(-1.0 /(2.0 * self.k)) - exp((-1.0 /(2.0 * self.k)) * (1.0 - x))) / (1 - exp(-1.0 /(2.0 * self.k)))) * (y ** 2.0 / (np.sqrt(3.0)) + 4.0 * self.k * y / 3.0 + (1.0 / (np.sqrt(3.0)) + 4.0 * self.k / 3.0) * (exp(- (np.sqrt(3.0)) /(2.0 * self.k)) - exp((- (np.sqrt(3.0)) /(2.0 * self.k)) * (1.0 - y))) / (1.0 - exp(- (np.sqrt(3.0)) /(2.0 * self.k))))

        return u

    def ExecuteInitializeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            f = self.source_term(node.X, node.Y, node.Z, time)
           # a = self.velocity(node.X, node.Y, node.Z, time)
            u = self.analytical_solution(node.X, node.Y, node.Z, time)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, f)
           # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, a)
            node.SetSolutionStepValue(KratosFluidTransport.NODAL_ANALYTIC_SOLUTION, u) # put the right variable

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        error = 0.0
        sum_error = 0.0
        sum_error_sqr = 0.0
        i = 0

        for node in self.model_part.Nodes:
            error = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.analytical_solution(node.X, node.Y, node.Z, time)
            sum_error += abs(error)
            sum_error_sqr += error ** 2
            i += 1

        #rms = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare(self.model_part, error)
        self.rms = np.sqrt(sum_error_sqr / i )

#        print ("***************")
#        print ("Error = ", sum_error)
#        print ("RMS = ",rms)
#        print ("***************")

        self.times.append (time)
        self.error_list.append (sum_error)
        self.rms_list.append (self.rms)

    def ExecuteFinalize(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        #f = plt.figure()

        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('error', color=color)
        ax1.plot(self.times, self.error_list, color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('rms', color=color)  # we already handled the x-label with ax1
        ax2.plot(self.times, self.rms_list, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        fig.savefig('plot_' + str(round(time, 2)) + '_delta_t_' + str(delta_time) + '_k_' + str(round(self.k, 2)) + '_s_' + str(round(self.s, 2)) + '.pdf', bbox_inches='tight')
        plt.show()

        np.savetxt('file_name_s'+ str(round(self.s, 2)) +'.csv', self.rms_list, delimiter=",", fmt='%s')


        print ("***************")
        print ("RMS = ",self.rms)
        print ("***************")