import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import matplotlib.pyplot as plt

from math import pi
from numpy import sin, cos, arctan, exp, arctan2
import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarManufacturedConstraintFunctionProcess3(Model, settings["Parameters"])


class ApplyScalarManufacturedConstraintFunctionProcess3(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # self.a1 = 2
        # self.a2 = 3
        self.k = 0.000001
        self.s = 1000.0
        self.rms = 0.0

        self.times = []
        self.error_list = []
        self.rms_list = []

    def source_term(self, x, y, z, t):

        #f = 4.0*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 6.0*x*y*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 16*x*y*self.s*(1 - x)*(1 - y)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + x*y*pi*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*cos(pi*t) + 3.0*x*y*(16*x - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 2.0*x*y*(16*y - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 3.0*x*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 2.0*y*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) - self.k*(-8*x*y*self.k**(-1.5)*(1.0 - 2*x)*(1 - y)*(2.0 - 4*x)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 8*x*y*self.k**(-1.5)*(1.0 - 2*y)*(1 - y)*(2.0 - 4*y)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 32*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*x)*(16*y - 16)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16*x - 16)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 8*x*y*self.k**(-0.5)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 4*x*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*(16*x - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 4*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*y*(16*y - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t))

        f = (x*y*self.k**0.5*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(-128.0*x - 192.0*y + 160.0)*sin(pi*t) - 16*x*y*pi*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*cos(pi*t) + 32*self.k*(-4*x*y*self.k**0.5*(x - 1)*(y - 1)*((2*x - 1.0)*(4*x - 2.0) + (2*y - 1.0)*(4*y - 2.0))*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625) + 2*self.k**0.5*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(2*x*y*(x - 1)*(y - 1) + x*y*(x - 1)*(2*y - 1.0) + x*y*(2*x - 1.0)*(y - 1) + x*(x - 1)*(y - 1)*(2*y - 1.0) + y*(x - 1)*(2*x - 1.0)*(y - 1)) + (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(x*(x - 1) + y*(y - 1))*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi))*sin(pi*t) - (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*(16*x*y*self.s*(x - 1)*(y - 1) + 48.0*x*y*(x - 1) + 32.0*x*y*(y - 1) + 48.0*x*(x - 1)*(y - 1) + 32.0*y*(x - 1)*(y - 1))*sin(pi*t))/(pi*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2)

        if np.sqrt(((x-1/4)**2 + (y-1/2)**2)) < 0.001:
            x = x+0.001
            y = y-0.001
            f = (x*y*self.k**0.5*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(-128.0*x - 192.0*y + 160.0)*sin(pi*t) - 16*x*y*pi*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*cos(pi*t) + 32*self.k*(-4*x*y*self.k**0.5*(x - 1)*(y - 1)*((2*x - 1.0)*(4*x - 2.0) + (2*y - 1.0)*(4*y - 2.0))*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625) + 2*self.k**0.5*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(2*x*y*(x - 1)*(y - 1) + x*y*(x - 1)*(2*y - 1.0) + x*y*(2*x - 1.0)*(y - 1) + x*(x - 1)*(y - 1)*(2*y - 1.0) + y*(x - 1)*(2*x - 1.0)*(y - 1)) + (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(x*(x - 1) + y*(y - 1))*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi))*sin(pi*t) - (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*(16*x*y*self.s*(x - 1)*(y - 1) + 48.0*x*y*(x - 1) + 32.0*x*y*(y - 1) + 48.0*x*(x - 1)*(y - 1) + 32.0*y*(x - 1)*(y - 1))*sin(pi*t))/(pi*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2)

        if np.sqrt(((x-3/4)**2 + (y-1/2)**2)) < 0.001:
            x = x-0.001
            y = y-0.001
            f = (x*y*self.k**0.5*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(-128.0*x - 192.0*y + 160.0)*sin(pi*t) - 16*x*y*pi*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*cos(pi*t) + 32*self.k*(-4*x*y*self.k**0.5*(x - 1)*(y - 1)*((2*x - 1.0)*(4*x - 2.0) + (2*y - 1.0)*(4*y - 2.0))*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625) + 2*self.k**0.5*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(2*x*y*(x - 1)*(y - 1) + x*y*(x - 1)*(2*y - 1.0) + x*y*(2*x - 1.0)*(y - 1) + x*(x - 1)*(y - 1)*(2*y - 1.0) + y*(x - 1)*(2*x - 1.0)*(y - 1)) + (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(x*(x - 1) + y*(y - 1))*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi))*sin(pi*t) - (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*(16*x*y*self.s*(x - 1)*(y - 1) + 48.0*x*y*(x - 1) + 32.0*x*y*(y - 1) + 48.0*x*(x - 1)*(y - 1) + 32.0*y*(x - 1)*(y - 1))*sin(pi*t))/(pi*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2)

        if np.sqrt(((x-1/2)**2 + (y-1/4)**2)) < 0.001:
            x = x+0.001
            y = y-0.001
            f = (x*y*self.k**0.5*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(-128.0*x - 192.0*y + 160.0)*sin(pi*t) - 16*x*y*pi*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*cos(pi*t) + 32*self.k*(-4*x*y*self.k**0.5*(x - 1)*(y - 1)*((2*x - 1.0)*(4*x - 2.0) + (2*y - 1.0)*(4*y - 2.0))*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625) + 2*self.k**0.5*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(2*x*y*(x - 1)*(y - 1) + x*y*(x - 1)*(2*y - 1.0) + x*y*(2*x - 1.0)*(y - 1) + x*(x - 1)*(y - 1)*(2*y - 1.0) + y*(x - 1)*(2*x - 1.0)*(y - 1)) + (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(x*(x - 1) + y*(y - 1))*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi))*sin(pi*t) - (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*(16*x*y*self.s*(x - 1)*(y - 1) + 48.0*x*y*(x - 1) + 32.0*x*y*(y - 1) + 48.0*x*(x - 1)*(y - 1) + 32.0*y*(x - 1)*(y - 1))*sin(pi*t))/(pi*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2)

        if np.sqrt(((x-1/2)**2 + (y-3/4)**2)) < 0.001:
            x = x+0.001
            y = y+0.001
            f = (x*y*self.k**0.5*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(-128.0*x - 192.0*y + 160.0)*sin(pi*t) - 16*x*y*pi*(x - 1)*(y - 1)*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*cos(pi*t) + 32*self.k*(-4*x*y*self.k**0.5*(x - 1)*(y - 1)*((2*x - 1.0)*(4*x - 2.0) + (2*y - 1.0)*(4*y - 2.0))*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625) + 2*self.k**0.5*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)*(2*x*y*(x - 1)*(y - 1) + x*y*(x - 1)*(2*y - 1.0) + x*y*(2*x - 1.0)*(y - 1) + x*(x - 1)*(y - 1)*(2*y - 1.0) + y*(x - 1)*(2*x - 1.0)*(y - 1)) + (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(x*(x - 1) + y*(y - 1))*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi))*sin(pi*t) - (self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2*(arctan2(self.k**(-0.5)*(2*(x - 0.5)**2 + 2*(y - 0.5)**2 - 0.125),1) - 0.5*pi)*(16*x*y*self.s*(x - 1)*(y - 1) + 48.0*x*y*(x - 1) + 32.0*x*y*(y - 1) + 48.0*x*(x - 1)*(y - 1) + 32.0*y*(x - 1)*(y - 1))*sin(pi*t))/(pi*(self.k**1.0 + 4*((x - 0.5)**2 + (y - 0.5)**2 - 0.0625)**2)**2)
        return f

   # def velocity(self, x, y ,z ,t):
   #     return  # TODO: paste here the function

    def analytical_solution(self, x, y, z, t):
        u = 16 * sin(pi * t) * x * (1 - x) * y*(1 - y) * (1/2 + (arctan2(2 * self.k ** (-1/2) * (0.25**2 - (x - 0.5)**2 - (y - 0.5)**2),1))/pi)

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

        print ("***************")
        print ("RMS = ",self.rms)
        print ("***************")