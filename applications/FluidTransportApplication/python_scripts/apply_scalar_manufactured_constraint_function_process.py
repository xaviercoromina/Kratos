import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

from math import pi
from numpy import sin, cos, arctan, exp
import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarManufacturedConstraintFunctionProcess(Model, settings["Parameters"])


class ApplyScalarManufacturedConstraintFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # self.a1 = 2
        # self.a2 = 3
        self.k = 0.000001
        self.s = 500.0

    def source_term(self, x, y, z, t):

        f = 4.0*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 6.0*x*y*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 16*x*y*self.s*(1 - x)*(1 - y)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + x*y*pi*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*cos(pi*t) + 3.0*x*y*(16*x - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 2.0*x*y*(16*y - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 3.0*x*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 2.0*y*(1 - y)*(16 - 16*x)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) - self.k*(-8*x*y*self.k**(-1.5)*(1.0 - 2*x)*(1 - y)*(2.0 - 4*x)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 8*x*y*self.k**(-1.5)*(1.0 - 2*y)*(1 - y)*(2.0 - 4*y)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 32*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*x)*(16*y - 16)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16*x - 16)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 8*x*y*self.k**(-0.5)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 4*x*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*(16*x - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t) + 4*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)*sin(pi*t)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*y*(16*y - 16)*(arctan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*t))

        # sense terme temporal
        # f = 4.0*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 6.0*x*y*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 16*x*y*self.s*(1 - x)*(1 - y)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) + 3.0*x*y*(16*x - 16)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) + 2.0*x*y*(16*y - 16)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) + 3.0*x*(1 - y)*(16 - 16*x)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) + 2.0*y*(1 - y)*(16 - 16*x)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) - self.k*(-8*x*y*self.k**(-1.5)*(1.0 - 2*x)*(1 - y)*(2.0 - 4*x)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 8*x*y*self.k**(-1.5)*(1.0 - 2*y)*(1 - y)*(2.0 - 4*y)*(16 - 16*x)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)**2) - 32*x*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*x)*(16*y - 16)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*y*self.k**(-0.5)*(1.0 - 2*y)*(16*x - 16)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) - 8*x*y*self.k**(-0.5)*(1 - y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 4*x*self.k**(-0.5)*(1.0 - 2*y)*(1 - y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*x*(16*x - 16)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5) + 4*y*self.k**(-0.5)*(1.0 - 2*x)*(1 - y)*(16 - 16*x)/(pi*(4*self.k**(-1.0)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625)**2 + 1)) + 2*y*(16*y - 16)*(atan(2*self.k**(-0.5)*(-(x - 0.5)**2 - (y - 0.5)**2 + 0.0625))/pi + 0.5))

        return f

   # def velocity(self, x, y ,z ,t):
   #     return  # TODO: paste here the function

    def analytical_solution(self, x, y, z, t):
        u = 16 * sin(pi * t) * x * (1 - x) * y*(1 - y) * (1/2 + (arctan(2 * self.k ** (-1/2) * (0.25**2 - (x - 0.5)**2 - (y - 0.5)**2)))/pi)

        # sense terme temporal:
        # u = 16 * x * (1 - x) * y*(1 - y) * (1/2 + (arctan(2 * self.k ** (-1/2) * (0.25**2 - (x - 0.5)**2 - (y - 0.5)**2)))/pi)

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
