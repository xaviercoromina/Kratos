import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
from KratosMultiphysics.FluidTransportApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *
import main_script as Main

model = KratosMultiphysics.Model()
solution = Main.Solution(model)
solution.Run()
