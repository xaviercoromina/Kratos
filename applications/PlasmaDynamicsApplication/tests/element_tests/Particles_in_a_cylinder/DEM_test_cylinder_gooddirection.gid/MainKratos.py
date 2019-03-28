import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
import KratosMultiphysics.DEMApplication
from KratosMultiphysics.FluidTransportApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *
import main_script as Main

model = KratosMultiphysics.Model()
solution = Main.Solution(model)
solution.Run()
