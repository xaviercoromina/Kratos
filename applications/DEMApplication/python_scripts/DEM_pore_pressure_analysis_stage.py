import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import math

class DEMPorePressureAnalysisStage(DEMAnalysisStage):

    def InitializeSolutionStep(self):
        
        super().InitializeSolutionStep()
        
        # TODO: Well data, depends on the particular GiD project
        # It is assumed that the axis of symmetry of the well matches the Z axis
        self.r_w = 1.0 # radius of the well
        self.h_w = 4.0 # water level near the well
        # Water level and corresponding radius at other location:
        self.r_i = 5.0
        self.h_i = 5.0
        
        self.ComputeGradientOfPorePressure()
    
    def ComputeWaterLevelSlope(self):
        
        self.slope = (self.h_i - self.h_w) / (self.radius * math.log(self.r_i / self.r_w))
        
    def ComputeGradientOfPorePressure(self):
        
        values = Array3()
        for element in self.spheres_model_part.Elements:
            
            coord_X = element.GetNode(0).X
            coord_Y = element.GetNode(0).Y
            self.radius = math.sqrt(coord_X * coord_X + coord_Y * coord_Y)
            
            self.ComputeWaterLevelSlope()

            if self.radius < self.r_w: # particles are inside the well
                values[0] = values[1] = 0.0
            else:
                values[0] = coord_X * self.slope / self.radius
                values[1] = coord_Y * self.slope / self.radius
            values[2] = 0.0
            
            element.GetNode(0).SetSolutionStepValue(PORE_PRESSURE_GRADIENT, values)
            

if __name__ == "__main__":

    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        
    project_parameters.AddEmptyValue("PostPorePressureGradientForce")
    project_parameters["PostPorePressureGradientForce"].SetBool(True)

    model = KratosMultiphysics.Model()
    DEMPorePressureAnalysisStage(model, project_parameters).Run()
