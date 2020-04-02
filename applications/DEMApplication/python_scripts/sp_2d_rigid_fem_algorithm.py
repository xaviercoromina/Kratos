import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import mesh_creator_sphere_2D
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import math

class DEMAnalysisStage2DSpRigidFem(DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super(DEMAnalysisStage2DSpRigidFem, self).__init__(model, project_parameters)
        # TEST NUMBER:
        # 1. CTW16
        # 2: CTW10
        # 3. CTW13
        # 4: BLIND
        # self.test_number = 4
        self.test_number = project_parameters["test_number"].GetInt()

    def Initialize(self):
        super(DEMAnalysisStage2DSpRigidFem, self).Initialize()

        self.CreateSPMeasuringRingSubmodelpart(self.spheres_model_part)
        porosity = self.ComputePorosity(self.spheres_model_part)
        print(porosity)
        dises

    def ComputeSpecimenFullVolume(self):

        outer_radius = 0.0508 # This is 2 inches
        blind_side = 0.3048 # This is 12 inches
        unit_depth = 1.0
        inner_radius = 0.0
        if self.test_number == 1:
            inner_radius = 0.00381 # CTW16 specimen
        elif self.test_number == 2:
            inner_radius = 0.0127 # CTW10 specimen
        elif self.test_number == 3:
            inner_radius = 0.00635 # CTW13 specimen
        else: #self.test_number == 4:
            inner_radius = 0.0381 # Blind Test

        if self.test_number < 4:
            return math.pi * (outer_radius * outer_radius - inner_radius * inner_radius) * unit_depth
        else:
            return (blind_side * blind_side - math.pi * inner_radius * inner_radius) * unit_depth

    def ComputePorosity(self, spheres_model_part):

        total_circles_volume = 0.0 # Remember we are here in 2D
        unit_depth = 1.0 # Remember we are here in 2D
        for node in spheres_model_part.Nodes:
            radius = node.GetSolutionStepValue(Kratos.RADIUS)
            total_circles_volume +=  math.pi * radius * radius * unit_depth

        total_spheres_volume = math.pi * total_circles_volume / 6.0

        specimen_porosity = 1.0 - total_spheres_volume / self.ComputeSpecimenFullVolume()
        print(specimen_porosity)
        sescalla
        # Porosity values in sandstones are in the rage of 25%
        target_porosity = 0.25
        SP_multiplier = (1.0 - target_porosity) / (1.0 - specimen_porosity)

        return SP_multiplier

    def CreateSPMeasuringRingSubmodelpart(self, spheres_model_part):

        if not self.spheres_model_part.HasSubModelPart("RingSubmodelPart"):
            self.spheres_model_part.CreateSubModelPart('RingSubmodelPart')
        self.ring_submodelpart = self.spheres_model_part.GetSubModelPart('RingSubmodelPart')

        if self.test_number == 1:
            self.zone_radius_to_measure_2d_sp = 0.02 # In meters. Maybe too much
        elif self.test_number == 2:
            self.zone_radius_to_measure_2d_sp = 0.02 # Maybe too much before, it was 0.03...
        elif self.test_number == 3:
            self.zone_radius_to_measure_2d_sp = 0.02 # In meters.
        else: #self.test_number == 4:
            self.zone_radius_to_measure_2d_sp = 0.06

        nodes_in_zone_radius_list = []
        elements_in_zone_radius_list = []

        for element in self.spheres_model_part.Elements:
            node = element.GetNode(0)
            x = node.X
            y = node.Y

            if (x * x + y * y) < self.zone_radius_to_measure_2d_sp * self.zone_radius_to_measure_2d_sp:
                nodes_in_zone_radius_list.append(node.Id)
                elements_in_zone_radius_list.append(element.Id)

        self.ring_submodelpart.AddNodes(nodes_in_zone_radius_list)
        self.ring_submodelpart.AddElements(elements_in_zone_radius_list)

    def PrintResultsForGid(self, time):
        super(DEMAnalysisStage2DSpRigidFem, self).PrintResultsForGid(time)

        DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.ring_submodelpart)

        self.SettingGeometricalSPValues()

        self.creator_destructor.MarkParticlesForErasingGivenCylinder(self.ring_submodelpart, self.center, self.axis, self.radius)

        DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.rigid_face_model_part, self.time)
        DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.rigid_face_model_part, self.time)

    def SettingGeometricalSPValues(self):

        self.center = KratosMultiphysics.Array3()
        self.center[0] = 0; # self.sp_parameters["problem_data"]["center"][0].GetDouble()
        self.center[1] = 0; # self.sp_parameters["problem_data"]["center"][1].GetDouble()
        self.center[2] = 0; # self.sp_parameters["problem_data"]["center"][2].GetDouble()
        self.axis = KratosMultiphysics.Array3()
        self.axis[0] = 0; # self.sp_parameters["problem_data"]["axis"][0].GetDouble()
        self.axis[1] = 0; # self.sp_parameters["problem_data"]["axis"][1].GetDouble()
        self.axis[2] = 1; # self.sp_parameters["problem_data"]["axis"][2].GetDouble()

        self.radius = 0.0
        if self.test_number == 1:
            self.radius = 0.0015; #0.0036195; #95% of the real hole. CTW16 specimen
        elif self.test_number == 2:
            self.radius = 0.005; #0.012065;  #95% of the real hole. CTW10 specimen
        elif self.test_number == 3:
            self.radius = 0.0060325; #95% of the real hole. CTW13 specimen
        else: #self.test_number == 4:
            self.radius = 0.015; #0.036195;  #95% of the real hole. Blind Test

    def AdditionalFinalizeOperations(self):

        spheres_mp_filename_post = self.problem_name + 'DEM_Post'

        if self.write_mdpa_from_results:
            mesh_creator_sphere_2D.WriteSphereMdpaFromResults(self.problem_name + 'DEM', self.main_path, spheres_mp_filename_post, self.file_msh, self.file_res, self.post_path)

if __name__ == "__main__":
    DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
