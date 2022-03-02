import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Importing other libraries
import math

class RVEAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):

        # input parameters of the analysis
        self.boundary_mp_name = project_parameters["rve_settings"]["boundary_mp_name"].GetString()
        self.averaging_mp_name = project_parameters["rve_settings"]["averaging_mp_name"].GetString()
        self.print_rve_post = project_parameters["rve_settings"]["print_rve_post"].GetBool()
        self.perturbation = project_parameters["rve_settings"]["perturbation"].GetDouble()

        self.averaging_volume = -1.0  # it will be computed in initialize
        domain_size = project_parameters["solver_settings"]["domain_size"].GetInt()
        if(domain_size == 2):
            self.strain_size = 3
        else:
            self.strain_size = 6

        # Pseudo time to be used for output
        self.time = 0.0

        self.populate_search_eps = 1e-4 ##tolerance in finding which conditions belong to the surface (will be multiplied by the lenght of the diagonal)
        self.geometrical_search_tolerance = 1e-4 #tolerance to be used in the search of the condition it falls into

        super().__init__(model, project_parameters)
        print("Constructor finalized")

    # Here populate the submodelparts to be used for periodicity
    def ModifyInitialGeometry(self):
        super().ModifyInitialGeometry()

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]

        # Construct auxiliary modelparts
        self.min_corner, self.max_corner = self._DetectBoundingBox(averaging_mp)
        self._ConstructFaceModelParts(self.min_corner, self.max_corner, boundary_mp)

        self.averaging_volume = (self.max_corner[0]-self.min_corner[0]) * (self.max_corner[1]-self.min_corner[1]) * (self.max_corner[2]-self.min_corner[2])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "RVE undeformed averaging volume = ", self.averaging_volume)
        print("Initial Geometry Modified")

    def InitializeSolutionStep(self):
        raise Exception("Should use the _CustomInitializeSolutionStep instead of this")

    def __CustomInitializeSolutionStep(self, strain, boundary_mp, averaging_mp):
        #reset position
        KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(averaging_mp.Nodes)

        self.ApplyBoundaryConditions()  # here the processes are called
        print("appied BCs")
        # construct MPCs according to the provided strain
        print("Custom init")
        self._ApplyPeriodicity(strain, averaging_mp, boundary_mp)
        print("Periodicity applied")
        # apply BCs for RVE according to the provided strain
        self._ApplyMinimalConstraints(
            averaging_mp, strain, self.min_corner, self.max_corner)

        self.ChangeMaterialProperties()  # this is normally empty
        self._GetSolver().InitializeSolutionStep()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def RunSolutionLoop(self):

        perturbation = self.perturbation

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]
        
        print("before strain ")
        stress_and_strain = []
        if(self.strain_size == 3):  # 2D case - ordering s00 s11 s01            
            stress_and_strain.append(self._ComputeEquivalentStress(0, 0, perturbation, boundary_mp, averaging_mp))
            print("first application of f.12, 2D")
            # stress_and_strain.append(self._ComputeEquivalentStress(1, 1, perturbation, boundary_mp, averaging_mp))
            # stress_and_strain.append(self._ComputeEquivalentStress(0, 1, perturbation, boundary_mp, averaging_mp))
        elif(self.strain_size == 6):  # 3D case - ordering:  s00 s11 s22 s01 s12 s02
            stress_and_strain.append(self._ComputeEquivalentStress(0, 0, perturbation, boundary_mp, averaging_mp)) 
            '''When ComputeEquivalentStress is called the i and j arguments, in the non-commented line [0, 0],
            are used to fill the strain matrix in the given indices'''
            print("First application of f.12, 3D")
            # stress_and_strain.append(self._ComputeEquivalentStress(1, 1, perturbation, boundary_mp, averaging_mp))
            # stress_and_strain.append(self._ComputeEquivalentStress(2, 2, perturbation, boundary_mp, averaging_mp))
            # stress_and_strain.append(self._ComputeEquivalentStress(0, 1, perturbation, boundary_mp, averaging_mp))
            # stress_and_strain.append(self._ComputeEquivalentStress(1, 2, perturbation, boundary_mp, averaging_mp))
            # stress_and_strain.append(self._ComputeEquivalentStress(0, 2, perturbation, boundary_mp, averaging_mp)) 

        C = self._ComputeEquivalentElasticTensor(stress_and_strain, perturbation)
        averaging_mp.SetValue(KratosMultiphysics.StructuralMechanicsApplication.ELASTICITY_TENSOR, C)
        self._MatrixOutput(C)
        print("Finished RunSolutionLoop")

    def _DetectBoundingBox(self, mp):
        print("start f.5")
        min_corner = KratosMultiphysics.Array3()
        min_corner[0] = 1e20
        min_corner[1] = 1e20
        min_corner[2] = 1e20

        max_corner = KratosMultiphysics.Array3()
        max_corner[0] = -1e20
        max_corner[1] = -1e20
        max_corner[2] = -1e20

        for node in mp.Nodes:
            x = node.X
            min_corner[0] = min(min_corner[0], x)
            max_corner[0] = max(max_corner[0], x)

            y = node.Y
            min_corner[1] = min(min_corner[1], y)
            max_corner[1] = max(max_corner[1], y)

            z = node.Z
            min_corner[2] = min(min_corner[2], z)
            max_corner[2] = max(max_corner[2], z)

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Boundng box detected")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Min. corner = ", min_corner)
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Max. corner = ", max_corner)

        print("End f.5")
        return min_corner, max_corner

    def __PopulateMp(self, face_name, coordinate, component, eps, mp):
        print("start f.6")
        if mp.NumberOfConditions() == 0:
            raise Exception("Boundary_mp is expected to have conditions and has none")

        mp = mp.GetRootModelPart()

        if not mp.HasSubModelPart(face_name):
            mp.CreateSubModelPart(face_name)
        face_mp = mp.GetSubModelPart(face_name)

        for cond in mp.Conditions:
            xc = cond.GetGeometry().Center()
            if abs(xc[component]-coordinate) < eps:
                face_mp.AddCondition(cond)

        node_ids = set()
        for cond in face_mp.Conditions:
            for node in cond.GetNodes():
                if(not node.Is(KratosMultiphysics.SLAVE)):
                    node_ids.add(node.Id)
                    node.Set(KratosMultiphysics.SLAVE)

        face_mp.AddNodes(list(node_ids))
        print("End f.6")
        return face_mp

    def _ConstructFaceModelParts(self, min_corner, max_corner, mp):
        
        print("Start f.7")
        diag_vect = max_corner - min_corner
        diag_lenght = math.sqrt(diag_vect[0]**2+diag_vect[1]**2+diag_vect[2]**2 )
        eps = self.populate_search_eps*diag_lenght

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLAVE, False, mp.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.MASTER, False, mp.Nodes)

        # Populate the slave faces
        print("Call to f.6")
        self.max_x_face = self.__PopulateMp("max_x_face", max_corner[0], 0, eps, mp)
        self.max_y_face = self.__PopulateMp("max_y_face", max_corner[1], 1, eps, mp)
        self.max_z_face = self.__PopulateMp("max_z_face", max_corner[2], 2, eps, mp)

        # First populate the master faces (min)
        print("More calls to f.6")
        self.min_x_face = self.__PopulateMp("min_x_face", min_corner[0], 0, eps, mp)
        self.min_y_face = self.__PopulateMp("min_y_face", min_corner[1], 1, eps, mp)
        self.min_z_face = self.__PopulateMp("min_z_face", min_corner[2], 2, eps, mp)

        if self.min_x_face.NumberOfConditions() == 0:
            raise Exception("min_x_face has 0 conditions")
        if self.min_y_face.NumberOfConditions() == 0:
            raise Exception("min_y_face has 0 conditions")
        if self.min_z_face.NumberOfConditions() == 0:
            raise Exception("min_z_face has 0 conditions")
        print("End of f.7")

    def _SelectClosestNode(self, mp, coords):
        print("Start of f.8")
        min_distance = 1e30
        selected_node = 0
        for node in mp.Nodes:
            dx = node.X0 - coords[0]
            dy = node.Y0 - coords[1]
            dz = node.Z0 - coords[2]
            d = dx**2 + dy**2 + dz**2

            if(d < min_distance):
                selected_node = node
                min_distance = d

        print("end of f.8")
        return selected_node

    # prescribed conditions to avoid rigid body motions
    def _ApplyMinimalConstraints(self, mp, strain, min_corner, max_corner):
        print("Start of f.9")
        aux = KratosMultiphysics.Array3()

        # point coinciding with the min_corner
        node = self._SelectClosestNode(mp, min_corner)
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        coords_min_corner = KratosMultiphysics.Array3(node)
        coords_min_corner[0] = node.X0
        coords_min_corner[1] = node.Y0
        coords_min_corner[2] = node.Z0

        disp_min_corner = strain*coords_min_corner
        node.SetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT, 0, disp_min_corner)
        print("End of f.9")

    def _ComputeEquivalentElasticTensor(self, stress_and_strain, perturbation):
        print("Start f.10")
        C = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)

        for j in range(len(stress_and_strain)):
            stress = stress_and_strain[j][0]
            for i in range(self.strain_size):
                C[i, j] = stress[i]

        inverse_perturbation = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)
        inverse_perturbation.fill(0.0)
        if(self.strain_size == 3):
            inverse_perturbation[0, 0] = 1.0/perturbation
            inverse_perturbation[1, 1] = 1.0/perturbation
            inverse_perturbation[2, 2] = 0.5/perturbation
        else:
            inverse_perturbation[0, 0] = 1.0/perturbation
            inverse_perturbation[1, 1] = 1.0/perturbation
            inverse_perturbation[2, 2] = 1.0/perturbation
            inverse_perturbation[3, 3] = 0.5/perturbation
            inverse_perturbation[4, 4] = 0.5/perturbation
            inverse_perturbation[5, 5] = 0.5/perturbation

        C = C*inverse_perturbation
        print("End of f.10")
        return C

    def _MatrixOutput(self, C, filename="rve_elasticity_tensor.txt"):
        print("Start of f.11")
        f = open(filename, 'w')

        if(self.strain_size == 3):  # 2D
            f.write(str(C[0, 0]) + " " + str(C[0, 1]) +
                    " " + str(C[0, 2]) + "\n")
            f.write(str(C[1, 0]) + " " + str(C[1, 1]) +
                    " " + str(C[1, 2]) + "\n")
            f.write(str(C[2, 0]) + " " + str(C[2, 1]) +
                    " " + str(C[2, 2]) + "\n")

        elif(self.strain_size == 6):
            for i in range(6):
                f.write(str(C[i, 0]) + " " + str(C[i, 1]) + " " + str(C[i, 2]) +
                        " " + str(C[i, 3]) + " " + str(C[i, 4]) + " " + str(C[i, 5]) + "\n")

        print("End of f.11")
        f.close()

    def _ComputeEquivalentStress(self, i, j, perturbation, boundary_mp, averaging_mp):
        print("Start of f.12, PROBLEMATIC")
        # Here use a pseudotime for output
        self.time = self.time + 1.0
        averaging_mp.GetRootModelPart().CloneTimeStep(self.time)

        print("create Strain matrix")
        strain = KratosMultiphysics.Matrix(3, 3)  #Strain tensor for the given step, in this test all of them are fixed with values = perturbation
        strain.fill(0.0)
        strain[1, 1] = perturbation     # Despite fixating strain in one of the main axis the post-process does not change.(Needed to re-compile Kratos)
        # strain[i, j] = perturbation
        # strain[j, i] = perturbation
        # print("foo")   #HOW TO PRINT FROM 'HERE'? GO TO script CALLING THIS ANALYSIS STAGE AND CHANGE WARNING TO INFO

        strain_vector = KratosMultiphysics.Vector(self.strain_size)  # Strain vector in Voigt notation
        if(self.strain_size == 3):
            print("2D strain selected, populate strain vector")
            strain_vector[0] = strain[0, 0]
            strain_vector[1] = strain[1, 1]
            strain_vector[2] = 2.0*strain[1, 2]
        elif(self.strain_size == 6):
            print("3D strain selected, populate strain vector")
            strain_vector[0] = strain[0, 0]
            strain_vector[1] = strain[1, 1]
            strain_vector[2] = strain[2, 2]
            strain_vector[3] = 2.0*strain[0, 1]
            strain_vector[4] = 2.0*strain[1, 2]
            strain_vector[5] = 2.0*strain[0, 2]
            # ORDER in Voigt notation???? Shouldn't 3 be 4, 4 be 5, and 5 be 3? (Kratos orders Voigt notation in this way, let's be consistent with it)
        print("pre-CustomInitializeSolutionStep (f.3)")
        self.__CustomInitializeSolutionStep(strain, boundary_mp, averaging_mp)
        print("Custom init done or solution step")
        self._GetSolver().Predict()

        self._GetSolver().SolveSolutionStep()
        process_info = averaging_mp.ProcessInfo
        avg_stress = KratosMultiphysics.Vector(self.strain_size)
        avg_stress.fill(0.0)
        measured_volume = 0.0

        for elem in averaging_mp.Elements:
            tmp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, process_info)
            ngauss = len(tmp)
            A = elem.GetGeometry().Area()
            measured_volume += A
            # TODO: this is only valid for gauss points with the same weight. should be generalized
            Agauss = A/ngauss
            for item in tmp:
                avg_stress = avg_stress + item*Agauss

        self._GetSolver().Clear()

        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Measured volume = ", measured_volume)

        avg_stress /= self.averaging_volume

        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Applied strain = ", strain_vector)
        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Average stress = ", avg_stress)

        if self.print_rve_post:
            self.OutputSolutionStep()

        # Reset position of nodes
        KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(averaging_mp.Nodes)
        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, zero, averaging_mp.Nodes)
        return avg_stress, strain_vector

    def _ApplyPeriodicity(self, strain, volume_mp, boundary_mp):
        # clear
        print("start of periodicity application")
        for constraint in volume_mp.GetRootModelPart().MasterSlaveConstraints:
            constraint.Set(KratosMultiphysics.TO_ERASE)
        volume_mp.GetRootModelPart().RemoveMasterSlaveConstraintsFromAllLevels(
            KratosMultiphysics.TO_ERASE)
        print("start of periodicity function")
        dx = self.max_corner[0] - self.min_corner[0]
        dy = self.max_corner[1] - self.min_corner[1]
        dz = self.max_corner[2] - self.min_corner[2]

        periodicity_utility = KratosMultiphysics.RVEPeriodicityUtility(self._GetSolver().GetComputingModelPart())
        print("before periodicity")
        # assign periodicity to faces
        search_tolerance = self.geometrical_search_tolerance # Is search tolerance the tolerance for finding closest node for applying periodicity?
        periodicity_utility.AssignPeriodicity(self.min_x_face, self.max_x_face, strain, KratosMultiphysics.Vector([dx, 0.0, 0.0]),search_tolerance)
        print("after first periodicity application")
        periodicity_utility.AssignPeriodicity(self.min_y_face, self.max_y_face, strain, KratosMultiphysics.Vector([0.0, dy, 0.0]),search_tolerance)
        periodicity_utility.AssignPeriodicity(self.min_z_face, self.max_z_face, strain, KratosMultiphysics.Vector([0.0, 0.0, dz]),search_tolerance)
        # This Finalize function does not work for the 2D element
        print("before finalize step")
        periodicity_utility.Finalize(KratosMultiphysics.DISPLACEMENT)
        print("after finalize step")
        # *****************This point is not reached************** 
        '''So, it's the Finalize function that crashes for a 2D element
           After commenting line 371 the test finishes but the GiD post-process shows that it is obviously wrong 
           Strain is rotated along Z axis'''
        
        # start from the exact solution in the case of a constant strain
        x = KratosMultiphysics.Array3()
        for node in volume_mp.Nodes:
            x[0] = node.X0
            x[1] = node.Y0
            x[2] = node.Z0
            d = strain*x
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, d)
        print("End of periodicity application")
        
