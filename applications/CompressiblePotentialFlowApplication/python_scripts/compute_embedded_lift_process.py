import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_lift_process import ComputeLiftProcess

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeEmbeddedLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeEmbeddedLiftProcess(ComputeLiftProcess):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the surface nodes",
            "moment_reference_point" : [0.0,0.0,0.0],
            "far_field_model_part_name": "",
            "is_infinite_wing": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)
        self.fluid_model_part=Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.reference_area =  self.fluid_model_part.ProcessInfo.GetValue(CPFApp.REFERENCE_CHORD)
        self.resultant_force=KratosMultiphysics.Vector(3)


        far_field_model_part_name = settings["far_field_model_part_name"].GetString()
        self.compute_far_field_forces = far_field_model_part_name != ""
        if self.compute_far_field_forces :
            self.far_field_model_part = Model[far_field_model_part_name]

        self.reference_area =  self.fluid_model_part.ProcessInfo.GetValue(CPFApp.REFERENCE_CHORD)
        self.moment_reference_point = settings["moment_reference_point"].GetVector()
        self.is_infinite_wing = settings["is_infinite_wing"].GetBool()

        if not self.reference_area > 0.0:
            raise Exception('The reference area should be larger than 0.')

    def _ComputeLiftFromPressure(self):
        if (self.fluid_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)==2):
            CPFApp.ComputeEmbeddedLiftProcess2D(self.fluid_model_part,self.resultant_force).Execute()
        elif (self.fluid_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)==3):
            CPFApp.ComputeEmbeddedLiftProcess3D(self.fluid_model_part,self.resultant_force).Execute()
        else:
            raise(Exception("Dimension of the problem is different than 2 or 3."))

        max_x = -1e10
        min_x = 1e10
        import math
        wake_direction = KratosMultiphysics.Vector(3, 0.0)
        wake_direction[0] = math.cos(math.radians(5.0))
        wake_direction[1] = math.sin(math.radians(5.0))
        for elem in self.fluid_model_part.Elements:
            npos = 0
            nneg = 0
            is_embedded = False
            for node in elem.GetNodes():
                distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                if distance>0.0:
                    npos += 1
                else:
                    nneg += 1

            if npos > 0 and nneg > 0:
                is_embedded = True
            if is_embedded:
                center = elem.GetGeometry().Center()
                projection = DotProduct(center, wake_direction)
                if projection > max_x:
                    max_x = projection
                    max_center = center
                if projection < min_x:
                    min_x = projection
                    min_center = center

        distance_vector = max_center - min_center
        print("CHORD COMPUTED", distance_vector.norm_2())
        # // Mark nodes of the furthest deactivated element and store its neighbour elements
        # std::size_t max_node_i = -1;
        # double max_x = -1e10;
        # auto& r_max_element_geometry=p_max_elem->GetGeometry();

        # for (unsigned int i_node= 0; i_node < r_max_element_geometry.size(); i_node++) {
        # // r_max_element_geometry[i_node].SetValue(AIRFOIL,true);
        # // BoundedVector<double, 3> nodal_position_vector = r_max_element_geometry[i_node].Coordinates() - wake_origin;
        # BoundedVector<double, 3> nodal_position_vector = r_max_element_geometry[i_node].Coordinates();
        # double projection = std::abs(inner_prod(nodal_position_vector, wake_direction));
        # // if (r_max_element_geometry[i_node].X()>max_x) {
        # if (projection>max_x) {
        # // if (projection>max_x) {
        #     max_x = projection;
        #     max_node_i = i_node;
        # }

        self.reference_area = distance_vector.norm_2()
        self._CalculateWakeTangentAndNormalDirections()

        self._ProjectForceToFreeStreamVelocity(self.resultant_force/self.reference_area)

        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cl = ', self.lift_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cd = ', self.drag_coefficient)

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT, self.lift_coefficient)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT, self.drag_coefficient)

    def _ComputeMomentFromPressure(self):
        # TODO Add implementation for embedded bodies
        pass

    def _GetTrailingEdgeNode(self):
        for node in self.fluid_model_part.GetSubModelPart("trailing_edge_sub_model_part").Nodes:
            self.te = node
            break

