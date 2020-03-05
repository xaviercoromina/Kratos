import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time
import math


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])


class DefineEmbeddedWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)
        self.model = Model

        self.main_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.wake_model_part=Model.CreateModelPart("wake")

        self.epsilon = settings["epsilon"].GetDouble()

    def ExecuteInitialize(self):
        ini_time = time.time()

        for node in self.main_model_part.Nodes:
            node.SetValue(CPFApp.AIRFOIL, False)


        self.list_of_failed_te_nodes =[]
        max_x_coordinate = -1e30
        restarted_search_maximum = 1e30
        trailing_edge_candidate = self.FindNode(max_x_coordinate, restarted_search_maximum)
        is_valid = self.CheckIfValid(trailing_edge_candidate)
        while(not is_valid):
            trailing_edge_candidate.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)
            self.list_of_failed_te_nodes.append(trailing_edge_candidate.Id)
            trailing_edge_candidate = self.FindNode(max_x_coordinate, trailing_edge_candidate.X)
            is_valid = self.CheckIfValid(trailing_edge_candidate)
            if not is_valid:
                self.list_of_failed_te_nodes.append(trailing_edge_candidate.Id)
                trailing_edge_candidate.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)

        avg_elem_num = 10
        avg_node_num = 10
        KratosMultiphysics.FindNodalNeighboursProcess(
            self.main_model_part, avg_elem_num, avg_node_num).Execute()
        self._DefineWakeModelPart()
        self._MoveAndRotateWake()

        # # Executing define wake process
        cpp_time = time.time()

        CPFApp.DefineEmbeddedWakeProcess(self.main_model_part, self.wake_model_part).Execute()
        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','CPP time: ',time.time()-cpp_time)

        list_of_structure_te_nodes= []

        self._RedefineWake()

#        self.main_model_part.GetNode(8003).SetValue(CPFApp.WING_TIP, True)

        # for elem in self.main_model_part.Elements:
        #     if elem.GetValue(CPFApp.WAKE):
        #         for node in elem.GetNodes():
        #             if node.GetValue(CPFApp.WING_TIP) and elem.IsNot(KratosMultiphysics.STRUCTURE):
        #                 elem.Set(KratosMultiphysics.STRUCTURE)
        #                 print("Setting element to structure", elem.Id)

        # self.main_model_part.GetElement(62304).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(62304).SetValue(CPFApp.WAKE, False)

        # self.main_model_part.GetElement(147531).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(75726).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(81767).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(81767).Set(KratosMultiphysics.STRUCTURE, True)
        # self.main_model_part.GetNode(74353).SetValue(CPFApp.TRAILING_EDGE, False)
        # self.main_model_part.GetNode(74353).SetValue(CPFApp.WING_TIP, False)
        # self.main_model_part.GetNode(43236).SetValue(CPFApp.WING_TIP, True)


        # self.main_model_part.GetElement(47403).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(36427).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(16270).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(49229).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(49229).Set(KratosMultiphysics.STRUCTURE, True)
        # self.main_model_part.GetNode(25841).SetValue(CPFApp.TRAILING_EDGE, False)
        # self.main_model_part.GetNode(25841).SetValue(CPFApp.WING_TIP, False)
        # self.main_model_part.GetNode(32563).SetValue(CPFApp.WING_TIP, True)

        # self.main_model_part.GetElement(21112).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(117327).Set(KratosMultiphysics.STRUCTURE, True)

        # self.main_model_part.GetElement(11562).Set(KratosMultiphysics.STRUCTURE, True)
        # self.main_model_part.GetElement(11562).SetValue(CPFApp.WAKE, False)
        # self.main_model_part.GetElement(11562).SetValue(CPFApp.KUTTA, True)
        # self.main_model_part.GetNode(43826).SetValue(CPFApp.WING_TIP, False)
        # self.main_model_part.GetNode(59302).SetValue(CPFApp.WING_TIP, False)

        # self.main_model_part.GetElement(80091).SetValue(CPFApp.WAKE, True)

        # self.main_model_part.GetElement(24090).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(23557).SetValue(CPFApp.WAKE, True)
        # self.main_model_part.GetElement(23557).Set(KratosMultiphysics.STRUCTURE, True)
        # # self.main_model_part.GetElement(36427).SetValue(CPFApp.KUTTA, True)
        list_of_strucutre_element_nodes_id = []
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.STRUCTURE):
                for node in elem.GetNodes():
                    list_of_strucutre_element_nodes_id.append(node.Id)
            if elem.IsNot(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if node.X>0.47:
                    # if node.X>0.495 and node.GetValue(CPFApp.LOWER_SURFACE):
                        node.SetValue(CPFApp.TRAILING_EDGE, True)
        for elem in self.main_model_part.Elements:
            lower_counter = 0
            is_trailing_edge = False

            for node in elem.GetNodes():
                if node.GetValue(CPFApp.LOWER_SURFACE):
                    lower_counter += 1
                if node.GetValue(CPFApp.TRAILING_EDGE):
                    is_trailing_edge = True
            if lower_counter == 3 and is_trailing_edge:
                elem.SetValue(CPFApp.KUTTA, True)

        for elem in self.main_model_part.Elements:
            geometry_elemental_distances = KratosMultiphysics.Vector(3,0.0)
            is_wing_tip = False
            for i_node, node in enumerate(elem.GetNodes()):
                if node.GetValue(CPFApp.TRAILING_EDGE) and elem.GetValue(CPFApp.KUTTA):
                    geometry_elemental_distances[i_node]=node.GetValue(KratosMultiphysics.TEMPERATURE)
                else:
                    geometry_elemental_distances[i_node]=node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                #if node.GetValue(CPFApp.WING_TIP):
                if node.Id in list_of_strucutre_element_nodes_id:
                    is_wing_tip = True
            elem.SetValue(CPFApp.GEOMETRY_ELEMENTAL_DISTANCES, geometry_elemental_distances)
            #if elem.GetValue(CPFApp.KUTTA):
            #    elem.Set(KratosMultiphysics.STRUCTURE,True)
            #    elem.SetValue(CPFApp.WAKE, True)
            wake_elemental_distances = elem.GetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES)
            wake_pos = 0
            wake_neg = 0
            geom_pos = 0
            geom_neg = 0
            for geom, wake in zip(geometry_elemental_distances, wake_elemental_distances):
                if geom>0:
                    geom_pos += 1
                else:
                    geom_neg += 1
                if wake>0:
                    wake_pos += 1
                else:
                    wake_neg += 1

            is_embedded = False
            if (geom_pos>0 and geom_neg > 0):
                is_embedded = True
            is_wake = False
            if (wake_pos>0 and wake_neg > 0):
                is_wake = True
            if not is_embedded:
                wake = elem.GetValue(CPFApp.WAKE)
              #  if (is_wake and not wake):
              #      if elem.GetValue(CPFApp.KUTTA):
                       # print("Elem non embedded #", elem.Id, -1*wake_elemental_distances)
                   #     elem.SetValue(CPFApp.GEOMETRY_ELEMENTAL_DISTANCES, -1*wake_elemental_distances)
               #     else:
                       # print("Elem non embedded #", elem.Id, -1*wake_elemental_distances)
                       # elem.SetValue(CPFApp.GEOMETRY_ELEMENTAL_DISTANCES, wake_elemental_distances)
            #if is_embedded and is_wake and is_wing_tip:
            if is_embedded and is_wake:
                if elem.GetValue(CPFApp.KUTTA):
                    print("Elem #", elem.Id, -1*wake_elemental_distances)
                    # elem.SetValue(CPFApp.GEOMETRY_ELEMENTAL_DISTANCES, -1*wake_elemental_distances)
                else:
                    pass



        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)
    def _RedefineWake(self):
        ini_time = time.time()

        initial_point = self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN)
        max_inactive_x = initial_point[0]
        lower_surface_elem_ids=[]
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.TO_SPLIT) and elem.Is(KratosMultiphysics.ACTIVE):

                wake_origin = self.moving_parameters["origin"].GetVector()
                wake_angle = self.moving_parameters["rotation_angle"].GetDouble()
                wake_normal = KratosMultiphysics.Vector(3, 0.0)
                wake_normal[0] = math.sin(-wake_angle)
                wake_normal[1] = math.cos(-wake_angle)
                skin_normal = -1*elem.GetValue(CPFApp.VELOCITY_LOWER)
                projection = wake_normal[0]*(skin_normal[0]) + wake_normal[1]*skin_normal[1]

                if projection > 0.0:
                    for node in elem.GetNodes():
                        node.SetValue(CPFApp.UPPER_SURFACE, True)

                elif projection < 0.0:
                    for node in elem.GetNodes():
                        node.SetValue(CPFApp.LOWER_SURFACE, True)

                # for node in elem.GetNodes():
                #     distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                #     if distance > 0.0:
                #         wake_origin = self.moving_parameters["origin"].GetVector()
                #         wake_angle = self.moving_parameters["rotation_angle"].GetDouble()
                #         wake_normal = KratosMultiphysics.Vector(3, 0.0)
                #         wake_normal[0] = math.sin(-wake_angle)
                #         wake_normal[1] = math.cos(-wake_angle)
                #         projection = wake_normal[0]*(node.X-wake_origin[0]) + wake_normal[1]*(node.Y-wake_origin[1])

                #         if projection > 0.0:
                #             elem.SetValue(CPFApp.KUTTA, False)
                #             for node in elem.GetNodes():
                #                 node.SetValue(CPFApp.UPPER_SURFACE, True)
                #         elif projection < 0.0:
                #             for node in elem.GetNodes():
                #                 node.SetValue(CPFApp.LOWER_SURFACE, True)
        for elem in self.main_model_part.Elements:
            counter = 0
            for node in elem.GetNodes():
                if node.GetValue(CPFApp.LOWER_SURFACE):
                    counter += 1
            if counter == 3:
                lower_surface_elem_ids.append(elem.Id)



                # for node in elem.GetNodes():
                #     distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                #     if distance > 0.0 and not node.GetValue(CPFApp.TRAILING_EDGE):
                #         if not node.GetValue(CPFApp.WAKE_DISTANCE) == 0.0:
                #             defining_node = node
                #             not_found = False
                #         else:
                #             null_wake_distance = True
                #     if null_wake_distance:
                #         for node in elem.GetNodes():
                #             distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                #             if distance < 0.0 and not node.GetValue(CPFApp.TRAILING_EDGE):
                #                 if not node.GetValue(CPFApp.WAKE_DISTANCE) == 0.0:
                #                     defining_node = node
                #                     not_found = False
                # if not_found:
                #     defining_node = node

                # if defining_node.GetValue(CPFApp.WAKE_DISTANCE) > 0.0:
                #     elem.SetValue(CPFApp.KUTTA, False)
                #     for node in elem.GetNodes():
                #         node.SetValue(CPFApp.UPPER_SURFACE, True)
                # elif defining_node.GetValue(CPFApp.WAKE_DISTANCE) < 0.0:
                #     lower_surface_elem_ids.append(elem.Id)
                #     for node in elem.GetNodes():
                #         node.SetValue(CPFApp.LOWER_SURFACE, True)

                # if node_pos.Y - self.trailing_edge_node.Y>= 0.0:
                #     elem.SetValue(CPFApp.KUTTA, False)
                #     for node in elem.GetNodes():
                #         node.SetValue(CPFApp.UPPER_SURFACE, True)
                # else:
                #     lower_surface_elem_ids.append(elem.Id)
                #     for node in elem.GetNodes():
                #         node.SetValue(CPFApp.LOWER_SURFACE, True)

        if self.main_model_part.HasSubModelPart("lower_surface_sub_model_part"):
            self.main_model_part.RemoveSubModelPart("lower_surface_sub_model_part")
        self.lower_surface_element_sub_model_part = self.main_model_part.CreateSubModelPart("lower_surface_sub_model_part")
        self.lower_surface_element_sub_model_part.AddElements(lower_surface_elem_ids)

        for node in self.main_model_part.Nodes:
            if node.GetValue(CPFApp.UPPER_SURFACE) and node.GetValue(CPFApp.LOWER_SURFACE) and node.X > max_inactive_x-0.1 and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
            # if node.GetValue(CPFApp.UPPER_SURFACE) and node.GetValue(CPFApp.LOWER_SURFACE) and node.X < self.trailing_edge_node.X  and node.X > max_inactive_x-0.1 and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
                node.SetValue(CPFApp.TRAILING_EDGE, True)
                print("SETTING NODE TO TRAILING EDGE", node.Id)

        for elem in self.lower_surface_element_sub_model_part.Elements:
            if not elem.GetValue(CPFApp.WAKE) and not elem.GetValue(CPFApp.KUTTA) and elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if node.GetValue(CPFApp.TRAILING_EDGE) and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
                        elem.SetValue(CPFApp.KUTTA, True)
                    if node.GetValue(CPFApp.WING_TIP) and node.GetValue(CPFApp.WAKE_DISTANCE) > 0.0:
                        elem.SetValue(CPFApp.KUTTA, True)
                        node.SetValue(CPFApp.TRAILING_EDGE, True)

        for node in self.main_model_part.Nodes:
            geometry_distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
            node.SetValue(KratosMultiphysics.TEMPERATURE, geometry_distance)
            if node.GetValue(CPFApp.TRAILING_EDGE) and geometry_distance< 0.0:
                node.Set(KratosMultiphysics.INSIDE) # is negative?
                node.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, -1e30)
                node.SetValue(KratosMultiphysics.TEMPERATURE, -1e30)
        for elem in self.main_model_part.Elements:
            is_trailing_edge = False
            for node in elem.GetNodes():
                if node.GetValue(CPFApp.TRAILING_EDGE):
                    is_trailing_edge = True
            if is_trailing_edge:
                elem_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                if elem.GetValue(CPFApp.KUTTA):
                    for node, elem_distance in zip(elem.GetNodes(), elem_distances):
                        if node.GetValue(CPFApp.TRAILING_EDGE):
                            old_distance = node.GetValue(KratosMultiphysics.TEMPERATURE)
                            if old_distance < 0.0:
                                new_distance = max(-abs(old_distance),-abs(elem_distance))
                            else:
                                new_distance = min(abs(old_distance),abs(elem_distance))
                            node.SetValue(KratosMultiphysics.TEMPERATURE, new_distance)
                            print(elem.Id, node.Id, elem_distance, node.Is(KratosMultiphysics.INSIDE), new_distance)
                else:
                    for node, elem_distance in zip(elem.GetNodes(), elem_distances):
                        if node.GetValue(CPFApp.TRAILING_EDGE):
                            old_distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
                            if old_distance < 0.0:
                                new_distance = max(-abs(old_distance),-abs(elem_distance))
                            else:
                                new_distance = min(abs(old_distance),abs(elem_distance))
                            node.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, new_distance)
                            print(elem.Id, node.Id, elem_distance, node.Is(KratosMultiphysics.INSIDE), new_distance)
        print("List of failed te nodes", self.list_of_failed_te_nodes)
        for node_id in self.list_of_failed_te_nodes:
            self.main_model_part.GetNode(node_id).SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)

        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Redefine time: ',time.time()-ini_time)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''
        # self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        # self.wake_model_part.CreateNewNode(2, 200.0, 0.0, 0.0)
        # self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))

    def _MoveAndRotateWake(self):
        ''' This function moves and rotates the wake with the same parameters as the geometry.
        '''
        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("origin")
        self.moving_parameters["origin"].SetVector(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        self.moving_parameters.AddEmptyValue("rotation_angle")
        angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        self.moving_parameters["rotation_angle"].SetDouble(angle)
        # print(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        # print(angle)

        # inital_point = self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN)
        # vector = KratosMultiphysics.Vector(3, 0.0)
        # vector[0] =  math.cos(angle)
        # vector[1] = -math.sin(-angle)

        # max_projection = -1e30
        # for element in self.main_model_part.Elements:
        #     if element.IsNot(KratosMultiphysics.ACTIVE):
        #         x_center = element.GetGeometry().Center().X
        #         y_center = element.GetGeometry().Center().X
        #         center_vector = KratosMultiphysics.Vector(3, 0.0)
        #         center_vector[0] = x_center - inital_point[0]
        #         center_vector[1] = y_center - inital_point[1]
        #         product = center_vector[0]*vector[0] + center_vector[1]*vector[1]
        #         if product > max_projection:
        #             max_projection = product
        #             max_y = element.GetGeometry().Center().Y
        #             max_x = element.GetGeometry().Center().X
        #             print(element.Id, max_projection, max_x, max_y)
        # wake_origin = KratosMultiphysics.Vector(3, 0.0)
        # wake_origin[0] = max_x - 0.00
        # wake_origin[1] = max_y


        max_x = -1e30
        for element in self.main_model_part.Elements:
            if element.IsNot(KratosMultiphysics.ACTIVE):
                x_center = element.GetGeometry().Center().X
                if x_center > max_x:
                    max_x = x_center
                    max_y = element.GetGeometry().Center().Y
        wake_origin = KratosMultiphysics.Vector(3, 0.0)

        # max_x_negative = -1e30
        # for node in self.main_model_part.Nodes:
        #     distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
        #     if distance < 0.0 and node.X > max_x_negative:
        #         max_x_negative = node.X
        #         max_y_negative = node.Y + 1e-7
        # wake_origin = KratosMultiphysics.Vector(3, 0.0)

        wake_origin[0] = max_x
        # wake_origin[1] = max_y
        self.wake_model_part.CreateNewNode(1, max_x, max_y, 0.0)
        self.wake_model_part.CreateNewNode(2, self.model.GetModelPart("skin").GetNode(1).X, self.model.GetModelPart("skin").GetNode(1).Y, 0.0)
        self.wake_model_part.CreateNewNode(3, 200.0, self.model.GetModelPart("skin").GetNode(1).Y, 0.0)
        # self.wake_model_part.CreateNewNode(1, max_x, max_y, 0.0)
        # self.wake_model_part.CreateNewNode(2, max_x_negative, max_y_negative, 0.0)
        # self.wake_model_part.CreateNewNode(3, 200.0, max_y_negative, 0.0)
        self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element2D2N", 2, [2,3], KratosMultiphysics.Properties(0))


        # wake_origin[0] = self.model.GetModelPart("skin").GetNode(5000).X - 0.0010
        # wake_origin[1] = self.model.GetModelPart("skin").GetNode(5000).Y
        print(self.model.GetModelPart("skin").GetNode(1).X, self.model.GetModelPart("skin").GetNode(1).Y)
        # self.moving_parameters["origin"].SetVector(wake_origin)
        # self.moving_parameters.AddEmptyValue("rotation_angle")
        # self.moving_parameters["rotation_angle"].SetDouble(0.0)


        # CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

    def ExecuteFinalizeSolutionStep(self):
        self.wake_sub_model_part = self.main_model_part.CreateSubModelPart("wake_sub_model_part")
        for elem in self.main_model_part.Elements:
            if elem.GetValue(CPFApp.WAKE) and elem.Is(KratosMultiphysics.ACTIVE):
                self.wake_sub_model_part.Elements.append(elem)

        absolute_tolerance = 1e-6
        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled2D(self.wake_sub_model_part, absolute_tolerance, 2)
        self.main_model_part.RemoveSubModelPart("wake_sub_model_part")

    def CheckIfValid(self, trailing_edge_candidate):
        for elem in self.main_model_part.Elements:
            is_neighbour = False
            for node in elem.GetNodes():
                if (node.Id == trailing_edge_candidate.Id):
                    is_neighbour = True
            if is_neighbour:
                if elem.IsNot(KratosMultiphysics.ACTIVE):
                    return True
                for node in elem.GetNodes():
                    if node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0 and node.X < trailing_edge_candidate.X :
                        is_valid = self.CheckIfValid(node)
                        if is_valid:
                            return True
        return False

    def FindNode(self, max_x_coordinate, restarted_search_maximum):
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.TO_SPLIT) and elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if(node.X > max_x_coordinate) and (node.X<restarted_search_maximum) and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
                        max_x_coordinate = node.X
                        trailing_edge_node = node
        return trailing_edge_node
