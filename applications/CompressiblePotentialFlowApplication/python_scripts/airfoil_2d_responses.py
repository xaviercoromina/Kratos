import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import math

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

class AngleOfAttackResponseFunction(ResponseFunctionInterface):
    def __init__(self, response_id, response_settings, model):
        self.model = model
        if not response_settings.Has("trailing_edge_model_part")  or not response_settings.Has("leading_edge_model_part"):
            raise(Exception("Please define a model part with the trailing edge node and a model part with the leading edge node"))

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_sub_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]
        # self.free_stream_velocity = self.main_model_part.ProcessInfo[KCPFApp.FREE_STREAM_VELOCITY]
        self.free_stream_velocity = KratosMultiphysics.Vector(3, 0.0)
        self.free_stream_velocity[0] = 1.0
        for node in self.model[self.trailing_edge_model_part_name].Nodes:
            self.te_node = node
            break
        for node in self.model[self.leading_edge_sub_model_part_name].Nodes:
            self.le_node = node
            break

    def CalculateValue(self):
        self.aoa = self._CalculateAOA(self.te_node.X,self.te_node.Y, self.le_node.X, self.le_node.Y)

    def CalculateGradient(self):
        epsilon = 1e-9
        self.te_x_gradient = (self._CalculateAOA(self.te_node.X + epsilon,self.te_node.Y, self.le_node.X, self.le_node.Y) - self.aoa) / epsilon
        self.te_y_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y+ epsilon, self.le_node.X, self.le_node.Y) - self.aoa) / epsilon
        self.le_x_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y, self.le_node.X+ epsilon, self.le_node.Y) - self.aoa) / epsilon
        self.le_y_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y, self.le_node.X, self.le_node.Y+ epsilon) - self.aoa) / epsilon

    def GetValue(self):
        return self.aoa

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient

    def _CalculateAOA(self, te_x, te_y, le_x, le_y):

        velocity_norm = self.free_stream_velocity.norm_2()
        chord_vector = KratosMultiphysics.Vector(3, 0.0)
        chord_vector[0] = le_x-te_x
        chord_vector[1] = le_y-te_y
        chord_norm = chord_vector.norm_2()
        aoa = math.acos(DotProduct(self.free_stream_velocity, -1*chord_vector)/(velocity_norm*chord_norm))

        return aoa

class ChordLengthResponseFunction(ResponseFunctionInterface):


    def __init__(self, response_id, response_settings, model):
        self.model = model

        if not response_settings.Has("trailing_edge_model_part")  or not response_settings.Has("leading_edge_model_part"):
            raise(Exception("Please define a model part with the trailing edge node and a model part with the leading edge node"))

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_sub_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]

        for node in self.model[self.trailing_edge_model_part_name].Nodes:
            self.te_node = node
            break
        for node in self.model[self.leading_edge_sub_model_part_name].Nodes:
            self.le_node = node
            break

    def _ComputeChord(self, te_x, te_y, le_x, le_y):
        chord = math.sqrt((te_x-le_x)**2+(te_y-le_y)**2)

        return chord

    def CalculateValue(self):
        self.chord = self._ComputeChord(self.te_node.X,self.te_node.Y, self.le_node.X,self.le_node.Y)

    def CalculateGradient(self):
        epsilon = 1e-9
        self.te_x_gradient = (self._ComputeChord(self.te_node.X + epsilon,self.te_node.Y, self.le_node.X, self.le_node.Y) - self.chord) / epsilon
        self.te_y_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y+ epsilon, self.le_node.X, self.le_node.Y) - self.chord) / epsilon
        self.le_x_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y, self.le_node.X+ epsilon, self.le_node.Y) - self.chord) / epsilon
        self.le_y_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y, self.le_node.X, self.le_node.Y+ epsilon) - self.chord) / epsilon

    def GetValue(self):
        return self.chord

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient



class PerimeterResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]
        self.body_model_part = self.main_model_part.GetSubModelPart("Body2D_Body")

    def _ComputePerimeter(self,  model_part):
        return KCPFApp.PotentialFlowUtilities.CalculateArea(model_part.Conditions)

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):

        return self._ComputePerimeter(self.body_model_part)

    def GetNodalGradient(self, variable):
        gradient = {}
        epsilon = 1e-6
        initial_perimeter =  self._ComputePerimeter(self.body_model_part)

        for node in self.body_model_part.Nodes:
            shape_gradient = KratosMultiphysics.Vector(3, 0.0)


            node.X += epsilon
            current_perimeter =  self._ComputePerimeter(self.body_model_part)
            shape_gradient[0] = (current_perimeter-initial_perimeter)/epsilon
            node.X -= epsilon

            node.Y += epsilon
            current_perimeter =  self._ComputePerimeter(self.body_model_part)
            shape_gradient[1] = (current_perimeter-initial_perimeter)/epsilon
            node.Y -= epsilon

            gradient[node.Id] = shape_gradient

        return gradient

class FixPointResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model
        self.traced_node_id = response_settings["traced_node"].GetInt()

    def Initialize(self):
        self.body_model_part = self.model["MainModelPart.Body2D_Body"]
        self.traced_node = self.model["MainModelPart"].GetNode(self.traced_node_id)
        self.initial_x = self.traced_node.X
        self.initial_y = self.traced_node.Y

    def CalculateValue(self):
        self.current_x =  self.traced_node.X
        self.current_y = self.traced_node.Y
        # self.value = (self.current_x-self.initial_x)**2 + (self.current_y-self.initial_y)**2
        self.value = (self.current_x-self.initial_x)**2
        # self.value = (self.current_y-self.initial_y)**2

    def CalculateGradient(self):
        self.shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        self.shape_gradient[0] = 2*(self.current_x-self.initial_x)
        # self.shape_gradient[1] = 2*(self.current_y-self.initial_y)

        return self.shape_gradient

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.body_model_part.Nodes}
        gradient[self.traced_node_id] = self.shape_gradient
        return gradient



# class FixPointResponseFunction(ResponseFunctionInterface):

#     def __init__(self, response_id, response_settings, model):
#         self.model = model
#         self.traced_node_id = response_settings["traced_node"].GetInt()
#         self.dim = response_settings["dim"].GetInt()

#     def Initialize(self):
#         self.body_model_part = self.model["MainModelPart.Body2D_Body"]
#         self.traced_node = self.model["MainModelPart"].GetNode(self.traced_node_id)
#         if self.dim == 0:
#             self.initial_value = self.traced_node.X
#         elif self.dim == 1:
#             self.initial_value = self.traced_node.Y

#     def CalculateValue(self):
#         if self.dim == 0:
#             self.current_value = self.traced_node.X
#         elif self.dim == 1:
#             self.current_value = self.traced_node.Y
#         self.value = (self.current_value-self.initial_value)**2

#     def CalculateGradient(self):
#         self.shape_gradient = KratosMultiphysics.Vector(3, 0.0)
#         self.shape_gradient[self.dim] = 2*(self.current_value-self.initial_value)

#         return self.shape_gradient

#     def GetValue(self):
#         return self.value

#     def GetNodalGradient(self, variable):
#         zero_vector = KratosMultiphysics.Vector(3, 0.0)
#         gradient ={node.Id : zero_vector for node in self.body_model_part.Nodes}
#         gradient[self.traced_node_id] = self.shape_gradient
#         return gradient