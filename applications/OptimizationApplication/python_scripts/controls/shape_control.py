import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator

class ShapeControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names"            : [""],
            "control_update_variable_name": "ARRAY3_CONTROL_UPDATE",
            "mesh_moving_analysis_name"   : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]

        control_update_variable_name = parameters["control_update_variable_name"].GetString()
        control_update_variable_type = Kratos.KratosGlobals.GetVariableType(control_update_variable_name)
        if control_update_variable_type != "Array":
            raise RuntimeError(f"{control_update_variable_name} with {control_update_variable_type} type is not supported. Only supports array variables")

        self.control_update_variable = Kratos.KratosGlobals.GetVariable(control_update_variable_name)
        self.mesh_moving_execution_policy_wrapper: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["mesh_moving_analysis_name"].GetString())

    def UpdateControl(self, control_values: KratosOA.NodalContainerVariableDataHolder):
        historical_container = KratosOA.HistoricalContainerVariableDataHolder(control_values)
        historical_container.AssignDataToContainerVariable(Kratos.MESH_DISPLACEMENT)
        self.mesh_moving_execution_policy_wrapper.Execute()

    def GetModelParts(self) -> 'list[Kratos.ModelPart]':
        return self.model_parts

    def CreateContainerVariableDataHolder(self, model_part: Kratos.ModelPart) -> KratosOA.NodalContainerVariableDataHolder:
        return KratosOA.NodalContainerVariableDataHolder(model_part)

    def GetControlSensitivityVariable(self) -> any:
        return Kratos.SHAPE_SENSITIVITY

    def GetControlUpdateVariable(self) -> any:
        return self.control_update_variable

