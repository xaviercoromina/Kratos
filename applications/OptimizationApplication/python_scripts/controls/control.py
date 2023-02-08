import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion

class Control(Kratos.Process):
    def __init__(self):
        super().__init__()

    def CreateContainerVariableDataHolder(self, model_part: Kratos.ModelPart) -> ContainerVariableDataHolderUnion:
        raise NotImplementedError("Calling base class Control::CreateContainerVariableDataHolder method. Please implement it in the derrived class.")

    def GetControlUpdateVariable(self) -> any:
        raise NotImplementedError("Calling base class Control::GetControlUpdateVariable method. Please implement it in the derrived class.")

    def GetControlSensitivityVariable(self) -> any:
        raise NotImplementedError("Calling base class Control::GetControlSensitivityVariable method. Please implement it in the derrived class.")

    def GetModelParts(self) -> 'list[Kratos.ModelPart]':
        raise NotImplementedError("Calling base class Control::GetModelParts method. Please implement it in the derrived class.")

    def UpdateControl(self, control_data: ContainerVariableDataHolderUnion):
        raise NotImplementedError("Calling base class Control::UpdateControl method. Please implement it in the derrived class.")
