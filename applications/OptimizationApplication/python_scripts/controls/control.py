from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion


class Control(OptimizationRoutine, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

    @abstractmethod
    def CreateContainerVariableDataHolder(self, model_part: Kratos.ModelPart) -> ContainerVariableDataHolderUnion:
        """Returns a new container variable data holder object for model_part

        This returns a ContainerVariableDataHolder object of correct type for this control
        with the model_part.

        Returns:
            Container data holder with the correct type
        """
        pass

    @abstractmethod
    def GetControlUpdateVariable(self) -> any:
        """Returns the control update variable

        This method returns the update value holding variable for the control.
        Eg: For DENSITY control variable, one can use DENSITY_UPDATE

        Returns:
            any: Control update variable
        """
        pass

    @abstractmethod
    def GetControlSensitivityVariable(self) -> any:
        """Returns the sensitivity computation variable

        This method returns the variable on which the computed sensitivities
        for a respective response is stored.
        Eg: For DENSITY control variable, DENSITY_SENSITIVITY can be used.

        Returns:
            any: Returns the sensitivity storage variable
        """
        pass

    @abstractmethod
    def GetModelParts(self) -> 'list[Kratos.ModelPart]':
        """_summary_

        Returns:
            list[Kratos.ModelPart]: _description_
        """
        pass

    @abstractmethod
    def UpdateControl(self, control_data: ContainerVariableDataHolderUnion):
        """Updates the corresponding controls with given control values

        Args:
            control_values (Kratos.Vector): Given updated control values
        """
        pass