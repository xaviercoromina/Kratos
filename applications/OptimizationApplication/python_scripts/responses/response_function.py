from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos

class ResponseFunction(Kratos.Process, ABC):
    def __init__(self):
        super().__init__()

    @abstractmethod
    def Check(self):
        pass

    @abstractmethod
    def CalculateValue(self) -> float:
        pass

    @abstractmethod
    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):
        pass

