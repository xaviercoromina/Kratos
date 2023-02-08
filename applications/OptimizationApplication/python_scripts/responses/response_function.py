import KratosMultiphysics as Kratos

class ResponseFunction(Kratos.Process):
    def __init__(self):
        super().__init__()

    def Check(self):
        raise RuntimeError("Calling ResponseFunction::Check. This should be implemented in the derrived class.")

    def CalculateValue(self) -> float:
        raise RuntimeError("Calling ResponseFunction::CalculateValue. This should be implemented in the derrived class.")

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):
        raise RuntimeError("Calling ResponseFunction::CalculateSensitivity. This should be implemented in the derrived class.")

