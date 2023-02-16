import KratosMultiphysics as KM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintInternalForcesNormProcess(model, settings["Parameters"])

class PrintInternalForcesNormProcess(KM.Process):


    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "model_part_name"           : "",
            "output_file_name"          : "trial.txt",
            "print_interval"            : 0.0
        }""")
        return default_parameters

    def __init__(self, model, settings):
        """The constructor of the WaveHeightOutputProcess"""

        KM.Process.__init__(self)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        model_part = model.GetModelPart(settings["model_part_name"].GetString())

        settings.RemoveValue("model_part_name")
        self.print_process = PFEM.PrintInternalForcesNormProcess(model_part, settings)

    def ExecuteFinalizeSolutionStep(self):
        self.print_process.ExecuteFinalizeSolutionStep()
