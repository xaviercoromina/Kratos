import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver

from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor

from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper

class OptimizationSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "model_part_name": "optimization_model_part",
            "model_parts"    : [],
            "analyses"       : [],
            "responses"      : [],
            "controls"       : [],
            "algorithms"     : [],
            "echo_level"     : 0
        }""")

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        super().__init__(model, settings)

        # creates the optimization info data holder
        self.optimization_info = OptimizationInfo()
        self.__is_converged = False
        self.__list_of_algorithms: 'list[Algorithm]' = []
        self.__list_of_algorithm_properties: 'list[(list[str], Kratos.Parameters, Kratos.Parameters)]' = []

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateAlgorithms()

        # set the optimization info buffer size
        self.optimization_info.SetBufferSize(self.GetMinimumBufferSize())

    def GetOptimizationInfo(self):
        return self.optimization_info

    def GetMinimumBufferSize(self):
        buffer_size = -1000
        for algorithm in self.__list_of_algorithms:
            buffer_size = max(buffer_size, algorithm.GetMinimumBufferSize())
        return buffer_size

    def AddVariables(self):
        CallOnAll(self.__list_of_algorithms, Algorithm.AddVariables)

    def AddDofs(self):
        CallOnAll(self.__list_of_algorithms, Algorithm.AddDofs)

    def ImportModelPart(self):
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.ImportModelPart)

        # now we create other types because, responses and
        # controls can be used on submodel parts which are
        # only available after importing the whole model part
        self._CreateResponses()
        self._CreateControlTechniques()

        # now we can assign responses and control techniques to algorithms
        self._AssignAlgorithmProperties()

    def Check(self):
        CallOnAll(self.__list_of_algorithms, Algorithm.Check)

    def Initialize(self):
        # set the current step to 0
        self.optimization_info["step"] = 0

        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.Initialize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ExecutionPolicyWrapper), ExecutionPolicyWrapper.Initialize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ResponseFunction), ResponseFunction.Initialize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique), ControlTransformationTechnique.Initialize)
        CallOnAll(self.__list_of_algorithms, Algorithm.Initialize)

    def InitializeSolutionStep(self):
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.InitializeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ExecutionPolicyWrapper), ExecutionPolicyWrapper.InitializeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ResponseFunction), ResponseFunction.InitializeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique), ControlTransformationTechnique.InitializeSolutionStep)
        CallOnAll(self.__list_of_algorithms, Algorithm.InitializeSolutionStep)

    def SolveSolutionStep(self):
        CallOnAll(self.__list_of_algorithms, Algorithm.SolveSolutionStep)
        self.__is_converged = True
        for algorithm in self.__list_of_algorithms:
            self.__is_converged = self.__is_converged and algorithm.IsConverged()
        return self.__is_converged

    def IsConverged(self):
        return self.__is_converged

    def FinalizeSolutionStep(self):
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.FinalizeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ExecutionPolicyWrapper), ExecutionPolicyWrapper.FinalizeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ResponseFunction), ResponseFunction.FinalizeSolutionStep)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique), ControlTransformationTechnique.FinalizeSolutionStep)
        CallOnAll(self.__list_of_algorithms, Algorithm.FinalizeSolutionStep)

    def Finalize(self):
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.Finalize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ExecutionPolicyWrapper), ExecutionPolicyWrapper.Finalize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ResponseFunction), ResponseFunction.Finalize)
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique), ControlTransformationTechnique.Finalize)
        CallOnAll(self.__list_of_algorithms, Algorithm.Finalize)

    def GetComputingModelPart(self):
        model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(model_part_name):
            self.model.CreateModelPart(model_part_name)

        return self.model[model_part_name]

    def AdvanceInTime(self, _):
        self.optimization_info.AdvanceSolutionStep()
        self.optimization_info["step"] = self.optimization_info.GetSolutionStepData(1)["step"] + 1

        # now apply the controls
        if self.optimization_info["step"] > 1:
            CallOnAll(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique), ControlTransformationTechnique.ApplyControlUpdate)

        self.GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_info["step"]

        return self.optimization_info["step"]

    def _CreateModelPartControllers(self):
        for model_part_controller_settings in self.settings["model_parts"]:
            default_model_part_controller_settings = Kratos.Parameters("""{
                "name"    : "",
                "module"  : "KratosMultiphysics.OptimizationApplication.model_part_controllers",
                "type"    : "MdpaModelPartController",
                "settings": {}
            }""")
            model_part_controller_settings.ValidateAndAssignDefaults(default_model_part_controller_settings)
            routine: ModelPartController = Factory(model_part_controller_settings["module"].GetString(), model_part_controller_settings["type"].GetString(), self.model, model_part_controller_settings["settings"], self.optimization_info, ModelPartController)
            self.optimization_info.AddOptimizationRoutine(ModelPartController, model_part_controller_settings["name"].GetString(), routine)

    def _CreateAnalyses(self):
        for analyses_settings in self.settings["analyses"]:
            routine = ExecutionPolicyWrapper(self.model, analyses_settings)
            self.optimization_info.AddOptimizationRoutine(ExecutionPolicyWrapper, routine.GetName(), routine)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "name"     : "",
            "module"   : "KratosMultiphysics.OptimizationApplication.responses",
            "type"     : "",
            "settings" : {}
        }""")
        for response_settings in self.settings["responses"]:
            response_settings.ValidateAndAssignDefaults(default_settings)
            routine: ResponseFunction = Factory(response_settings["module"].GetString(), response_settings["type"].GetString(), self.model, response_settings["settings"], self.optimization_info, ResponseFunction)
            self.optimization_info.AddOptimizationRoutine(ResponseFunction, response_settings["name"].GetString(), routine)

    def _CreateControlTechniques(self):
        for control_settings in self.settings["controls"]:
            routine = ControlTransformationTechnique(self.model, control_settings, self.optimization_info)
            self.optimization_info.AddOptimizationRoutine(ControlTransformationTechnique, routine.GetName(), routine)

    def _CreateAlgorithms(self):
        default_settings = Kratos.Parameters("""{
            "name"         : "",
            "module"       : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"         : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "control_names": [],
            "objectives"   : [],
            "constraints"  : [],
            "settings"     : {}
        }""")
        for algorithm_settings in self.settings["algorithms"]:
            algorithm_settings.ValidateAndAssignDefaults(default_settings)
            algorithm: Algorithm = Factory(algorithm_settings["module"].GetString(), algorithm_settings["type"].GetString(), self.model, algorithm_settings["settings"], self.optimization_info, Algorithm)
            algorithm.SetName(algorithm_settings["name"].GetString())
            self.__list_of_algorithms.append(algorithm)
            self.__list_of_algorithm_properties.append(
                (
                    algorithm_settings["control_names"].GetStringArray(),
                    algorithm_settings["objectives"],
                    algorithm_settings["constraints"]
                ))

    def _AssignAlgorithmProperties(self):
        for algorithm, (control_names, objectives_settings, contraints_settings) in zip(self.__list_of_algorithms, self.__list_of_algorithm_properties):
            algorithm.SetControllers([self.optimization_info.GetOptimizationRoutine(ControlTransformationTechnique, control_name) for control_name in control_names])
            algorithm.SetObjectives([ObjectiveResponseFunctionImplementor(objective_settings, self.optimization_info) for objective_settings in objectives_settings])
            algorithm.SetConstraints([ConstraintResponseFunctionImplementor(constraint_settings, self.optimization_info) for constraint_settings in contraints_settings])