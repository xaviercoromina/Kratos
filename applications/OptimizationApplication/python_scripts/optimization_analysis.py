import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationProcessFactory
from KratosMultiphysics.OptimizationApplication.utilities.container_variable_data_holder_output_decorator import ContainerVariableDataHolderOutputDecorator

class OptimizationAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        default_parametetrs = Kratos.Parameters("""{
            "model_parts"       : [],
            "analyses"          : [],
            "responses"         : [],
            "controls"          : [],
            "output_processes"  : [],
            "algorithm_settings": {}
        }""")

        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.AddMissingParameters(default_parametetrs)

        self.optimization_info = OptimizationInfo(self.project_parameters["problem_data"]["echo_level"].GetInt())

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateResponses()
        self._CreateControlTechniques()

        super().__init__(model, project_parameters)

    def _CreateSolver(self):
        default_algorithm_settings = Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_algorithm_settings)

        return OptimizationProcessFactory(algorithm_settings["module"].GetString(), algorithm_settings["type"].GetString(), self.model, algorithm_settings, self.optimization_info, Algorithm)

    def _GetSimulationName(self):
        return "OptimizationAnalysis"

    def _GetListOfProcesses(self):
        if not hasattr(self, '_modified_list_of_processes'):
            self._modified_list_of_processes = super()._GetListOfProcesses()

            # adding all the optimization processes
            for optimization_processes_type in self._GetOptimizationProcessTypesOrder():
                for optimization_process in self.optimization_info.GetOptimizationProcesses(optimization_processes_type):
                    if not isinstance(optimization_process, ControlTransformationTechnique):
                        self._modified_list_of_processes.append(optimization_process)
                    else:
                        # adds the control
                        self._modified_list_of_processes.append(optimization_process.GetControl())

                        # adds the transformation techniques
                        for transformation_technique in optimization_process.GetTransformationTechniques():
                            self._modified_list_of_processes.append(transformation_technique)

                        # adds the control tranformation technique
                        self._modified_list_of_processes.append(optimization_process)

        return self._modified_list_of_processes

    def _GetListOfOutputProcesses(self):
        if not hasattr(self, '_modified__list_of_output_processes'):
            self._modified__list_of_output_processes = super()._GetListOfOutputProcesses()
            self._CreateContainerVariableDataHolderOutputDecorators()
            self._modified__list_of_output_processes.extend(self.optimization_info.GetOptimizationProcesses(ContainerVariableDataHolderOutputDecorator))

        return self._modified__list_of_output_processes

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self.optimization_info["step"])

    def KeepAdvancingSolutionLoop(self):
        return self.optimization_info["step"] < self.end_time and not self._GetSolver().IsConverged()

    def _CreateModelPartControllers(self):
        for model_part_controller_settings in self.project_parameters["model_parts"]:
            default_model_part_controller_settings = Kratos.Parameters("""{
                "name"    : "",
                "module"  : "KratosMultiphysics.OptimizationApplication.model_part_controllers",
                "type"    : "MdpaModelPartController",
                "settings": {}
            }""")
            model_part_controller_settings.ValidateAndAssignDefaults(default_model_part_controller_settings)
            routine: ModelPartController = OptimizationProcessFactory(model_part_controller_settings["module"].GetString(), model_part_controller_settings["type"].GetString(), self.model, model_part_controller_settings["settings"], self.optimization_info, ModelPartController)
            self.optimization_info.AddOptimizationProcess(ModelPartController, model_part_controller_settings["name"].GetString(), routine)

    def _CreateAnalyses(self):
        for analyses_settings in self.project_parameters["analyses"]:
            routine = ExecutionPolicyDecorator(self.model, analyses_settings)
            self.optimization_info.AddOptimizationProcess(ExecutionPolicyDecorator, routine.GetExecutionPolicyName(), routine)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "name"     : "",
            "module"   : "KratosMultiphysics.OptimizationApplication.responses",
            "type"     : "",
            "settings" : {}
        }""")
        for response_settings in self.project_parameters["responses"]:
            response_settings.ValidateAndAssignDefaults(default_settings)
            routine: ResponseFunction = OptimizationProcessFactory(response_settings["module"].GetString(), response_settings["type"].GetString(), self.model, response_settings["settings"], self.optimization_info, ResponseFunction)
            self.optimization_info.AddOptimizationProcess(ResponseFunction, response_settings["name"].GetString(), routine)

    def _CreateControlTechniques(self):
        for control_settings in self.project_parameters["controls"]:
            routine = ControlTransformationTechnique(self.model, control_settings, self.optimization_info)
            self.optimization_info.AddOptimizationProcess(ControlTransformationTechnique, routine.GetName(), routine)

    def _CreateContainerVariableDataHolderOutputDecorators(self):
        for output_process_settings in self.project_parameters["output_processes"]:
            self.optimization_info.AddOptimizationProcess(ContainerVariableDataHolderOutputDecorator, ContainerVariableDataHolderOutputDecorator(self.model, output_process_settings, self.optimization_info))

    def _GetOptimizationProcessTypesOrder(self):
        return [ModelPartController, ExecutionPolicyDecorator, ResponseFunction, ControlTransformationTechnique]
