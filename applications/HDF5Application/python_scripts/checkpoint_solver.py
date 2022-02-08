__all__ = ["Factory"]

# Core imports
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# HDF5 imports
from KratosMultiphysics.HDF5Application.user_defined_io_process import Factory as UserDefinedIOProcessFactory
from KratosMultiphysics.HDF5Application import ModelPartPattern

# STL imports
import pathlib
import functools


def RequiresInitialized(*object_names: str):
    """Member function decorator that checks whether input names exist in the class and have values other than None."""
    def Decorator(function):
        @functools.wraps(function)
        def WrappedFunction(this, *args, **kwargs):
            for object_name in object_names:
                if not hasattr(this, object_name) or getattr(this, object_name) == None:
                    raise RuntimeError("'{object_name}' is required to be initialized by '{function_name}::{class_name}'".format(
                                        object_name = object_name,
                                        function_name = function.__name__,
                                        class_name = type(this).__name__))
            return function(this, *args, **kwargs)
        return WrappedFunction
    return Decorator


def RecursivelyRemoveExtraParameters(subject: KratosMultiphysics.Parameters, reference: KratosMultiphysics.Parameters) -> None:
    for key in subject.keys():
        if not reference.Has(key):
            subject.RemoveValue(key)
        else:
            if reference[key].IsSubParameter():
                if subject[key].IsSubParameter():
                    RecursivelyRemoveExtraParameters(subject[key], reference[key])
                else:
                    raise ValueError("Type mismatch for key: '{}'".format(key))


class CheckpointIOProcessBase(KratosMultiphysics.Process):
    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model = model
        self.model_part = self.model.GetModelPart(parameters["model_part_name"].GetString())
        self.parameters = parameters["checkpoint_settings"]

        # Placeholder for initialization later
        self.process = None

    def ExecuteInitialize(self) -> None:
        self.process = self.__ProcessFactory()
        self.process.ExecuteInitialize()
        self.process.ExecuteBeforeSolutionLoop()

    def Execute(self) -> None:
        raise RuntimeError("Attempt to call a pure virtual member")

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "checkpoint_settings" {
                "prefix" : "",
                "check_mesh_consistency" : false,
                "file_settings" : {
                    "file_name" : "checkpoints/<model_part_name>_checkpoint_<step>.h5",
                    "file_access_mode" : "truncate",
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters()
        parameters.AddEmptyArray("Parameters")

        parameters["Parameters"].Append(KratosMultiphysics.Parameters())
        parameters["Parameters"][0].AddEmptyArray("list_of_operations")

        # Register operations
        variable_parameters = self.__GetVariablesForIOParameters()
        for operation_type, variable_type in self._GetProcessMap().items():
            operation_parameters = KratosMultiphysics.Parameters()
            operation_parameters.AddString("operation_type", operation_type)
            operation_parameters.AddValue("prefix", self.parameters["prefix"])
            if variable_parameters.Has(variable_type):
                operation_parameters.AddValue("list_of_variables", variable_parameters[variable_type]["list_of_variables"])
            parameters["Parameters"][0]["list_of_operations"].Append(operation_parameters)

        for process_parameters in parameters["Parameters"]:
            process_parameters.AddValue("io_settings", self.parameters["file_settings"])
            process_parameters.AddString("model_part_name", self.model_part.Name)
            process_parameters.AddString("process_step", "")

        return parameters

    @staticmethod
    def _GetProcessMap() -> dict:
        raise RuntimeError("Attempt to call a pure virtual member")

    def __ProcessFactory(self) -> KratosMultiphysics.Process:
        return UserDefinedIOProcessFactory(self._GetIOParameters(), self.model)

    def __GetVariablesForIOParameters(self) -> KratosMultiphysics.Parameters:
        """Collect all variables present in the nodes of the model part."""
        parameters = KratosMultiphysics.Parameters("""{
            "nodal_solution_step_data" : {},
            "nodal_data" : {},
            "element_data" : {},
            "condition_data" : {}
        }""")

        for item in parameters:
            item.AddEmptyArray("list_of_variables")

        check_mesh_consistency = self.parameters["check_mesh_consistency"].GetBool()

        parameters["nodal_solution_step_data"]["list_of_variables"].SetStringArray(
            list(self.model_part.GetHistoricalVariablesNames()))
        parameters["nodal_data"]["list_of_variables"].SetStringArray(
            list(self.model_part.GetNonHistoricalVariablesNames(self.model_part.Nodes, check_mesh_consistency)))
        parameters["element_data"]["list_of_variables"].SetStringArray(
            list(self.model_part.GetNonHistoricalVariablesNames(self.model_part.Elements, check_mesh_consistency)))
        #parameters["condition_data"]["list_of_variables"].SetStringArray(
        #    list(self.model_part.GetNonHistoricalVariablesNames(self.model_part.Conditions, check_mesh_consistency)))

        return parameters


class CheckpointInputProcess(CheckpointIOProcessBase):
    @RequiresInitialized("process")
    def Execute(self) -> None:
        self.process.ExecuteInitializeSolutionStep()

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "checkpoint_settings" : {
                "prefix" : "",
                "check_mesh_consistency" : false,
                "file_settings" : {
                    "file_name" : "checkpoints/<model_part_name>_checkpoint_<step>.h5",
                    "file_access_mode" : "read_only",
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    @staticmethod
    def _GetProcessMap() -> dict:
        return {
            "nodal_solution_step_data_input" : "nodal_solution_step_data",
            "nodal_data_value_input" : "nodal_data",
            "element_data_value_input" : "element_data",
            "condition_data_value_input" : "condition_data",
            "process_info_input" : ""
        }

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        parameters = super()._GetIOParameters()
        parameters["Parameters"][0]["process_step"].SetString("initialize_solution_step")
        parameters["Parameters"][0].AddValue("io_settings", self.parameters["file_settings"])
        return parameters


class CheckpointOutputProcess(CheckpointIOProcessBase):
    @RequiresInitialized("process")
    def Execute(self) -> None:
        self.process.ExecuteFinalizeSolutionStep()

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "checkpoint_settings" : {
                "prefix" : "",
                "check_mesh_consistency" : false,
                "file_settings" : {
                    "file_name" : "checkpoints/<model_part_name>_checkpoint_<step>.h5",
                    "file_access_mode" : "truncate",
                    "max_files_to_keep" : "unlimited",
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    @staticmethod
    def _GetProcessMap() -> dict:
        return {
            "nodal_solution_step_data_output" : "nodal_solution_step_data",
            "nodal_data_value_output" : "nodal_data",
            "element_data_value_output" : "element_data",
            "condition_data_value_output" : "condition_data",
            "process_info_output" : ""
        }

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        parameters = super()._GetIOParameters()
        parameters["Parameters"][0]["process_step"].SetString("finalize_solution_step")

        # Adjust file limit to accomodate buffer checkpoints
        file_parameters = self.parameters["file_settings"].Clone()
        if file_parameters["max_files_to_keep"].GetString() != "unlimited":
            file_parameters["max_files_to_keep"].SetString(str(int(file_parameters["max_files_to_keep"].GetString()) * self.model_part.GetBufferSize()))
        parameters["Parameters"][0].AddValue("io_settings", file_parameters)

        return parameters


class Factory:
    def __init__(self, solver_type: PythonSolver):
        class WrappedSolver(solver_type):
            __name__ = solver_type.__name__ + "WithCheckpoints"
            __qualname__ = solver_type.__qualname__ + "WithCheckpoints"

            def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters, *arguments, **keyword_arguments):
                """

                Notes:
                 - the output file pattern must include a "<step>" placeholder to identify which
                   step the checkpoint belongs to.
                 - steps are assumed to begin at 0.
                 - the buffer size is assumed (and required) to be constant.
                 - the number of kept files depends on the model part's buffer size. Checkpoints
                   need to restore the state of the buffer as well, so the true number of files is:
                   max_files_to_keep * (buffer_size - 1) + 1
                   The "+1" comes from static information written to the first checkpoint that needs
                   to be kept.
                """
                # Set legacy attributes (occasionally used in kratos)
                self.__class__.__name__ = self.__name__

                # Extract checkpoint parameters
                parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
                self.checkpoint_parameters = parameters.Clone()
                RecursivelyRemoveExtraParameters(self.checkpoint_parameters, self.GetDefaultParameters())

                # Initialize base
                super().__init__(model, parameters, *arguments, **keyword_arguments)

                # Parse input parameters
                self.model = model
                self.model_part = model.GetModelPart(parameters["model_part_name"].GetString())
                self.checkpoint_frequency = self.checkpoint_parameters["checkpoint_settings"]["output_time_settings"]["step_frequency"].GetInt()
                self.checkpoint_parameters["checkpoint_settings"].RemoveValue("output_time_settings")

                # Check arguments
                file_pattern = self.checkpoint_parameters["checkpoint_settings"]["file_settings"]["file_name"].GetString()
                if not "<step>" in file_pattern:
                    raise ValueError("Output file pattern ('{}') must include a '<step>'".format(file_pattern))

                # Placeholders for members instantiated later
                self.checkpoint_input_process = None
                self.checkpoint_output_process = None

            def Initialize(self) -> None:
                """
                Collect all variables in the model part and initialize the IO processes.
                Note: initialization must be executed on a populated model part.
                """
                super().Initialize()
                self.__InitializeIOProcesses()

            #@RequiresInitialized("checkpoint_output_process")
            #def InitializeSolutionStep(self) -> None:
            #    super().InitializeSolutionStep()
            #    self.checkpoint_output_process.ExecuteInitializeSolutionStep()

            @RequiresInitialized("checkpoint_output_process")
            def FinalizeSolutionStep(self) -> None:
                super().FinalizeSolutionStep()
                if self.IsCheckpointOutputStep():
                    self.checkpoint_output_process.Execute()

            #@RequiresInitialized("checkpoint_output_process")
            #def Finalize(self) -> None:
            #    super().Finalize()
            #    self.checkpoint_output_process.ExecuteFinalize()

            @RequiresInitialized("checkpoint_input_process")
            def LoadCheckpoint(self, target_step: int) -> None:
                """
                Restore the model part's state at the specified step from checkpoint files.

                Note: the model part is rolled back to the input step, and its buffer
                is overwritten accordingly. Every variable in the process info is
                overwritten except STEP, that's the responsibility of the caller.
                """
                original_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                checkpoints = self.GetCheckpoints()

                # Fill buffers with data from the checkpoints
                step_begin = target_step - self.model_part.GetBufferSize() + 1
                step_end = target_step + 1
                for current_step in range(step_begin, step_end):
                    # Get the relevant checkpoint
                    checkpoint_path = next((checkpoint["path"] for checkpoint in checkpoints if checkpoint["<step>"]==current_step), False)
                    if not checkpoint_path:
                        raise FileNotFoundError("No checkpoint found for step {}. Detected checkpoint files:\n{}".format(
                            current_step,
                            '\n'.join(path for path in (checkpoint["path"] for checkpoint in checkpoints))))

                    # Cycle the buffer and load data from the checkpoint
                    if current_step != step_begin: # the buffer need not be cycled on the first step
                        self.model_part.CloneSolutionStep() # TODO: ModelPart::CreateSolutionStep would suffice but throws an exception for now
                    self.model_part.ProcessInfo[KratosMultiphysics.STEP] = current_step
                    self.checkpoint_input_process.Execute()

                # Set process info
                self.model_part.ProcessInfo[KratosMultiphysics.STEP] = original_step

            def GetCheckpoints(self) -> list[dict]:
                file_pattern_absolute = pathlib.Path(self.checkpoint_parameters["checkpoint_settings"]["file_settings"]["file_name"].GetString()).absolute()
                pattern = ModelPartPattern(str(file_pattern_absolute))
                checkpoints = []

                # Collect files and extract their placeholders' values
                for path in pattern.Glob():
                    matches = pattern.Match(path)
                    checkpoint = {"path" : pathlib.Path(path)}
                    for key, value_type in (("<model_part_name>", str), ("<step>", int), ("<time>", float)):
                        values = matches.get(key, [])
                        if values:
                            checkpoint[key] = value_type(values[0]) # assume all placeholders share the same value

                    checkpoints.append(checkpoint)

                return sorted(checkpoints, key=lambda item: item["<step>"])

            def IsCheckpointOutputStep(self) -> bool:
                buffer_size = self.model_part.GetBufferSize()
                step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                return self.checkpoint_frequency - buffer_size <= (step-1) % self.checkpoint_frequency

            def __InitializeIOProcesses(self) -> None:
                output_parameters = self.checkpoint_parameters.Clone()
                RecursivelyRemoveExtraParameters(output_parameters, CheckpointOutputProcess.GetDefaultParameters())
                output_parameters["checkpoint_settings"]["file_settings"]["file_access_mode"].SetString("truncate")
                self.checkpoint_output_process = CheckpointOutputProcess(self.model, output_parameters)
                self.checkpoint_output_process.ExecuteInitialize()

                input_parameters = self.checkpoint_parameters.Clone()
                RecursivelyRemoveExtraParameters(input_parameters, CheckpointInputProcess.GetDefaultParameters())
                self.checkpoint_input_process = CheckpointInputProcess(self.model, input_parameters)
                self.checkpoint_input_process.ExecuteInitialize()

            @classmethod
            def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
                parameters = super().GetDefaultParameters()
                parameters.RecursivelyAddMissingParameters(CheckpointInputProcess.GetDefaultParameters())
                parameters.RecursivelyAddMissingParameters(CheckpointOutputProcess.GetDefaultParameters())
                parameters["checkpoint_settings"].AddEmptyValue("output_time_settings")
                parameters["checkpoint_settings"]["output_time_settings"].AddInt("step_frequency", 0)
                return parameters

        self.solver = WrappedSolver

    def Create(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        return self.Type(model, parameters)

    def __call__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        return self.Create(model, parameters)

    @property
    def Type(self) -> type:
        return self.solver