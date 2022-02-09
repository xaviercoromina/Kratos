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
import abc


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


# Resolve metaclass conflicts
if type(KratosMultiphysics.Process) != type(abc.ABC):
    class AbstractProcessMetaClass(type(abc.ABC), type(KratosMultiphysics.Process)):
        pass
else:
    class AbstractProcessMetaClass(type(KratosMultiphysics.Process)):
        pass


class CheckpointIOProcessBase(abc.ABC, KratosMultiphysics.Process, metaclass=AbstractProcessMetaClass):
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

    @abc.abstractmethod
    def Execute(self) -> None:
        """Pure virtual function for executing the io process."""
        pass

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "checkpoint_settings" {
                "prefix" : "",
                "check_mesh_consistency" : false,
                "file_settings" : {
                    "file_name" : "checkpoints/<model_part_name>_checkpoint_<step>.h5",
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        """Get parameters that can be used to create a user defined io process (HDF5Application.user_defined_io_process.Factory)."""
        parameters = KratosMultiphysics.Parameters()
        parameters.AddEmptyArray("Parameters")

        parameters["Parameters"].Append(KratosMultiphysics.Parameters())
        parameters["Parameters"][0].AddEmptyArray("list_of_operations")

        # Register operations
        for operation_type, parameter_generator in self._GetProcessMap().items():
            operation_parameters = parameter_generator()
            operation_parameters.AddString("operation_type", operation_type)
            parameters["Parameters"][0]["list_of_operations"].Append(operation_parameters)

        for process_parameters in parameters["Parameters"]:
            process_parameters.AddValue("io_settings", self.parameters["file_settings"])
            process_parameters.AddString("model_part_name", self.model_part.FullName())
            process_parameters.AddString("process_step", "")

        return parameters

    @abc.abstractmethod
    def _GetProcessMap(self) -> dict:
        """Pure virtual function for getting a map that associates parameter generators to operations."""
        return {}

    def _GetNodalSolutionStepDataParameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters(r"""{"list_of_variables" : []}""")
        variable_names = list(self.model_part.GetHistoricalVariablesNames())
        parameters["list_of_variables"].SetStringArray(variable_names)
        return parameters

    def _GetNodalDataNames(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters(r"""{"list_of_variables" : []}""")
        variable_names = list(self.model_part.GetNonHistoricalVariablesNames(
            self.model_part.Nodes,
            self.parameters["check_mesh_consistency"].GetBool()))
        parameters["list_of_variables"].SetStringArray(variable_names)
        return parameters

    def _GetElementDataParameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters(r"""{"list_of_variables" : []}""")
        variable_names = list(self.model_part.GetNonHistoricalVariablesNames(
            self.model_part.Elements,
            self.parameters["check_mesh_consistency"].GetBool()))
        parameters["list_of_variables"].SetStringArray(variable_names)
        return parameters

    def _GetConditionDataParameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters(r"""{"list_of_variables" : []}""")
        variable_names = list(self.model_part.GetNonHistoricalVariablesNames(
            self.model_part.Conditions,
            self.parameters["check_mesh_consistency"].GetBool()))
        # TODO: handle legacy variables
        #parameters["list_of_variables"].SetStringArray(variable_names)
        return parameters

    def _GetProcessInfoParameters(self) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters()

    def __ProcessFactory(self) -> KratosMultiphysics.Process:
        return UserDefinedIOProcessFactory(self._GetIOParameters(), self.model)


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
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    def _GetProcessMap(self) -> dict:
        return {
            "nodal_solution_step_data_input" :  self._GetNodalSolutionStepDataParameters,
            "nodal_data_value_input" :          self._GetNodalDataNames,
            "element_data_value_input" :        self._GetElementDataParameters,
            "condition_data_value_input" :      self._GetConditionDataParameters,
            "process_info_input" :              self._GetProcessInfoParameters
        }

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        parameters = super()._GetIOParameters()
        parameters["Parameters"][0]["process_step"].SetString("initialize_solution_step")

        # Set file access mode
        file_parameters = self.parameters["file_settings"].Clone()
        file_parameters.AddString("file_access_mode", "read_only")

        parameters["Parameters"][0].AddValue("io_settings", file_parameters)
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
                    "max_files_to_keep" : "unlimited",
                    "time_format" : "0.4f",
                    "echo_level" : 0
                }
            }
        }""")

    def _GetProcessMap(self) -> dict:
        return {
            "nodal_solution_step_data_output" : self._GetNodalSolutionStepDataParameters,
            "nodal_data_value_output" :         self._GetNodalDataNames,
            "element_data_value_output" :       self._GetElementDataParameters,
            "condition_data_value_output" :     self._GetConditionDataParameters,
            "process_info_output" :             self._GetProcessInfoParameters
        }

    def _GetIOParameters(self) -> KratosMultiphysics.Parameters:
        parameters = super()._GetIOParameters()
        parameters["Parameters"][0]["process_step"].SetString("finalize_solution_step")

        # Adjust file limit to accomodate buffer checkpoints
        file_parameters = self.parameters["file_settings"].Clone()
        if file_parameters["max_files_to_keep"].GetString() != "unlimited":
            file_parameters["max_files_to_keep"].SetString(str(int(file_parameters["max_files_to_keep"].GetString()) * self.model_part.GetBufferSize()))

        # Set file access mode
        file_parameters.AddString("file_access_mode", "truncate")

        parameters["Parameters"][0].AddValue("io_settings", file_parameters)
        return parameters


class Factory:
    """
    Create a wrapped solver that has the additional functionality of writing/loading checkpoints.

    """

    def __init__(self, solver_type: PythonSolver):
        # Resolve metaclass conflicts
        if type(solver_type) != type(abc.ABC):
            class AbstractSolverMetaClass(type(abc.ABC), type(solver_type)):
                pass
        else:
            class AbstractSolverMetaClass(type(solver_type)):
                pass

        class WrappedSolver(abc.ABC, solver_type, metaclass=AbstractSolverMetaClass):
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

            def AdvanceInTime(self, time: float) -> float:
                """
                Calls the base class' method first, then loads a checkpoint and resets TIME if necessary.

                Note: loading a checkpoint is only possible in AdvanceInTime because PythonSolver
                keeps track of TIME separately, based on the return value of this function. Override
                this method in derived classes if you need more control over when checkpoints are loaded.
                """
                super().AdvanceInTime(time)
                step_to_load = self._ShouldLoadCheckpoint()
                if step_to_load != self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]:
                    self.LoadCheckpoint(step_to_load)
                return self.GetMainModelPart().ProcessInfo[KratosMultiphysics.TIME]

            @RequiresInitialized("checkpoint_output_process")
            def FinalizeSolutionStep(self) -> None:
                super().FinalizeSolutionStep()
                if self._IsCheckpointOutputStep():
                    self.checkpoint_output_process.Execute()

            @RequiresInitialized("checkpoint_input_process")
            def LoadCheckpoint(self, target_step: int) -> None:
                """
                Restore the model part's state at the specified step from checkpoint files.

                Note: the model part is rolled back to the input step, and its buffer
                is overwritten accordingly. Every variable in the process info is
                overwritten except STEP, that's the responsibility of the caller.
                """
                original_step = self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]
                checkpoints = self.GetCheckpoints()

                # Fill buffers with data from the checkpoints
                step_begin = target_step - self.GetMainModelPart().GetBufferSize() + 1
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
                        self.GetMainModelPart().CloneSolutionStep() # TODO: ModelPart::CreateSolutionStep would suffice but throws an exception for now
                    self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP] = current_step
                    self.checkpoint_input_process.Execute()

                # Set process info
                self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP] = original_step

            def GetCheckpoints(self) -> list[dict]:
                """
                Return a list of dictionaries containing all file paths matching the
                checkpoint file pattern as well as all placeholders' values in each file
                (except '<model_part_name>', which is inferred directly from the solver's model part).

                Example:
                    pattern: '/<model_part_name>/<step>.h5'
                    matching files: /FluidModelPart/0.h5
                                    /FluidModelPart/10.h5
                    returns:
                            [
                                {'path' : pathlib.Path('/FluidModelPart/0.h5'), '<step>' : 0},
                                {'path' : pathlib.Path('/FluidModelPart/10.h5'), '<step>' : 10}
                            ]
                """
                # TODO: The HDF5Application incorrectly uses ModelPart::Name instead of the standard ModelPart::FullName
                file_pattern_absolute = str(pathlib.Path(self.checkpoint_parameters["checkpoint_settings"]["file_settings"]["file_name"].GetString()).absolute())
                file_pattern_absolute = file_pattern_absolute.replace("<model_part_name>", self.GetMainModelPart().Name)

                pattern = ModelPartPattern(file_pattern_absolute)
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

            def GetMainModelPart(self) -> KratosMultiphysics.ModelPart:
                if hasattr(super(), "GetMainModelPart"):
                    return super().GetMainModelPart()
                else:
                    return self.GetComputingModelPart().GetRootModelPart()

            @classmethod
            def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
                """Static virtual function for getting a complete set of input parameters with default values."""
                # Get solver parameters
                parameters = super().GetDefaultParameters()

                # Get io parameters
                io_parameters = CheckpointInputProcess.GetDefaultParameters()
                io_parameters.RecursivelyAddMissingParameters(CheckpointOutputProcess.GetDefaultParameters())
                io_parameters.RemoveValue("model_part_name") # Don't force 'model_part_name' because PythonSolver::GetComputingModelPart is used
                parameters.RecursivelyAddMissingParameters(io_parameters)

                # Add temporal control parameters
                parameters["checkpoint_settings"].AddEmptyValue("output_time_settings")
                parameters["checkpoint_settings"]["output_time_settings"].AddInt("step_frequency", 0)

                return parameters

            def _IsCheckpointOutputStep(self) -> bool:
                """Virtual function determining whether a checkpoint should be written at the end of the current step."""
                buffer_size = self.GetMainModelPart().GetBufferSize()
                step = self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]
                return self.checkpoint_frequency - buffer_size <= (step-1) % self.checkpoint_frequency

            @abc.abstractmethod
            def _ShouldLoadCheckpoint(self) -> int:
                """
                Pure virtual function for determining which step to revert to, if necessary.
                Return the current step if no checkpoint needs to be loaded.
                """
                return self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]

            def _CheckpointInputProcessFactory(self) -> CheckpointIOProcessBase:
                """Virtual function instantiating a checkpoint loader process."""
                parameters = self.checkpoint_parameters.Clone()
                parameters.AddString("model_part_name", self.GetMainModelPart().FullName())
                RecursivelyRemoveExtraParameters(parameters, CheckpointInputProcess.GetDefaultParameters())
                return CheckpointInputProcess(self.model, parameters)

            def _CheckpointOutputProcessFactory(self) -> CheckpointIOProcessBase:
                """Virtual function instantiating a checkpoint writer process."""
                parameters = self.checkpoint_parameters.Clone()
                parameters.AddString("model_part_name", self.GetMainModelPart().FullName())
                RecursivelyRemoveExtraParameters(parameters, CheckpointOutputProcess.GetDefaultParameters())
                return CheckpointOutputProcess(self.model, parameters)

            def __InitializeIOProcesses(self) -> None:
                """Instantiate and initialize checkpoint io processes."""
                self.checkpoint_input_process = self._CheckpointInputProcessFactory()
                self.checkpoint_input_process.ExecuteInitialize()

                self.checkpoint_output_process = self._CheckpointOutputProcessFactory()
                self.checkpoint_output_process.ExecuteInitialize()

        self.solver = WrappedSolver

    def Create(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        return self.Type(model, parameters)

    def __call__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        return self.Create(model, parameters)

    @property
    def Type(self) -> type:
        return self.solver