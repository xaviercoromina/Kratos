__all__ = ["Factory"]

# Core imports
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.python_solver import PythonSolver

# HDF5 imports
from KratosMultiphysics.HDF5Application.user_defined_io_process import Factory as UserDefinedIOProcessFactory
from KratosMultiphysics.HDF5Application import ModelPartPattern

# STL imports
import pathlib
import functools
import abc
import typing


def RequiresInitialized(*object_names: str) -> typing.Callable:
    """@brief Member function decorator that checks whether input names exist in the class and have values other than @c None."""
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
    """@brief Recursively remove all entries from @p subject that aren't in @p reference."""
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
    """ @brief Base class for managing checkpoint input/output operations.

        @details This class provides a flexible interface for generating and managing a user defined
                 io process (via @ref KratosMultiphysics.HDF5Application.user_defined_io_process.Factory).
                 Public methods (new or overridden):
                 - @ref CheckpointIOProcessBase.__init__
                 - @ref CheckpointIOProcessBase.Execute
                 Pure virtual methods:
                 - @ref CheckpointIOProcessBase.Execute
                 - @ref CheckpointIOProcessBase._GetProcessMap
                 @notes - Derived classes must override @c Execute and @c _GetProcessMap.
                        - Operations can be defined by overriding @c _GetProcessMap, which associates parameter generator functions to
                          operation names they create parameters for. See @ref CheckpointInputProcess and @ref CheckpointOutputProcess
                          for examples.
                        - As the logic of when to load/write checkpoints tends to be more complex, the process
                          gets carried out must be delegated to @c Execute rather than @c ExecuteInitialize, @c ExecuteFinalize, etc.
                          Deciding on when to call @c Execute and thus performing the IO process is the user's task.
                        - @c _GetIOParameters can be overridden to control what parameters get passed on to the HDF5 factory,
                          but the output must be compatible with @ref KratosMultiphysics.HDF5Application.user_defined_io_process.Factory.
    """

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model = model
        self.model_part = self.model.GetModelPart(parameters["model_part_name"].GetString())
        self.parameters = parameters["checkpoint_settings"].Clone()

        # Placeholder for initialization later
        self.process = None

    def ExecuteInitialize(self) -> None:
        self.process = self.__ProcessFactory()
        self.process.ExecuteInitialize()
        self.process.ExecuteBeforeSolutionLoop()

    @abc.abstractmethod
    def Execute(self) -> None:
        """@brief Pure virtual function for executing the io process."""
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
    """@brief Default input process for a checkpoint solver."""

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


class CheckpointOutputProcess(KratosMultiphysics.OutputProcess, CheckpointIOProcessBase):

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        """@brief Defer construction to the base class @ref CheckpointIOProcessBase"""
        KratosMultiphysics.OutputProcess.__init__(self)
        CheckpointIOProcessBase.__init__(self, model, parameters)

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

        # Set file access mode
        file_parameters = self.parameters["file_settings"]
        file_parameters.AddString("max_files_to_keep", "unlimited")
        file_parameters.AddString("file_access_mode", "truncate")

        parameters["Parameters"][0].AddValue("io_settings", file_parameters)
        return parameters


class Factory:
    """@brief Create a wrapped solver that has the additional functionality of writing/loading checkpoints.
    @details The factory exposes the class of the wrapped solver via the Type property, and can
             instantiate it by calling Create or __call__.
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
            """ @brief A wrapped PythonSolver with additional checkpoint writing/loading functionality.
                @details - The underlying solver must manage <b>EXACTLY ONE</b> root model part that can be
                           accessed via @c GetComputingModelPart().GetRootModelPart(). This requirement
                           is not checked at runtime and satisfying it is the responsibility of the user.
                         - The output file pattern must include a "<step>" placeholder to identify which
                           step the checkpoint belongs to.
                         - Steps are assumed to begin at 1.
                         - Derived classes must implement @c _ShouldLoadCheckpoint, that determines which checkpoint
                           to load if necessary. By default, no checkpoint is loaded if the returned step is equal
                           to the current step of the model part.
                         - The buffer size is assumed (and required) to be constant, though this requirement is not checked.
                         - The number of kept files depends on the model part's buffer size. Checkpoints
                           need to restore the state of the buffer as well, so the true number of files is
                           \f( max_files_to_keep \cdot buffer_size \f)
                         - By default, checkpoints are written during the first few steps while the buffer gets
                           filled up, and then, beginning at that step, periodically at every @c checkpoint_frequency
                           (with additional buffer checkpoints written when necessary; if the buffer size matches
                           or exceeds @c checkpoint_frequency, a checkpoint is written at each step).
                           If you need a different logic, override @c _IsCheckpointOutputStep.
                         - Historical and non-historical variables of nodes, elements, and conditions in the
                           solver's model part are <b>detected at solver initialization time</b>. All of the mentioned variables
                           are written/loaded to/from checkpoints along with the model part's process info.
                           If you need to customize what data gets stored in checkpoints, derive new classes
                           from @ref CheckpointIOProcessBase and override @c _CheckpointInputProcessFactory and
                           @c _CheckpointOutputProcessFactory accordingly. Depending on whether you need to touch
                           the input parameter structure, you may need to override @c GetDefaultParameters as well.
                @note - @ref Initialize must be called before invoking solution loop methods.
            """

            __name__ = solver_type.__name__ + "WithCheckpoints"
            __qualname__ = solver_type.__qualname__ + "WithCheckpoints"

            def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters, *arguments, **keyword_arguments):
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
                """ @brief Collect all variables in the @ref ModelPart and initialize the IO processes.
                @warning Initialization must be executed on a populated @ref ModelPart.
                """
                super().Initialize()
                self.__InitializeIOProcesses()

            def AdvanceInTime(self, time: float) -> float:
                """ @brief Calls the base class' method first, then loads a checkpoint and resets @c TIME if necessary.
                @note Loading a checkpoint is only possible in @c AdvanceInTime because @ref PythonSolver
                      keeps track of @c TIME separately, based on the return value of this function. Override
                      this method in derived classes if you need more control over when checkpoints are loaded.
                """
                super().AdvanceInTime(time)
                step_to_load = self._ShouldLoadCheckpoint()
                if step_to_load != self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]:
                    self.LoadCheckpoint(step_to_load)
                return self.GetMainModelPart().ProcessInfo[KratosMultiphysics.TIME]

            @RequiresInitialized("checkpoint_output_process")
            def FinalizeSolutionStep(self) -> None:
                """ @brief Write checkpoints if necessary, and delete obsolete ones."""
                super().FinalizeSolutionStep()
                if self._IsCheckpointOutputStep():
                    # Write checkpoint
                    self.checkpoint_output_process.Execute()

                    # Delete obsolete checkpoints
                    file_limit = self.checkpoint_parameters["checkpoint_settings"]["max_files_to_keep"].GetString()
                    if file_limit != "unlimited":
                        try:
                            # Adjust file limit to accomodate buffer checkpoints
                            max_files = int(file_limit) * self.GetComputingModelPart().GetBufferSize()
                        except Exception as exception:
                            raise ValueError("Invalid value for 'max_files_to_keep' ('{}')! Options are 'unlimited' or positive integers (in string format)".format(file_limit))

                        checkpoints = self.GetCheckpoints()
                        if max_files < len(checkpoints):
                            for checkpoint in checkpoints[:len(checkpoints) - max_files]:
                                DeleteFileIfExisting(str(checkpoint["path"]))

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
                            '\n'.join(str(path) for path in (checkpoint["path"] for checkpoint in checkpoints))))

                    # Cycle the buffer and load data from the checkpoint
                    if current_step != step_begin: # the buffer need not be cycled on the first step
                        self.GetMainModelPart().CloneSolutionStep() # TODO: ModelPart::CreateSolutionStep would suffice but throws an exception for now
                    self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP] = current_step
                    self.checkpoint_input_process.Execute()

                # Set process info
                self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP] = original_step

            def GetCheckpoints(self) -> list:
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

                # Add file limit (hdf5 controllers doesn't handle this correctly all the time)
                parameters["checkpoint_settings"].AddString("max_files_to_keep", "unlimited")

                return parameters

            def _IsCheckpointOutputStep(self) -> bool:
                """ @brief Virtual function determining whether a checkpoint should be written at the end of the current step.
                @details By default check"""
                buffer_size = self.GetMainModelPart().GetBufferSize()
                step = self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]
                return self.checkpoint_frequency - buffer_size <= (step-1-buffer_size) % self.checkpoint_frequency

            @abc.abstractmethod
            def _ShouldLoadCheckpoint(self) -> int:
                """ @brief Pure virtual function for determining which step to load a checkpoint from, if necessary.
                @return Step index whose checkpoint is to be loaded, or the current step index if no loading is necessary.
                """
                return self.GetMainModelPart().ProcessInfo[KratosMultiphysics.STEP]

            def _CheckpointInputProcessFactory(self) -> CheckpointIOProcessBase:
                """ @brief Virtual function instantiating a checkpoint loader process."""
                parameters = self.checkpoint_parameters.Clone()
                parameters.AddString("model_part_name", self.GetMainModelPart().FullName())
                RecursivelyRemoveExtraParameters(parameters, CheckpointInputProcess.GetDefaultParameters())
                return CheckpointInputProcess(self.model, parameters)

            def _CheckpointOutputProcessFactory(self) -> CheckpointIOProcessBase:
                """ @brief Virtual function instantiating a checkpoint writer process."""
                parameters = self.checkpoint_parameters.Clone()
                parameters.AddString("model_part_name", self.GetMainModelPart().FullName())
                RecursivelyRemoveExtraParameters(parameters, CheckpointOutputProcess.GetDefaultParameters())
                return CheckpointOutputProcess(self.model, parameters)

            def __InitializeIOProcesses(self) -> None:
                """ @brief Instantiate and initialize checkpoint io processes."""
                self.checkpoint_input_process = self._CheckpointInputProcessFactory()
                self.checkpoint_input_process.ExecuteInitialize()

                self.checkpoint_output_process = self._CheckpointOutputProcessFactory()
                self.checkpoint_output_process.ExecuteInitialize()

        self.solver = WrappedSolver

    def Create(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        """ @brief Construct an instance of the augmented solver."""
        return self.Type(model, parameters)

    def __call__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> PythonSolver:
        """ @brief Construct an instance of the augmented solver."""
        return self.Create(model, parameters)

    @property
    def Type(self) -> type:
        return self.solver
