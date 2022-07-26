# Core imports
import KratosMultiphysics

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5
from KratosMultiphysics.HDF5Application.utils import OpenHDF5File
import KratosMultiphysics.HDF5Application.core.operations.model_part as Operations
from KratosMultiphysics.HDF5Application import ModelPartPattern

# Core imports
import pathlib
import abc
import typing


class Snapshot(abc.ABCMeta):
    """@brief Class representing a snapshot of a @ref ModelPart state.
       @details A snapshot is uniquely defined by its path ID and step index
                for a specific analysis. The path ID indicates how many times
                the solution loop jumped back and continued from an earlier @ref Checkpoint,
                while the step index counts the number of steps since the analysis
                began, disregarding steps that branched off the current analysis path.
       @note Specialized for keeping data in memory or on disk.
    """

    def __init__(self, path_id: int, step: int):
        self.__path_id = path_id
        self.__step = step

    @abc.abstractmethod
    def Load(self, model_part: KratosMultiphysics.ModelPart) -> None:
        """@brief Load data from this snapshot to the specified model part."""
        pass

    @abc.abstractmethod
    def Write(self, model_part: KratosMultiphysics.ModelPart) -> None:
        """@brief Write data from the current state of the specified model part to the snapshot."""
        pass

    @staticmethod
    def GetSolutionPath(snapshots: list) -> list:
        """@brief Pick snapshots from the provided list that are part of the solution path.
           @param snapshots: list of snapshots of the analysis.
           @return A sorted list of snapshots that make up the solution path.
           @details A path is assembled backtracking from the last snapshot, recreating the
                    solution path iff the input list contains the solution path. Otherwise
                    the assembled path is the one that has a dead-end at the last snapshot."""
        solution_path = []

        # Assemble the reversed solution path
        if snapshots:
            snapshots = sorted(snapshots)[::-1]
            solution_path.append(snapshots.pop(0))

            while snapshots:
                last = solution_path[-1]
                current = snapshots.pop(0)

                if last.path_id == current.path_id: # last and current snapshots are on the same branch
                    solution_path.append(current)
                    continue
                else: # the last snapshot opened a new branch
                    # Looking for the snapshot with branch ID and
                    # step index strictly LOWER than those of the last snapshot.
                    if current.step < last.step:
                        solution_path.append(current)

        return solution_path[::-1]

    @property
    def path_id(self) -> int:
        return self.__path_id

    @property
    def step(self) -> int:
        return self.__step

    def __lt__(self, other: "Snapshot") -> bool:
        if self.__path_id == other.__path_id:
            return self.__step < other.__step
        else:
            return self.__path_id < other.__path_id

    def __gt__(self, other: "Snapshot") -> bool:
        if self.__path_id == other.__path_id:
            return self.__step > other.__step
        else:
            return self.__path_id > other.__path_id

    def __eq__(self, other: "Snapshot") -> bool:
        return not (self < other) and not (self > other)

    def __ne__(self, other: "Snapshot") -> bool:
        return not (self == other)

    def __str__(self) -> str:
        return f"Snapshot on path {self.path_id} and step {self.step}"

    def __repr__(self) -> str:
        return self.__str__()


class SnapshotInMemory(Snapshot):
    """@brief Class representing a snapshot of a @ref ModelPart state in memory (stores a deep copy of the model part).
       @todo Implement in-memory snapshots if necessary (@matekelemen).
    """
    pass


class SnapshotIOBase(abc.ABCMeta):
    """@brief Base class with common functionality to writing/loading snapshots to/from disk."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        self.__parameters = parameters
        self.__parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

    def __call__(self, model_part: KratosMultiphysics.ModelPart) -> None:
        """@brief Execute all defined IO operations."""
        with OpenHDF5File(self.__parameters["io_settings"], model_part.IsDistributed()) as file:
            for operation in self._GetOperations(self):
                operation(model_part, file)

    def ReadPathID(self) -> int:
        with OpenHDF5File(self.__GetInputParameters(), KratosMultiphysics.IsDistributedRun()) as file:
            return file.ReadIntAttribute(self.__parameters["operation_settings"]["prefix"].GetString(), "path_id")

    def WritePathID(self, id: int, file: HDF5.File = None) -> None:
        if file == None:
            with OpenHDF5File(self.__GetOutputParameters(), KratosMultiphysics.IsDistributedRun()) as file:
                self.__WritePathID(id, file)
        else:
            self.__WritePathID(id, file)

    def ReadStep(self) -> int:
        """@brief Read the step index from the snapshot file.
           @note The step index is read from the process info, which requires instantiating
                 a @ref Model and a @ModelPart, then reading the entire @ref ProcessInfo only
                 to get @ref STEP. This could be improved but would probably involve introducing
                 some redundancy.
        """
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("temporary")
        with OpenHDF5File(self.__GetInputParameters(), KratosMultiphysics.IsDistributedRun()) as file:
            Operations.ReadProcessInfo(self.__parameters["operation_settings"])(model_part, file)
        return model_part.ProcessInfo[KratosMultiphysics.STEP]

    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        return self.__parameters.Copy()

    @staticmethod
    def _ExtractNodalSolutionStepDataNames(model_part: KratosMultiphysics.ModelPart) -> "list[str]":
        return list(model_part.GetHistoricalVariableNames())

    @staticmethod
    def _ExtractNodalDataNames(model_part: KratosMultiphysics.ModelPart, check_mesh_consistency: bool = False) -> "list[str]":
        return list(model_part.GetNonHistoricalVariablesNames(model_part.Nodes,
                                                              check_mesh_consistency))

    @staticmethod
    def _ExtractElementDataNames(model_part: KratosMultiphysics.ModelPart, check_mesh_consistency: bool = False) -> "list[str]":
        return list(model_part.GetNonHistoricalVariablesNames(model_part.Elements,
                                                              check_mesh_consistency))

    @staticmethod
    def _ExtractConditionDataNames(model_part: KratosMultiphysics.ModelPart, check_mesh_consistency: bool = False) -> "list[str]":
        return list(model_part.GetNonHistoricalVariablesNames(model_part.Conditions,
                                                              check_mesh_consistency))

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters(R"""{
            "operation_settings" : {
                "prefix" : "/",
                "time_format" : "0.4f"
            }
        }""")
        parameters.AddValue("io_settings", cls.GetDefaultIOParameters())
        return parameters

    @abc.abstractstaticmethod
    def GetDefaultIOParameters() -> KratosMultiphysics.Parameters:
        pass

    @abc.abstractmethod
    def _GetOperations(self, model_part: KratosMultiphysics.ModelPart) -> typing.Iterable:
        """@brief Get all operations to be performed on the input model part."""
        pass

    def __WritePathID(self, id: int, file: HDF5.File) -> None:
        file.WriteIntAttribute(self.__parameters["operation_settings"]["prefix"].GetString(), "path_id", id)

    def __GetInputParameters(self) -> KratosMultiphysics.Parameters:
        """@brief Get IO parameters for reading a file regardless of whether the class is meant for reading or writing."""
        io_parameters = self.GetDefaultIOParameters()
        io_parameters["file_access_mode"].SetString("read_only")
        return io_parameters

    def __GetOutputParameters(self) -> KratosMultiphysics.Parameters:
        """@brief Get IO parameters for writing to a file regardless of whether the class is meant for reading or writing."""
        io_parameters = self.GetDefaultIOParameters()
        io_parameters["file_access_mode"].SetString("read_write")
        return io_parameters


class DefaultSnapshotOutput(SnapshotIOBase):
    """@brief Output class for writing most data in the model part to an HDF5 snapshot.
       @details Data written: - nodal solution step data
                              - nodal data value
                              - condition data value
                              - process info
    """

    @staticmethod
    def GetDefaultIOParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
        {
            "file_name" : "checkpoints/<model_part_name>_snapshot_<step>.h5",
            "file_access_mode" : "read_write"
        }
        """)

    def _GetOperations(self, model_part: KratosMultiphysics.ModelPart) -> list:
        operations = []

        # Variables
        for operation, variable_names in ((Operations.NodalSolutionStepDataOutput, self._ExtractNodalSolutionStepDataNames(model_part)),
                                          (Operations.NodalDataValueOutput, self._ExtractNodalDataNames(model_part)),
                                          (Operations.ElementDataValueOutput, self._ExtractElementDataNames(model_part)),
                                          (Operations.ConditionDataValueOutput), self._ExtractConditionDataNames(model_part)):
            parameters = self.parameters["operation_settings"]
            parameters.AddStringArray("list_of_variables", variable_names)
            operations.append(operation(parameters))

        # ProcessInfo
        operations.append(Operations.ProcessInfoOutput(self.parameters["operation_settings"]))

        return operations


class DefaultSnapshotInput(SnapshotIOBase):
    """@brief Input class for reading most data from an HDF5 snapshot to a model part.
       @details Data read: - nodal solution step data
                           - nodal data value
                           - condition data value
                           - process info
    """

    @staticmethod
    def GetDefaultIOParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
        {
            "file_name" : "checkpoints/<model_part_name>_snapshot_<step>.h5",
            "file_access_mode" : "read_only"
        }
        """)

    def _GetOperations(self, model_part: KratosMultiphysics.ModelPart) -> list:
        operations = []

        # Variables
        for operation, variable_names in ((Operations.NodalSolutionStepDataInput, self._ExtractNodalSolutionStepDataNames(model_part)),
                                          (Operations.NodalDataValueInput, self._ExtractNodalDataNames(model_part)),
                                          (Operations.ElementDataValueInput, self._ExtractElementDataNames(model_part)),
                                          (Operations.ConditionDataValueInput), self._ExtractConditionDataNames(model_part)):
            parameters = self.parameters["operation_settings"]
            parameters.AddStringArray("list_of_variables", variable_names)
            operations.append(operation(parameters))

        # ProcessInfo
        operations.append(Operations.ProcessInfoInput(self.parameters["operation_settings"]))

        return operations


class SnapshotOnDisk(Snapshot):
    """@brief Class representing a snapshot of a @ref ModelPart state and its associated output file."""

    def __init__(self,
                 path_id: int,
                 step: int,
                 file_path: pathlib.Path,
                 input: SnapshotIOBase = DefaultSnapshotInput(),
                 output: SnapshotIOBase = DefaultSnapshotOutput()):
        super().__init__(path_id, step)
        self.__file_path = file_path
        self.__input = input
        self.__output = output

    def Write(self, model_part: KratosMultiphysics.ModelPart) -> None:
        if self.__file_path.exists():
            if self.__file_path.is_dir():
                raise FileNotFoundError(f"{self.__file_path} is a directory")
            else:
                raise FileExistsError(f"File exists: {self.__file_path}")
        self.__output(model_part)

    def Load(self, model_part: KratosMultiphysics.ModelPart) -> None:
        if self.__file_path.is_file():
            self.__input(model_part)
        elif self.__file_path.is_dir():
            raise FileExistsError(f"{self.__file_path} is a directory")
        else:
            raise FileNotFoundError(f"File not found: {self.__file_path}")

    @staticmethod
    def FromFile(file_path: pathlib.Path,
                 input: SnapshotIOBase = DefaultSnapshotInput(),
                 output: SnapshotIOBase = DefaultSnapshotOutput()) -> "SnapshotOnDisk":
        """@brief Construct a @ref Snapshot instance from a snapshot file.
           @param file_path: path to a snapshot file to parse.
           @param input: IO instance responsible for reading snapshot files
           @param output: IO instance responsible for writing snapshot files
        """
        if file_path.is_file():
            path_id = input.ReadPathID()
            step = input.ReadStep()
            return SnapshotOnDisk(path_id, step, file_path, input, output)
        elif file_path.is_dir():
            raise FileExistsError(f"{file_path} is a directory")
        else:
            raise FileNotFoundError(f"File not found: {file_path}")

    @staticmethod
    def Collect(pattern: str,
                input: SnapshotIOBase = DefaultSnapshotInput(),
                output: SnapshotIOBase = DefaultSnapshotOutput()) -> list:
        """@brief Find and read all snapshot files that match the provided file name pattern.
           @param pattern: the file name pattern compatible with @ref ModelPartPattern to search for.
           @param input: IO instance responsible for reading snapshot files
           @param output: IO instance responsible for writing snapshot files
           @return a list of @ref SnapsotOnDisk loaded from discovered snapshot files, sorted in
                   ascending order (comparison is performed lexicographically over {path_id, step_index}).
        """
        return sorted([SnapshotIOBase.FromFile(file_path, input = input, output = output) for file_path in ModelPartPattern(pattern).Glob()])

    @property
    def file_path(self) -> pathlib.Path:
        return self.__file_path

    def __str__(self) -> str:
        return f"{super().__str__()} in {self.file_path}"
