# Core imports
import KratosMultiphysics

# STD imports
import abc
#import sys
import platform


def ReplaceString(string: str, value, parameters: KratosMultiphysics.Parameters) -> None:
    """@brief Replace instances of the input string in the provided @ref Parameters with the specified value.
       @param string: string instance to replace.
       @param value: value to replace @a string with. Supported types: @a bool, @a int, @a float.
       @param parameters: @ref Parameters to modify.
       @throws @a TypeError if value is not an instance of any of the supported types.
    """
    if isinstance(value, bool):
        replace = lambda subparameters: subparameters.SetBool(value)
    elif isinstance(value, int):
        replace = lambda subparameters: subparameters.SetInt(value)
    elif isinstance(value, float):
        replace = lambda subparameters: subparameters.SetDouble(value)
    else:
        raise TypeError(f"Expecting a bool, int, or float for 'value' but got {type(value)}")

    for item in parameters:
        if item.IsString() and item.GetString() == string:
            replace(item)


class ModelPartPredicate(abc.ABC):
    """@brief Base class for model part predicates.
       @details The task of a model part predicate is to give a binary
                True/False answer based on the state of a @ref ModelPart.
       @details Example: check a model part whether @ref TIME is within
                a specified range.
       @note Any derived class must be constructed from a single
             @ref Parameters object.
    """

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        super().__init__()

    @abc.abstractmethod
    def __call__(self, model_part: KratosMultiphysics.ModelPart) -> bool:
        """@brief Perform the check and return a boolean answer."""
        pass

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters()


class ConstantPredicate(ModelPartPredicate):
    """@brief A predicate that always returns True or False based on what value it was constructed with."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        super().__init__(parameters)
        parameters.AddMissingParameters(self.GetDefaultParameters())
        self.__value = parameters["value"].GetBool()

    def __call__(self, model_part: KratosMultiphysics.ModelPart) -> bool:
        return self.__value

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "value" : true
        }""")


class IntervalPredicate(ModelPartPredicate):
    """@brief Base predicate providing an interface for checking whether a value is within the given range [begin,end).
       @details The interval can be specified with begin/end integers or "Begin"/"End" respectively.
       @note Default parameters: {
                 "interval" : ["Begin", "End"]
             }
    """

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        parameters.AddMissingParameters(self.GetDefaultParameters())
        self.__ParseBoundaries(parameters)
        super().__init__(parameters)
        self.parameters = parameters

    def __ParseBoundaries(self, parameters: KratosMultiphysics.Parameters) -> None:
        interval = parameters["interval"]
        if interval.IsArray():
            if interval.size() == 2:
                ReplaceString("Begin", self.min_begin, interval[0])
                ReplaceString("End", self.max_end, interval[1])
            else:
                raise ValueError(f"Invalid 'interval' size ({interval.size()}) in {interval}")
        else:
            raise ValueError(f"Expecting 'interval' to be an array of size 2, but got {interval}")

    @abc.abstractproperty
    def min_begin(self):
        """@brief Lowest possible begin of an interval."""
        return None

    @abc.abstractproperty
    def max_end(self):
        """@brief Highest possible end of an interval."""
        return None

    @abc.abstractmethod
    def _IsInInterval(self, value) -> bool:
        """@brief Check whether the input value is within the interval."""
        pass

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "interval" : ["Begin", "End"]
        }""")


class StepIntervalPredicate(IntervalPredicate):
    """@brief Predicate checking whether @STEP is within the given interval [begin, end)."""

    def __call__(self, model_part: KratosMultiphysics.ModelPart) -> bool:
        step = model_part.ProcessInfo[KratosMultiphysics.STEP]
        return self._IsInInterval(step)

    @property
    def min_begin(self) -> int:
        return 0

    @property
    def max_end(self) -> int:
        # Python automatically grows its int, to the point at which pybind can't
        # convert it to a C++ integer => we need a smaller value that
        # fits into the allocated bitset (and note that integers are always
        # signed in python).
        #return sys.maxsize
        bit_count = int(platform.architecture()[0][:-3])
        unsigned_bit_count = (bit_count - 1) // 2 - 1
        return 1 << unsigned_bit_count

    def _IsInInterval(self, value: int) -> bool:
        return self.parameters["interval"][0].GetInt() <= value and value < self.parameters["interval"][1].GetInt()


class TimeIntervalPredicate(IntervalPredicate):
    """@brief Predicate checking whether @TIME is within the given interval [begin, end)."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        super().__init__(parameters)
        self.__interval_utility = KratosMultiphysics.IntervalUtility(self.parameters)

    def __call__(self, model_part: KratosMultiphysics.ModelPart) -> bool:
        time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self._IsInInterval(time)

    @property
    def min_begin(self) -> float:
        return float("-inf")

    @property
    def max_end(self) -> float:
        return float("inf")

    def _IsInInterval(self, value: float) -> bool:
        return self.__interval_utility.IsInInterval(value)


class CyclicIntervalPredicate(IntervalPredicate):
    """@brief Base class for recurring predicates with intervals."""

    @abc.abstractproperty
    def cycle_size(self):
        """@brief Modulo."""
        return None

    def _IsInInterval(self, value) -> bool:
        return super()._IsInInterval(value % self.cycle_size)


class CyclicStepIntervalPredicate(CyclicIntervalPredicate, StepIntervalPredicate):
    """@brief Predicate checking whether @ref STEP % cycle_size is within the given interval [begin, end)."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        super().__init__(parameters)
        parameters.AddMissingParameters(self.GetDefaultParameters())
        cycle_size = parameters["cycle_size"]
        if cycle_size.IsInt():
            self.__cycle_size = cycle_size.GetInt()
        else:
            raise TypeError(f"Expecting 'cycle_size' as an integer but got {cycle_size}")

    @property
    def cycle_size(self) -> int:
        return self.__cycle_size

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "cycle_size" : 1,
            "interval" : ["Begin", "End"]
        }""")


class CyclicTimeIntervalPredicate(CyclicIntervalPredicate, TimeIntervalPredicate):
    """@brief Predicate checking whether @ref TIME % cycle_size is within the provided interval [begin, end)."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        super().__init__(parameters)
        if parameters.Has("cycle_size"):
            cycle_size = parameters["cycle_size"]
            if cycle_size.IsDouble():
                self.__cycle_size = cycle_size.GetDouble()
            else:
                raise TypeError(f"Expecting 'cycle_size' as a float but got {cycle_size}")
        else:
            raise ValueError(f"Expecting parameters with 'cycle_size' but got {parameters}")

    @property
    def cycle_size(self) -> float:
        return self.__cycle_size

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "cycle_size" : 1.0,
            "interval" : ["Begin", "End"]
        }""")