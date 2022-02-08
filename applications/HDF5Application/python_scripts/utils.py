"""HDF5 custom IO process utilities.

license: HDF5Application/license.txt
"""


__all__ = [
    "ParametersWrapper",
    "IsDistributed",
    "CreateOperationSettings",
    "OpenHDF5File"
]


# Core imports
import KratosMultiphysics

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.core.utils import ParametersWrapper


def CreateOperationSettings(operation_type, user_settings):
    """Return core settings for an operation type.

    See core.operations.
    """
    core_settings = ParametersWrapper("""
    {
        "operation_type": "%s"
    }
    """ % operation_type)
    for key in user_settings:
        core_settings[key] = user_settings[key]
    return core_settings


class OpenHDF5File(object):
    def __init__(self, file_parameters: KratosMultiphysics.Parameters, is_distributed: bool):
        """Context manager for HDF5 files."""
        file_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        if is_distributed:
            file_parameters["file_driver"].SetString("mpio")
            self.file = HDF5Application.HDF5FileParallel(file_parameters)
        else:
            self.file = HDF5Application.HDF5FileSerial(file_parameters)

    def __enter__(self) -> HDF5Application.HDF5File:
        return self.file

    def __exit__(self, exit_type, exit_value, exit_traceback) -> None:
        # HDF5::File has RAII, so this is the best we can do to close it
        self.file = None
        # TODO: handle exceptions

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "file_name" : "",
            "file_access_mode" : "read",
            "file_driver" : "sec2",
            "echo_level" : 0
        }""")