"""HDF5 custom IO process utilities.

license: HDF5Application/license.txt
"""


__all__ = [
    "ParametersWrapper",
    "IsDistributed",
    "CreateOperationSettings",
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
