from importlib import import_module
from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

ContainerVariableDataHolderUnion = Union[
        KratosOA.HistoricalContainerVariableDataHolder,
        KratosOA.NodalContainerVariableDataHolder,
        KratosOA.ConditionContainerVariableDataHolder,
        KratosOA.ElementContainerVariableDataHolder,
        KratosOA.ConditionPropertiesContainerVariableDataHolder,
        KratosOA.ElementPropertiesContainerVariableDataHolder]

ContainerTypes = Union[
    Kratos.NodesArray,
    Kratos.ConditionsArray,
    Kratos.ElementsArray
]

def OptimizationProcessFactory(
    module_name: str,
    class_name: str,
    model: Kratos.Model,
    parameters: Kratos.Parameters,
    optimization_info: OptimizationInfo = None,
    required_object_type = Kratos.Process):

    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')
    full_module_name = f"{module_name}.{python_file_name}"

    module = import_module(full_module_name)
    if optimization_info is not None:
        retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)
    else:
        retrieved_object = getattr(module, class_name)(model, parameters)

    # check retrieved object is of the required type
    if not isinstance(retrieved_object, required_object_type):
        raise RuntimeError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not derrived from the \"{required_object_type.__name__}\".")
    else:
        return retrieved_object