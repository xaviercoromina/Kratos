from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# used union types for type hinting in python and for automated documentation generation
SupportedSensitivityFieldVariableTypes = Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]

ContainerExpressionTypes = Union[
                                Kratos.ContainerExpression.HistoricalExpression,
                                Kratos.ContainerExpression.NodalNonHistoricalExpression,
                                Kratos.ContainerExpression.ConditionNonHistoricalExpression,
                                Kratos.ContainerExpression.ElementNonHistoricalExpression,
                                KratosOA.ContainerExpression.ConditionPropertiesExpression,
                                KratosOA.ContainerExpression.ElementPropertiesExpression]