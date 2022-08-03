# Core imports
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
import KratosMultiphysics.HDF5Application.checkpoint.model_part_predicate as Predicates

# STD imports
import typing


def SetModelPartData(model_part: KratosMultiphysics.ModelPart, step: int = 0, path: int = 0, time: float = 0.0) -> None:
    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, path + step * time * (node.Id << 1)) # historical
        node[KratosMultiphysics.NODAL_H] = path + step * time * node.Id # non-historical


def MakeModel() -> "tuple[KratosMultiphysics.Model, KratosMultiphysics.ModelPart]":
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("test")
    model_part.SetBufferSize(2)
    model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
    model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
    model_part.ProcessInfo[HDF5Application.ANALYSIS_PATH] = 0
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

    # Nodes
    # 3-----4
    # |   / |
    # |  /  |
    # | /   |
    # 1-----2
    for i_node in range(4):
        model_part.CreateNewNode(i_node + 1, i_node % 2.0, i_node // 2.0, 0.0)

    # Elements
    # +-----+
    # | 2 / |
    # |  /  |
    # | / 1 |
    # +-----+
    model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], KratosMultiphysics.Properties(0))
    model_part.CreateNewElement("Element2D3N", 2, [1, 4, 3], KratosMultiphysics.Properties(0))

    # Conditions
    # 1-----1
    # |   / |
    # |  /  |
    # | /   |
    # +-----+
    model_part.CreateNewCondition("LineCondition2D2N", 1, [3, 4], KratosMultiphysics.Properties(0))

    SetModelPartData(model_part)
    return model, model_part


class TestModelPartConditions(KratosUnittest.TestCase):

    def test_ConstantPredicate(self) -> None:
        # Always true
        true_parameters = KratosMultiphysics.Parameters(R'{"value" : true}')
        true_predicate = Predicates.ConstantPredicate(true_parameters)
        self.CheckPredicate(true_predicate, lambda m: True, MakeModel()[1])

        # Always false
        false_parameters = KratosMultiphysics.Parameters(R'{"value" : false}')
        false_predicate = Predicates.ConstantPredicate(false_parameters)
        self.CheckPredicate(false_predicate, lambda m: False, MakeModel()[1])

    def test_StepIntervalPredicate(self) -> None:
        # Default parameters (should always return true)
        parameters = KratosMultiphysics.Parameters()
        predicate = Predicates.StepIntervalPredicate(parameters)
        self.CheckPredicate(predicate, lambda m: True, MakeModel()[1])

        # [0,5)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0,5]}')
        predicate = Predicates.StepIntervalPredicate(parameters)
        self.CheckPredicate(predicate,
                            lambda m: m.ProcessInfo[KratosMultiphysics.STEP] in (0,1,2,3,4),
                            MakeModel()[1])

        # [5,10)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [5,10]}')
        predicate = Predicates.StepIntervalPredicate(parameters)
        self.CheckPredicate(predicate,
                            lambda m: m.ProcessInfo[KratosMultiphysics.STEP] in (5,6,7,8,9),
                            MakeModel()[1])

    def test_TimeIntervalPredicate(self) -> None:
        step_and_time = [(index, 1.0 / (2 << (9-exponent))) for index, exponent in enumerate(range(10), 1)]

        # Default parameters (should always return true)
        parameters = KratosMultiphysics.Parameters()
        predicate = Predicates.TimeIntervalPredicate(parameters)
        self.CheckPredicate(predicate, lambda m: True, MakeModel()[1])

        # [0,0.5) - exactly representable boundaries
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0.0, 0.5]}')
        predicate = Predicates.TimeIntervalPredicate(parameters)
        self.CheckPredicate(predicate,
                            lambda m: self.IsInFloatInterval(0.0, m.ProcessInfo[KratosMultiphysics.TIME], 0.5),
                            MakeModel()[1],
                            step_and_time = step_and_time)

        # [0.5,1.0) - exactly representable boundaries
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0.5, 1.0]}')
        predicate = Predicates.TimeIntervalPredicate(parameters)
        self.CheckPredicate(predicate,
                            lambda m: self.IsInFloatInterval(0.5, m.ProcessInfo[KratosMultiphysics.TIME], 1.0),
                            MakeModel()[1],
                            step_and_time = step_and_time)

    def test_CyclicStepIntervalPredicate(self) -> None:
        # Default parameters (should always return true)
        parameters = KratosMultiphysics.Parameters()
        predicate = Predicates.CyclicStepIntervalPredicate(parameters)
        self.CheckPredicate(predicate, lambda m: True, MakeModel()[1])

        # mod 5 [0,3)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0,3], "cycle_size" : 5}')
        predicate = Predicates.CyclicStepIntervalPredicate(parameters)
        model_part = MakeModel()[1]

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 1
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 2
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 3
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 4
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 5
        predicate(model_part)
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 6
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 7
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 8
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 9
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 10
        self.assertTrue(predicate(model_part))

        # mod 5 [3,5)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [3,5], "cycle_size" : 5}')
        predicate = Predicates.CyclicStepIntervalPredicate(parameters)
        model_part = MakeModel()[1]

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 1
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 2
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 3
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 4
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 5
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 6
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 7
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 8
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 9
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.STEP] = 10
        self.assertFalse(predicate(model_part))

    def test_CyclicTimeIntervalPredicate(self) -> None:
        # Default parameters (should always return true)
        parameters = KratosMultiphysics.Parameters()
        predicate = Predicates.CyclicTimeIntervalPredicate(parameters)
        self.CheckPredicate(predicate, lambda m: True, MakeModel()[1])

        # mod 0.5 [0, 0.3)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0.0, 0.3], "cycle_size" : 0.5}')
        predicate = Predicates.CyclicTimeIntervalPredicate(parameters)
        model_part = MakeModel()[1]

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.1
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.2
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.3
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.4
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.5
        predicate(model_part)
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.6
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.7
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.8
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.9
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0
        self.assertTrue(predicate(model_part))

        # mod 0.5 [0.3, 0.5)
        parameters = KratosMultiphysics.Parameters(R'{"interval" : [0.3, 0.5], "cycle_size" : 0.5}')
        predicate = Predicates.CyclicTimeIntervalPredicate(parameters)
        model_part = MakeModel()[1]

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.1
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.2
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.3
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.4
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.5
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.6
        self.assertFalse(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.7
        self.assertFalse(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.8
        self.assertTrue(predicate(model_part))
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.9
        self.assertTrue(predicate(model_part))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0
        self.assertFalse(predicate(model_part))

    def CheckPredicate(self,
                       predicate: Predicates.ModelPartPredicate,
                       reference_predicate: typing.Callable,
                       model_part: KratosMultiphysics.ModelPart,
                       step_and_time: list = [(s, (1 if s%2 else -1) * (s<<1)) for s in range(20)]) -> None:
        for step, time in step_and_time:
            model_part.CloneTimeStep(time)
            model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            self.assertEqual(predicate(model_part), reference_predicate(model_part), (step, time))
            self.assertEqual(predicate(model_part), reference_predicate(model_part), (step, time))

    @staticmethod
    def IsInFloatInterval(lower: float, value: float, upper: float) -> bool:
        """@brief Check whether lower <= value < float assuming exactly representable boundaries and test value."""
        return lower <= value and value < upper


if __name__ == "__main__":
    KratosUnittest.main()
