
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class TestOptimizationInfo(kratos_unittest.TestCase):
    class TestRoutine(Kratos.Process):
        def __init__(self):
            super().__init__()

    def test_AdvanceSolutionStepGet(self):
        optimization_info = OptimizationInfo()
        optimization_info.SetBufferSize(3)

        optimization_info["step"] = 1
        self.assertEqual(optimization_info["step"], 1)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 2
        self.assertEqual(optimization_info["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 1)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 3
        self.assertEqual(optimization_info["step"], 3)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(2)["step"], 1)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 4
        self.assertEqual(optimization_info["step"], 4)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 3)
        self.assertEqual(optimization_info.GetSolutionStepData(2)["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(3)["step"], 4)

    def test_AdvanceSolutionStepSet(self):
        optimization_info = OptimizationInfo()
        optimization_info.SetBufferSize(3)

        optimization_info["step"] = 1
        self.assertEqual(optimization_info["step"], 1)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 2
        optimization_info.GetSolutionStepData(1)["step"] = 10
        self.assertEqual(optimization_info["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 10)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 3
        optimization_info.GetSolutionStepData(2)["step"] = 15
        self.assertEqual(optimization_info["step"], 3)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(2)["step"], 15)

        optimization_info.AdvanceSolutionStep()
        optimization_info["step"] = 4
        self.assertEqual(optimization_info["step"], 4)
        self.assertEqual(optimization_info.GetSolutionStepData(1)["step"], 3)
        self.assertEqual(optimization_info.GetSolutionStepData(2)["step"], 2)
        self.assertEqual(optimization_info.GetSolutionStepData(3)["step"], 4)

    def test_OptimizationProcesss(self):
        optimization_info = OptimizationInfo()
        optimization_info.SetBufferSize(1)
        temp = TestOptimizationInfo.TestRoutine()
        optimization_info.AddOptimizationProcess(TestOptimizationInfo.TestRoutine, "temp", temp)

        self.assertTrue(optimization_info.HasOptimizationProcess(TestOptimizationInfo.TestRoutine, "temp"))
        self.assertTrue(optimization_info.HasOptimizationProcessType(TestOptimizationInfo.TestRoutine))
        self.assertEqual(temp, optimization_info.GetOptimizationProcess(TestOptimizationInfo.TestRoutine, "temp"))
        self.assertEqual([temp], optimization_info.GetOptimizationProcesses(TestOptimizationInfo.TestRoutine))


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()