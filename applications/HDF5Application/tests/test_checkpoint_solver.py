# Core imports
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics import KratosUnittest

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.checkpoint_solver import Factory as CheckpointSolverFactory


class DummySolver(CheckpointSolverFactory(PythonSolver).Type):
    def GetDofsList(self) -> list:
        return [KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.VELOCITY]

    def ImportModelPart(self) -> None:
        self.model_part = self.model.CreateModelPart("MainModelPart")

        self.model_part.SetBufferSize(self.GetMinimumBufferSize())

        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_X)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        for index, (x, y) in enumerate(((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))):
            self.model_part.CreateNewNode(index + 1, x, y, 0.0)
        self.model_part.CreateNewElement("Element2D4N", 1, (1, 2, 3, 4), self.model_part.GetProperties()[1])
        self.model_part.CreateNewCondition("PointCondition3D1N", 1, [1], self.model_part.GetProperties()[1])

    def GetMinimumBufferSize(self) -> int:
        return 3

    def ExportModelPart(self) -> None:
        pass

    def AdvanceInTime(self, time: float) -> float:
        if self.GetComputingModelPart().GetBufferSize() < self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
            time += 1.0
        self.model_part.CloneTimeStep(time)
        self.model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        return super().AdvanceInTime(time)

    def SolveSolutionStep(self) -> bool:
        time = round(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, time)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [time, time << 2, time << 3])
        return True

    def GetComputingModelPart(self) -> KratosMultiphysics.ModelPart:
        return self.model_part

    def _ShouldLoadCheckpoint(self) -> int:
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        offset_step = step - self.GetComputingModelPart().GetBufferSize()
        time = round(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        if 0 < offset_step and offset_step - offset_step // self.checkpoint_frequency < time:
            return step - self.checkpoint_frequency
        else:
            return step


class DummyAnalysis(AnalysisStage):
    def _CreateSolver(self):
        parameters = DummySolver.GetDefaultParameters()
        parameters["checkpoint_settings"]["output_time_settings"]["step_frequency"].SetInt(5)
        return DummySolver(self.model, parameters)

class TestCheckpointSolver(KratosUnittest.TestCase):
    def setUp(self) -> None:
        DeleteDirectoryIfExisting("checkpoints")

    def tearDown(self) -> None:
        DeleteDirectoryIfExisting("checkpoints")

    def test_dummy_solver(self) -> None:
        model = KratosMultiphysics.Model()
        parameters = KratosMultiphysics.Parameters("""{
            "problem_data" : {
                "echo_level" : 0,
                "parallel_type" : "serial",
                "start_time" : 0.0,
                "end_time" : 20.0
            }
        }""")
        analysis = DummyAnalysis(model, parameters)
        analysis.Run()


if __name__ == "__main__":
    KratosUnittest.main()