# Core imports
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5Application
import KratosMultiphysics.HDF5Application.checkpoint.snapshot as Snapshots

# STD imports
import pathlib


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


class TestSnapshotOnDisk(KratosUnittest.TestCase):

    def setUp(self) -> None:
        DeleteDirectoryIfExisting(str(self.test_directory))
        (self.test_directory / "checkpoints").mkdir(parents = True, exist_ok = False)

    def tearDown(self) -> None:
        return
        DeleteDirectoryIfExisting(str(self.test_directory))

    def test_ReadWrite(self) -> None:
        input_parameters = Snapshots.DefaultSnapshotInput.GetDefaultParameters()
        output_parameters = Snapshots.DefaultSnapshotOutput.GetDefaultParameters()
        for parameters in (input_parameters, output_parameters):
            parameters["io_settings"]["file_name"].SetString(str(self.test_directory / parameters["io_settings"]["file_name"].GetString()))

        model, source_model_part = MakeModel()
        source_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        source_model_part.ProcessInfo[HDF5Application.ANALYSIS_PATH] = 0
        snapshot = Snapshots.SnapshotOnDisk(0, 0, input_parameters, output_parameters)
        snapshot.Write(source_model_part)

        target_model_part = model.CreateModelPart("read")
        snapshot.Load(target_model_part)
        print(source_model_part)
        print(target_model_part)

    @property
    def test_directory(self) -> pathlib.Path:
        return pathlib.Path("test_snapshot_on_disk")

if __name__ == "__main__":
    KratosUnittest.main()
