
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestVtuOutput(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.Id)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

    def test_WriteMesh(self):
        vtu_output = Kratos.VtuOutput()
        vtu_output.SetEntities(self.model_part, Kratos.VtuOutput.ELEMENTS)
        vtu_output.AddHistoricalVariable(Kratos.PRESSURE)
        vtu_output.AddHistoricalVariable(Kratos.VELOCITY)
        vtu_output.Write("test")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()