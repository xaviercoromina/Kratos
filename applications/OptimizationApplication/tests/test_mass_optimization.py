
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestMassOptimization(kratos_unittest.TestCase):
    def test_MassOptimization(self):
        model = Kratos.Model()
        with kratos_unittest.WorkFolderScope("mass_optimization", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("mass_optimization", __file__):
            DeleteFileIfExisting("Structure.time")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()