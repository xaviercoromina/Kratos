import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator
import os


class CompressibleNavierStokesSymbolicGeneratorTest(KratosUnittest.TestCase):

    @classmethod
    def _GetParameters(cls, geometry_name):
        return KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "geometry": "{geometry}",
                "template_filename" : "test_compressible_symbolic_{geometry}.template",
                "output_filename"   : "test_compressible_symbolic_{geometry}.result",
                "echo_level" : 0
            },
            "reference": "test_compressible_symbolic_{geometry}.reference"
        }
        """.replace("{geometry}", geometry_name)
        )


    def _RunTest(self, geometry_name):
        parameters = self._GetParameters(geometry_name)

        try:
            os.remove(parameters["Parameters"]["output_filename"].GetString())
        except FileNotFoundError:
            pass

        with KratosUnittest.WorkFolderScope(".", __file__):
            generator = CompressibleNavierStokesSymbolicGenerator(parameters["Parameters"])
            generator.Generate()
            generator.Write()

            with open(parameters["Parameters"]["output_filename"].GetString(), "r") as f:
                result = f.readlines()

            with open(parameters["reference"].GetString(), "r") as f:
                reference = f.readlines()

        self.assertEqual(result, reference)


    def testGeneratorTetra(self):
        self._RunTest("tetrahedron")

    def testGeneratorTriangle(self):
        self._RunTest("triangle")


if __name__ == '__main__':
    KratosUnittest.main()