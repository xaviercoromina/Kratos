
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestEntitySpecificPropertiesProcess(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        properties = cls.model_part.CreateNewProperties(1)

        properties[Kratos.DENSITY] = 10.0
        properties[Kratos.DIAMETER] = 15.0
        properties[Kratos.PRESSURE] = -5.0
        properties[Kratos.VELOCITY] = Kratos.Array3([1, 2, 3])

        number_of_entities = 10
        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewNode(i, 0, 0, 0)

        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewCondition("LineCondition2D2N", i, [(i % number_of_entities) + 1, ((i + 1) % number_of_entities) + 1], properties)

        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewElement("Element2D3N", i, [(i % number_of_entities) + 1, ((i + 1) % number_of_entities) + 1, ((i + 2) % number_of_entities) + 1], properties)


    def test_EntitySpecificPropertiesProcess(self):
        parameters = Kratos.Parameters("""{
            "model_part_name": "test",
            "container_type" : "all",
            "echo_level"     : 0
        }""")

        process = KratosOA.EntitySpecificPropertiesProcess(self.model, parameters)

        process.ExecuteInitializeSolutionStep()

        for condition in self.model_part.Conditions:
            condition.Properties[Kratos.DENSITY] = condition.Id + 1
            condition.Properties[Kratos.PRESSURE] = condition.Id + 2

        for element in self.model_part.Elements:
            element.Properties[Kratos.DENSITY] = element.Id + 3
            element.Properties[Kratos.PRESSURE] = element.Id + 4

        for condition in self.model_part.Conditions:
            self.assertEqual(condition.Properties[Kratos.DENSITY], condition.Id + 1)
            self.assertEqual(condition.Properties[Kratos.PRESSURE], condition.Id + 2)
            self.assertEqual(condition.Properties[Kratos.DIAMETER], 15.0)
            self.assertVectorAlmostEqual(condition.Properties[Kratos.VELOCITY], Kratos.Array3([1, 2, 3]), 9)

        # now checking for updated properties
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Id + 3)
            self.assertEqual(element.Properties[Kratos.PRESSURE], element.Id + 4)
            self.assertEqual(element.Properties[Kratos.DIAMETER], 15.0)
            self.assertVectorAlmostEqual(element.Properties[Kratos.VELOCITY], Kratos.Array3([1, 2, 3]), 9)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()