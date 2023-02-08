
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
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
            "echo_level"     : 1
        }""")

        process = KratosOA.EntitySpecificPropertiesProcess(self.model, parameters)

        process.ExecuteInitializeSolutionStep()
        # change the entity specific properties
        for element in self.model_part.Elements:
            element.Properties[Kratos.DENSITY] = element.Id
            element.Properties[Kratos.PRESSURE] = element.Id + 1

        for condition in self.model_part.Conditions:
            condition.Properties[Kratos.DENSITY] = condition.Id + 2
            condition.Properties[Kratos.PRESSURE] = condition.Id + 3

        # after the first iteration, following variables should be in the element/condition data value container with the
        # values from the properties
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Id)
            self.assertEqual(element.Properties[Kratos.PRESSURE], element.Id + 1)

        for condition in self.model_part.Conditions:
            self.assertEqual(condition.Properties[Kratos.DENSITY], condition.Id + 2)
            self.assertEqual(condition.Properties[Kratos.PRESSURE], condition.Id + 3)


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()