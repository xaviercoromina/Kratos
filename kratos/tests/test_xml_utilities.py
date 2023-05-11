
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestXmlUtilities(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))

        for condition in cls.model_part.Conditions:
            id = condition.Id
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))

        for element in cls.model_part.Elements:
            id = element.Id
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))

    def test_XmlNodes(self):
        xml_nodes = Kratos.XML_Utilities.XmlNodes([node for node in self.model_part.Nodes])
        # print(xml_nodes.Print())

    def test_XmlConditions(self):
        # first generate the list of nodes corresponding to conditions
        list_of_nodes = []
        for condition in self.model_part.Conditions:
            for node in condition.GetGeometry():
                list_of_nodes.append(node)
        xml_nodes = Kratos.XML_Utilities.XmlNodes(list(set(list_of_nodes)))
        xml_conditions = Kratos.XML_Utilities.XmlConditions(xml_nodes, [condition for condition in self.model_part.Conditions])
        print(xml_conditions.Print())


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()