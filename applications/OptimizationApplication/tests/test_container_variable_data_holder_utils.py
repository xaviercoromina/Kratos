import math

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerVariableDataHolderUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    def test_ContainerVariableDataHolderNormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 15)

    def test_ContainerVariableDataHolderNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        l2_norm = 0.0
        for node in self.model_part.Nodes:
            l2_norm += node.GetSolutionStepValue(Kratos.PRESSURE) ** 2
        l2_norm = math.sqrt(l2_norm)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(a), l2_norm)

        l2_norm = 0.0
        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            l2_norm += v[0] ** 2 + v[1] ** 2 + v[2] ** 2
        l2_norm = math.sqrt(l2_norm)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(a), l2_norm)

    def test_ContainerVariableDataHolderEntityMaxNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), math.sqrt(15**2 + 14**2 + 13**2))

    def test_ContainerVariableDataHolderInnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, b), 890)

    def test_CollectiveVariableDataHolderNormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(collective_1), max(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), KratosOA.ContainerVariableDataHolderUtils.NormInf(b)))

    def test_CollectiveVariableDataHolderNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(collective_1), math.sqrt(KratosOA.ContainerVariableDataHolderUtils.NormL2(a)**2 + KratosOA.ContainerVariableDataHolderUtils.NormL2(b)**2))


    def test_CollectiveVariableDataHolderInnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(collective_1, collective_1), KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, a) + KratosOA.ContainerVariableDataHolderUtils.InnerProduct(b, b))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()