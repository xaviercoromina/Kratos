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

    def test_ProductWithEntityMatrix(self):
        number_of_nodes = self.model_part.NumberOfNodes()

        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        m = Kratos.Matrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                m[i, j] = (i + 1) * (j + 1)

        b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(b, m, a)
        b.AssignDataToContainerVariable(Kratos.DENSITY)

        for i, node_b in enumerate(b.GetContainer()):
            v = 0
            for j, node_a in enumerate(a.GetContainer()):
                v += m[i, j] * node_a.GetSolutionStepValue(Kratos.PRESSURE)
            self.assertEqual(v, node_b.GetSolutionStepValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrixSparse(self):
        number_of_nodes = self.model_part.NumberOfNodes()

        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        # first build the normal matrix
        dense_m = Kratos.Matrix(
            [
                [5, 0, 0, 2, 0, 0, 0, 0, 0, 2],
                [0, 0, 3, 2, 3, 0, 0, 0, 0, 0],
                [0, 2, 0, 2, 0, 2, 0, 6, 0, 0],
                [0, 0, 0, 2, 0, 4, 0, 5, 0, 0],
                [0, 0, 2, 0, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 3, 3, 0],
                [0, 2, 0, 0, 3, 0, 4, 0, 0, 0],
                [3, 0, 0, 2, 0, 0, 0, 2, 0, 2],
                [9, 0, 0, 0, 0, 7, 0, 7, 0, 5],
                [0, 0, 7, 2, 0, 0, 0, 0, 0, 2]
            ]
        )

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        dense_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

        sparse_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

    def test_Transpose(self):
        number_of_nodes = self.model_part.NumberOfNodes()

        # first build the normal matrix
        dense_m = Kratos.Matrix(
            [
                [5, 0, 0, 2, 0, 0, 0, 0, 0],
                [0, 0, 3, 2, 3, 0, 0, 0, 0],
                [0, 2, 0, 2, 0, 2, 0, 6, 0],
                [0, 0, 0, 2, 0, 4, 0, 5, 0],
                [0, 0, 2, 0, 0, 2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 3, 3],
                [0, 2, 0, 0, 3, 0, 4, 0, 0],
                [3, 0, 0, 2, 0, 0, 0, 2, 0],
                [9, 0, 0, 0, 0, 7, 0, 7, 0],
                [0, 0, 7, 2, 0, 0, 0, 0, 0]
            ]
        )
        transpose_dense_m = Kratos.Matrix()
        KratosOA.ContainerVariableDataHolderUtils.Transpose(transpose_dense_m, dense_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_dense_m[j, i], dense_m[i, j])

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        transpose_sparse_m = Kratos.CompressedMatrix()
        KratosOA.ContainerVariableDataHolderUtils.Transpose(transpose_sparse_m, sparse_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()