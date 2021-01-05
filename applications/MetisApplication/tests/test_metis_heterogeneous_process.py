import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import pathlib
import itertools


def GetFilePath(file_name):
    return str(pathlib.Path( __file__ ).parent / pathlib.Path(file_name))


class KratosMetisHeterogeneousProcessTests(KratosUnittest.TestCase):

    def setUp(self):
        self.file_name = "test_metis_heterogeneous_process"
        self.import_export_flags = KratosMultiphysics.IO.READ | KratosMultiphysics.IO.WRITE | KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR
        self.import_flags = KratosMultiphysics.IO.READ | KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR
        self.export_flags = KratosMultiphysics.IO.WRITE | KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR
        self.communicator = KratosMultiphysics.DataCommunicator.GetDefault()

        if (self.communicator.Size() == 1):
            self.skipTest("Metis partitioning requires MPI parallelism")

    def tearDown(self):
        self.DeleteMDPAs()

    @property
    def is_main_rank(self):
        return self.communicator.Rank() == 0

    @property
    def partitioned_file_name(self):
        return self.GetPartitionedFileName(self.communicator.Rank())

    def DeleteMDPAs(self):
        """Delete all model part files"""
        if self.is_main_rank:
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(GetFilePath(self.file_name + ".mdpa"))
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(GetFilePath(self.file_name + ".time"))
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(GetFilePath(self.partitioned_file_name + ".mdpa"))
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(GetFilePath(self.partitioned_file_name + ".time"))

    def GetPartitionedFileName(self, rank):
        return self.file_name + "_rank_{rank}".format(rank=rank)

    def CheckModel(self):
        """
        Compare the current partitioned model (defined by self.partitioned_file_name)
        to the main model (defined by self.file_name).

        The following checks are performed on the partitioned model parts:
            - contains at least one vertex and element
            - stored elements have all their nodes stored as well
            - local and global nodes match
            - local and global elements match

        The following global checks are performed:
            - all nodes and elements are present
        """
        if self.is_main_rank:
            model = KratosMultiphysics.Model()
            model_parts = []

            # Read main model part
            main_model_part = model.CreateModelPart("Main")
            KratosMultiphysics.ModelPartIO(
                GetFilePath(self.file_name), self.import_flags).ReadModelPart(main_model_part)

            # Create dictionaries that count instances of nodes and elements in partitioned model parts
            # -> these will be checked at the end of this function
            node_instance_counter = dict([(node.Id, 0) for node in main_model_part.Nodes])
            element_instance_counter = dict([(element.Id, 0) for element in main_model_part.Elements])

            # Define helper functions for checks
            def IsNodeInModelPart(node, model_part):
                """Check for matching node Ids."""
                for local_node in model_part.Nodes:
                    if node.Id == local_node.Id:
                        return True
                return False

            def GetNodeInModelPart(node, model_part):
                """Return node from the model part with a matching Id."""
                for local_node in model_part.Nodes:
                    if node.Id == local_node.Id:
                        return local_node
                self.fail("Cannot find node in model part")

            def NodesWithIdenticalCoordinates(node0, node1, tolerance=0.0):
                """Check whether the nodes are sufficiently close to each other. (1-norm)"""
                norm = abs(node0.X0 - node1.X0) + abs(node0.Y0 - node1.Y0) + abs(node0.Z0 - node1.Z0)
                return norm <= tolerance


            def GetElementInModelPart(element, model_part):
                """Return element from the model part with a matching Id."""
                for local_element in model_part.Elements:
                    if element.Id == local_element.Id:
                        return local_element
                self.fail("Cannot find element in model part")

            def ElementsWithIdenticalNodeIds(element0,element1):
                """Check whether the elements have identical node Ids in the same order."""
                for node0, node1 in zip(element0.GetNodes(), element1.GetNodes()):
                    if node0.Id != node1.Id:
                        return False
                return True

            # Loop through partitioned model parts
            for rank in range(self.communicator.Size()):

                # Read partitioned model part
                model_parts.append(model.CreateModelPart("Rank{rank}".format(rank=rank)))
                current_model_part = model_parts[-1]

                KratosMultiphysics.ModelPartIO(
                    GetFilePath(self.GetPartitionedFileName(rank)), self.import_flags).ReadModelPart(current_model_part)

                # Perform checks on partitioned model part
                self.assertGreater(current_model_part.NumberOfNodes(), 0)
                self.assertGreater(current_model_part.NumberOfElements(),0)

                for node in current_model_part.Nodes:
                    # Check whether the node is in the main model part
                    self.assertTrue(IsNodeInModelPart(node, main_model_part), msg=node)

                    # Check whether the node has the same coordinates as its counterpart
                    node_in_main = GetNodeInModelPart(node, main_model_part)
                    self.assertTrue(NodesWithIdenticalCoordinates(node, node_in_main), msg=node)

                    # Update instance counter
                    node_instance_counter[node.Id] += 1

                for element in current_model_part.Elements:
                    # Check whether the element is in the main model part
                    element_in_main = GetElementInModelPart(element, main_model_part)
                    self.assertTrue(ElementsWithIdenticalNodeIds(element, element_in_main), msg=element)

                    # Check whether all nodes of the element are present in the model part
                    for node in element.GetNodes():
                        self.assertTrue(IsNodeInModelPart(node, current_model_part), msg=node)

                    # Update instance counter
                    element_instance_counter[element.Id] += 1

            # Check whether the partitions contain all nodes
            for iNode, number_of_instances in node_instance_counter.items():
                self.assertGreater(number_of_instances, 0, msg=iNode)

            # Check whether the partitions contain all elements
            for iElement, number_of_instances in element_instance_counter.items():
                self.assertGreater(number_of_instances, 0, msg=iElement)

        self.communicator.Barrier()

    def MakeLineModel(self):
        """Construct a model part consisting of line elements and write it to an mdpa file."""
        if self.is_main_rank:
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("Main")

            # Generate nodes on a line
            number_of_nodes = 25
            for iNode in range(number_of_nodes):
                model_part.CreateNewNode(iNode + 1, iNode, 0.0, 0.0)

            model_properties = model_part.CreateNewProperties(0)

            # Generate line mesh
            #element_source_node_ids = itertools.chain(
            #    range(int(number_of_nodes / 2) - 1),
            #    range(int(number_of_nodes / 2) + 1, number_of_nodes - 1)
            #)  # <-- this creates a hanging node in the middle

            element_source_node_ids = range(number_of_nodes - 1)

            for iElement, iSourceNode in enumerate(element_source_node_ids):
                model_part.CreateNewElement("Element2D2N", iElement + 1, (iSourceNode + 1, iSourceNode + 2), model_properties)
                """
                NOTE: elements can be created with 0-based Ids,
                but will be shifted to 1-based Ids in the partitioned model.
                In that case, all element-matching tests would fail.
                Is this behaviour intentional?
                """

            # Write model part to file
            model_part_output = KratosMultiphysics.ModelPartIO(
                GetFilePath(self.file_name), self.export_flags)
            model_part_output.WriteModelPart(model_part)
        self.communicator.Barrier()

    def test_LineModel(self):
        """Check metis partitioning on a line mesh."""
        self.MakeLineModel()

        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Main")

        # Get arguments for metis
        model_part_serial_io = KratosMultiphysics.ModelPartIO(
            GetFilePath(self.file_name), self.import_flags)
        model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(
            GetFilePath(self.file_name), self.import_flags)

        number_of_partitions = self.communicator.Size()
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        verbosity = 0
        synchronize_conditions = True

        # Partition in memory
        partitioner = KratosMetis.MetisDivideHeterogeneousInputInMemoryProcess(
            model_part_io, model_part_serial_io, number_of_partitions, domain_size, verbosity, synchronize_conditions)
        partitioner.Execute()

        model_part_serial_io.ReadModelPart(model_part)

        # Write partitioned model part to a file
        KratosMultiphysics.ModelPartIO(
            GetFilePath(self.partitioned_file_name), self.export_flags).WriteModelPart(model_part)
        self.communicator.Barrier()

        self.CheckModel()
        self.DeleteMDPAs()



if __name__ == "__main__":
    KratosUnittest.main()