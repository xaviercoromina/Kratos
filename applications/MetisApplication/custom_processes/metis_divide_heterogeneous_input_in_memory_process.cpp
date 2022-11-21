//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//


// System includes


// External includes


// Project includes
#include "mpi/includes/mpi_data_communicator.h"
#include "metis_divide_heterogeneous_input_in_memory_process.h"


namespace Kratos {

void MetisDivideHeterogeneousInputInMemoryProcess::Execute()
{
    MPI_Comm this_comm = MPIDataCommunicator::GetMPICommunicator(mrDataComm);

    int mpi_rank = mrDataComm.Rank();
    int mpi_size = mrDataComm.Size();

    int * msgSendSize = new int[mpi_size];
    int * msgRecvSize = new int[mpi_size];

    const char ** mpi_send_buffer = new const char * [mpi_size];
    char ** mpi_recv_buffer = new char * [mpi_size];
    std::string * str = new std::string[mpi_size];

    // Set size
    for(int i = 0; i < mpi_size; i++) {
        msgSendSize[i] = 0;
        msgRecvSize[i] = 0;
    }

    // Transfer Streams
    Kratos::shared_ptr<std::iostream> * streams = new Kratos::shared_ptr<std::iostream>[mpi_size];
    std::stringbuf * stringbufs = new std::stringbuf[mpi_size];

    for(auto i = 0; i < mpi_size; i++) {
        streams[i] = Kratos::shared_ptr<std::iostream>(new std::iostream(&stringbufs[i]));
    }

    // Main process calculates the partitions and writes the result into temporal streams
    if(mpi_rank == 0) {
        // Read nodal graph from input

        IO::ConnectivitiesContainerType KratosFormatNodeConnectivities;

        SizeType NumNodes = BaseType::mrIO.ReadNodalGraph(KratosFormatNodeConnectivities);

        // Write connectivity data in CSR format
        idxtype* NodeIndices = 0;
        idxtype* NodeConnectivities = 0;

        ConvertKratosToCSRFormat(KratosFormatNodeConnectivities, &NodeIndices, &NodeConnectivities);

        std::vector<idxtype> NodePartition;
        PartitionNodes(NumNodes,NodeIndices,NodeConnectivities,NodePartition);

        // Free some memory we no longer need
        delete [] NodeIndices;
        delete [] NodeConnectivities;

        // Partition elements
        IO::ConnectivitiesContainerType ElementConnectivities;
        SizeType NumElements =  BaseType::mrIO.ReadElementsConnectivities(ElementConnectivities);
        if (NumElements != ElementConnectivities.size())
        {
            std::stringstream Msg;
            Msg << std::endl;
            Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
            Msg << "Read " << NumElements << " elements, but element list has " << ElementConnectivities.size() << " entries." << std::endl;
            Msg << "Elements are most likely not correlatively numbered." << std::endl;

            KRATOS_ERROR << Msg.str();
        }

        std::vector<idxtype> ElementPartition;

        if (mSynchronizeConditions)
            PartitionElementsSynchronous(NodePartition,ElementConnectivities,ElementPartition);
        else
            PartitionMesh(NodePartition,ElementConnectivities,ElementPartition);

        // Partition conditions
        IO::ConnectivitiesContainerType ConditionConnectivities;
        SizeType NumConditions = BaseType::mrIO.ReadConditionsConnectivities(ConditionConnectivities);
        if (NumConditions != ConditionConnectivities.size())
        {
            std::stringstream Msg;
            Msg << std::endl;
            Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
            Msg << "Read " << NumConditions << " conditions, but condition list has " << ConditionConnectivities.size() << " entries." << std::endl;
            Msg << "Conditions are most likely not correlatively numbered." << std::endl;

            KRATOS_ERROR << Msg.str();
        }

        std::vector<idxtype> ConditionPartition;

        if (mSynchronizeConditions)
            PartitionConditionsSynchronous(NodePartition,ElementPartition,ConditionConnectivities,ElementConnectivities,ConditionPartition);
        else
            PartitionMesh(NodePartition,ConditionConnectivities,ConditionPartition);

        // Detect hanging nodes (nodes that belong to a partition where no local elements have them) and send them to another partition.
        // Hanging nodes should be avoided, as they can cause problems when setting the Dofs
        RedistributeHangingNodes(NodePartition,ElementPartition,ElementConnectivities,ConditionPartition,ConditionConnectivities);

        // Coloring
        GraphType DomainGraph = zero_matrix<int>(mNumberOfPartitions);
        CalculateDomainsGraph(DomainGraph,NumElements,ElementConnectivities,NodePartition,ElementPartition);
        CalculateDomainsGraph(DomainGraph,NumConditions,ConditionConnectivities,NodePartition,ConditionPartition);

        int NumColors;
        GraphType ColoredDomainGraph;
        GraphColoringProcess(mNumberOfPartitions,DomainGraph,ColoredDomainGraph,NumColors).Execute();

        if (mVerbosity > 0) {
            KRATOS_INFO("NumColors") << NumColors << std::endl;
        }

        if (mVerbosity > 2) {
            KRATOS_INFO("ColoredDomainGraph") << ColoredDomainGraph << std::endl;
        }

        // Write partition info into separate input files
        IO::PartitionIndicesContainerType nodes_all_partitions;
        IO::PartitionIndicesContainerType elements_all_partitions;
        IO::PartitionIndicesContainerType conditions_all_partitions;

        // Create lists containing all nodes/elements/conditions known to each partition
        DividingNodes(nodes_all_partitions, ElementConnectivities, ConditionConnectivities, NodePartition, ElementPartition, ConditionPartition);
        DividingElements(elements_all_partitions, ElementPartition);
        DividingConditions(conditions_all_partitions, ConditionPartition);

        if (mVerbosity > 1) {
            std::cout << "Final list of nodes known by each partition" << std::endl;
            for(SizeType i = 0 ; i < NumNodes ; i++) {
                std::cout << "Node #" << i+1 << "->";
                for(std::vector<std::size_t>::iterator j = nodes_all_partitions[i].begin() ; j != nodes_all_partitions[i].end() ; j++) {
                    std::cout << *j << ",";
                }
                std::cout << std::endl;
            }
        }

        IO::PartitionIndicesType io_nodes_partitions(NodePartition.begin(), NodePartition.end());
        IO::PartitionIndicesType io_elements_partitions(ElementPartition.begin(), ElementPartition.end());
        IO::PartitionIndicesType io_conditions_partitions(ConditionPartition.begin(), ConditionPartition.end());

        // Write files
        mrIO.DivideInputToPartitions(
            streams, mNumberOfPartitions, ColoredDomainGraph,
            io_nodes_partitions, io_elements_partitions, io_conditions_partitions,
            nodes_all_partitions, elements_all_partitions, conditions_all_partitions
        );
    }

    // Calculate the message and prepare the buffers
    if(mpi_rank == 0) {
        for(auto i = 0; i < mpi_size; i++) {
            str[i] = stringbufs[i].str();
            msgSendSize[i] = str[i].size();
            mpi_send_buffer[i] = str[i].c_str();
        }
    }

    // Send the message size to all processes
    MPI_Scatter(msgSendSize,1,MPI_INT,&msgRecvSize[mpi_rank],1,MPI_INT,0,this_comm);

    // Calculate the number of events:
    auto NumberOfCommunicationEvents = 1 + mpi_size * !mpi_rank;
    auto NumberOfCommunicationEventsIndex = 0;

    // Prepare the communication events
    MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
    MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

    // Set up all receive and send events
    if( mpi_rank == 0) {
        for(auto i = 0; i < mpi_size; i++) {
            char* aux_char = const_cast<char*>(mpi_send_buffer[i]);
            MPI_Isend(aux_char,msgSendSize[i],MPI_CHAR,i,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);
        }
    }

    // Recieve the buffers
    mpi_recv_buffer[mpi_rank] = (char *)malloc(sizeof(char) * msgRecvSize[mpi_rank]);
    MPI_Irecv(mpi_recv_buffer[mpi_rank],msgRecvSize[mpi_rank],MPI_CHAR,0,0,this_comm,&reqs[NumberOfCommunicationEventsIndex++]);

    // Wait untill all communications finish
    if( MPI_Waitall(NumberOfCommunicationEvents, reqs, stats) != MPI_SUCCESS ) {
        KRATOS_ERROR << "Error in metis_partition_mem" << std::endl;
    }

    mrDataComm.Barrier();

    if(mpi_rank != 0) {
        streams[mpi_rank]->write(mpi_recv_buffer[mpi_rank], msgRecvSize[mpi_rank]);
    }

    if (mVerbosity > 1) {
        std::ofstream debug_ofstream("MetisDivideHeterogeneousInputInMemoryProcess_debug_modelpart_"+std::to_string(mpi_rank)+".mdpa");
        debug_ofstream << stringbufs[mpi_rank].str() << std::endl;
    }

    // TODO: Try to come up with a better way to change the buffer.
    mrSerialIO.SwapStreamSource(streams[mpi_rank]);

    // Free buffers
    free(mpi_recv_buffer[mpi_rank]);

    delete [] reqs;
    delete [] stats;

    delete [] mpi_recv_buffer;
    delete [] mpi_send_buffer;
    delete [] str;

    delete [] msgSendSize;
    delete [] msgRecvSize;
}

} // namespace Kratos