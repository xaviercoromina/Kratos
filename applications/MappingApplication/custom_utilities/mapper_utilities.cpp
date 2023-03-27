//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/stream_serializer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {
namespace MapperUtilities {

typedef std::size_t SizeType;
typedef std::size_t IndexType;

void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator)
{
    if (rModelPartCommunicator.GetDataCommunicator().IsNullOnThisRank()) {
        return;
    }

    const int num_nodes_local = rModelPartCommunicator.LocalMesh().NumberOfNodes();
    int num_nodes_accumulated = rModelPartCommunicator.GetDataCommunicator().ScanSum(num_nodes_local);
    const int start_equation_id = num_nodes_accumulated - num_nodes_local;
    const auto nodes_begin = rModelPartCommunicator.LocalMesh().NodesBegin();

    IndexPartition<unsigned int>(num_nodes_local).for_each(
        [nodes_begin, start_equation_id](unsigned int i){
            (nodes_begin + i)->SetValue(INTERFACE_EQUATION_ID, start_equation_id + i);
        }
    );

    rModelPartCommunicator.SynchronizeNonHistoricalVariable(INTERFACE_EQUATION_ID);
}

void CheckInterfaceModelParts(const int CommRank)
{
    // const int num_nodes_origin = MapperUtilities::ComputeNumberOfNodes(mrModelPartOrigin);
    // const int num_conditions_origin = MapperUtilities::ComputeNumberOfConditions(mrModelPartOrigin);
    // const int num_elements_origin = MapperUtilities::ComputeNumberOfElements(mrModelPartOrigin);

    // const int num_nodes_destination = MapperUtilities::ComputeNumberOfNodes(mrModelPartDestination);
    // const int num_conditions_destination = MapperUtilities::ComputeNumberOfConditions(mrModelPartDestination);
    // const int num_elements_destination = MapperUtilities::ComputeNumberOfElements(mrModelPartDestination);

    // // Check if the ModelPart contains entities
    // KRATOS_ERROR_IF(num_nodes_origin + num_conditions_origin + num_elements_origin < 1)
    //     << "Neither Nodes nor Conditions nor Elements found "
    //     << "in the Origin ModelPart" << std::endl;

    // KRATOS_ERROR_IF(num_nodes_destination + num_conditions_destination + num_elements_destination < 1)
    //     << "Neither Nodes nor Conditions nor Elements found "
    //     << "in the Destination ModelPart" << std::endl;

    // // Check if the inpt ModelParts contain both Elements and Conditions
    // // This is NOT possible, bcs the InterfaceObjects are constructed
    // // with whatever exists in the Modelpart (see the InterfaceObjectManagerBase,
    // // function "InitializeInterfaceGeometryObjectManager")
    // KRATOS_ERROR_IF(num_conditions_origin > 0 && num_elements_origin > 0)
    //     << "Origin ModelPart contains both Conditions and Elements "
    //     << "which is not permitted" << std::endl;

    // KRATOS_ERROR_IF(num_conditions_destination > 0 && num_elements_destination > 0)
    //     << "Destination ModelPart contains both Conditions and Elements "
    //     << "which is not permitted" << std::endl;

    // if (mEchoLevel >= 2) {
    //     std::vector<double> model_part_origin_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartOrigin);
    //     std::vector<double> model_part_destination_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartDestination);

    //     bool bbox_overlapping = MapperUtilities::ComputeBoundingBoxIntersection(
    //                                                 model_part_origin_bbox,
    //                                                 model_part_destination_bbox);
    //     if(CommRank == 0)
    //     {
    //         if (!bbox_overlapping) {
    //             std::cout << "MAPPER WARNING, the bounding boxes of the "
    //                         << "Modelparts do not overlap! "
    //                         << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
    //                                                                         model_part_destination_bbox)
    //                         << std::endl;
    //         } else if (mEchoLevel >= 3)
    //         {
    //             std::cout << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
    //                                                                         model_part_destination_bbox)
    //                         << std::endl;
    //         }
    //     }
    // }
}

void CreateSearchLocalSystemsFromNodes(const SearchLocalSystem& rSearchLocalSystemPrototype,
                                       const Communicator& rModelPartCommunicator,
                                       std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems)
{
    const std::size_t num_nodes = rModelPartCommunicator.LocalMesh().NumberOfNodes();
    const auto nodes_ptr_begin = rModelPartCommunicator.LocalMesh().Nodes().ptr_begin();

    if (rLocalSystems.size() != num_nodes) {
        rLocalSystems.resize(num_nodes);
    }

    IndexPartition<std::size_t>(num_nodes).for_each([&](std::size_t i){
        InterfaceObject::NodePointerType p_node = (nodes_ptr_begin + i)->get();
        rLocalSystems[i] = rSearchLocalSystemPrototype.Create(p_node);
    });

    if (rModelPartCommunicator.GetDataCommunicator().IsDefinedOnThisRank()) {
        int num_local_systems = rModelPartCommunicator.GetDataCommunicator().SumAll((int)(rLocalSystems.size())); // int bcs of MPI

        KRATOS_ERROR_IF_NOT(num_local_systems > 0) << "No mapper local systems were created" << std::endl;
    }
}

void CreateSearchLocalSystemsFromGeometries(const SearchLocalSystem& rSearchLocalSystemPrototype,
                                            const Communicator& rModelPartCommunicator,
                                            std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems)
{
    const std::size_t num_conditions = rModelPartCommunicator.LocalMesh().NumberOfConditions();
    const auto cond_begin = rModelPartCommunicator.LocalMesh().ConditionsBegin();

    if (rLocalSystems.size() != num_conditions) rLocalSystems.resize(num_conditions);

    IndexPartition<std::size_t>(num_conditions).for_each([&](std::size_t i){
        InterfaceObject::GeometryPointerType p_geom = &((cond_begin+i)->GetGeometry());
        rLocalSystems[i] = rSearchLocalSystemPrototype.Create(p_geom);
    });

    if (rModelPartCommunicator.GetDataCommunicator().IsDefinedOnThisRank()) {
        const int num_local_systems = rModelPartCommunicator.GetDataCommunicator().SumAll((int)(rLocalSystems.size())); // int bcs of MPI

        KRATOS_ERROR_IF_NOT(num_local_systems > 0) << "No mapper local systems were created" << std::endl;
    }
}

void SaveCurrentConfiguration(ModelPart& rModelPart)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.SetValue(CURRENT_COORDINATES, rNode.Coordinates());
    });

    KRATOS_CATCH("");
}

void RestoreCurrentConfiguration(ModelPart& rModelPart)
{
    KRATOS_TRY;

    if (rModelPart.NumberOfNodes() > 0) {
        KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->Has(CURRENT_COORDINATES)) << "Nodes do not have CURRENT_COORDINATES for restoring the current configuration!" << std::endl;

        block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
            noalias(rNode.Coordinates()) = rNode.GetValue(CURRENT_COORDINATES);
            rNode.GetData().Erase(CURRENT_COORDINATES);
        });
    }

    KRATOS_CATCH("");
}

} // namespace MapperUtilities
} // namespace Kratos.
