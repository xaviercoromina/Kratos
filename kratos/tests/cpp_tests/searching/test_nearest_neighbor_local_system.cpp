//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "includes/serializer.h"
#include "testing/testing.h"
#include "includes/stream_serializer.h"
#include "includes/variables.h"
#include "searching/search_local_system/nearest_neighbor_local_system.h"
#include "searching/search_interface_info/nearest_neighbor_interface_info.h"

namespace Kratos::Testing {

using MatrixType = typename SearchLocalSystem::MatrixType;
using EquationIdVectorType = typename SearchLocalSystem::EquationIdVectorType;

using NodeType = Node<3>;

KRATOS_TEST_CASE_IN_SUITE(SearchLocalSystem_BasicTests, KratosCoreFastSuite)
{
    // This test covers the basic functionalities provided by the "SearchLocalSystem"
    // A "NearestNeighborLocalSystem" is being used since "SearchLocalSystem" is a pure virtual class

    Point coords_1(1.0, 2.45, -23.8);
    auto node_local_sys(Kratos::make_intrusive<NodeType>(5, coords_1));

    NearestNeighborLocalSystem local_sys(node_local_sys.get());

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(local_sys.Coordinates()[i], coords_1[i]);

    KRATOS_CHECK_IS_FALSE(local_sys.HasInterfaceInfo());
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystem_BasicTests, KratosCoreFastSuite)
{
    auto node_local_sys(Kratos::make_intrusive<NodeType>(8, 1.0, 2.5, -5.0));

    NearestNeighborLocalSystem local_sys(node_local_sys.get());

    // Computing the local system
    // this should return nothing since no InterfaceInfos are available
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys.EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids.size(), 0);

    local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 0);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), 0);
    KRATOS_CHECK_EQUAL(origin_ids2.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 0);

    std::stringstream str_steam;
    local_sys.PairingInfo(str_steam, 4);
    KRATOS_CHECK_STRING_EQUAL(str_steam.str(), "NearestNeighborLocalSystem based on Node #8 at Coordinates 1 | 2.5 | -5");
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystem_ComputeLocalSystem, KratosCoreFastSuite)
{
    const int dest_id = 13;

    auto node_local_sys(Kratos::make_intrusive<NodeType>(5, 1.0, 2.5, -5.0));
    node_local_sys->SetValue(INTERFACE_EQUATION_ID, dest_id);

    NearestNeighborLocalSystem local_sys(node_local_sys.get());

    // Create the NearestNeighborInfos to be used by the NearestNeighborLocalSystem
    auto node_1(Kratos::make_intrusive<NodeType>(1, 18.0, 2.7, 30.0));
    auto node_2(Kratos::make_intrusive<NodeType>(3, 1.0, 2.5, -3.0)); // this is the nearest neighbor

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));

    const int expected_id_found = 67;

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    SearchInterfaceInfo::Pointer p_nearest_neighbor_info_1(Kratos::make_shared<NearestNeighborInterfaceInfo>(local_sys.Coordinates(), 0, 0));
    SearchInterfaceInfo::Pointer p_nearest_neighbor_info_2(Kratos::make_shared<NearestNeighborInterfaceInfo>(local_sys.Coordinates(), 0, 0));

    p_nearest_neighbor_info_1->ProcessSearchResult(*interface_node_1);
    p_nearest_neighbor_info_2->ProcessSearchResult(*interface_node_2);

    local_sys.AddInterfaceInfo(p_nearest_neighbor_info_1);
    local_sys.AddInterfaceInfo(p_nearest_neighbor_info_2);

    // Computing the local system
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys.EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), 1);
    KRATOS_CHECK_EQUAL(destination_ids.size(), 1);

    local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 1);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), 1);
    KRATOS_CHECK_EQUAL(origin_ids2.size(), 1);
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 1);

    KRATOS_CHECK_DOUBLE_EQUAL(local_mapping_matrix(0,0), 1.0);
    KRATOS_CHECK_EQUAL(origin_ids2[0], expected_id_found);
    KRATOS_CHECK_EQUAL(destination_ids2[0], dest_id);
}

}  // namespace Kratos::Testing