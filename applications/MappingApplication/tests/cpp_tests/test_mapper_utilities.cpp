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
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/stream_serializer.h"
#include "utilities/cpp_tests_utilities.h"
#include "utilities/variable_utils.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_mappers/nearest_neighbor_mapper.h"

namespace Kratos {
namespace Testing {
namespace {

void CreateNodesForMapping(ModelPart& rModelPart, const int NumNodes)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    const int size = rModelPart.GetCommunicator().TotalProcesses();

    const int start_id = NumNodes * rank + 1;

    // creating nodes with random coordinates
    for (int i=0; i< NumNodes; ++i)
        rModelPart.CreateNewNode(i+start_id, i*0.1*rank*size+0.134,
                                             i*0.2+rank*3.48*size,
                                             i*0.3*rank*6.13*size);
}

} //empty namespace

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceEquationIds, KratosMappingApplicationSerialTestSuite)
{
    const int num_nodes = 11;

    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    CreateNodesForMapping(model_part, num_nodes);

    MapperUtilities::AssignInterfaceEquationIds(model_part.GetCommunicator());

    int idx = 0;

    for (const auto& r_node : model_part/*.GetCommunicator().LocalMesh()*/.Nodes())
    {
        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }
}



KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_CreateSearchLocalSystemsFromNodes, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_CHECK_GREATER_EQUAL(model_part.NumberOfNodes(), 0);

    std::vector<Kratos::unique_ptr<SearchLocalSystem>> search_local_systems;

    MapperUtilities::CreateSearchLocalSystemsFromNodes(
        NearestNeighborLocalSystem(nullptr),
        model_part.GetCommunicator(),
        search_local_systems);

    KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), search_local_systems.size());
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_EraseNodalVariable, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_CHECK_GREATER_EQUAL(model_part.NumberOfNodes(), 0);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.Has(DISPLACEMENT_X));
        r_node[DISPLACEMENT_X] = 15.3;
        KRATOS_CHECK(r_node.Has(DISPLACEMENT_X));
    }

    MapperUtilities::EraseNodalVariable(model_part, DISPLACEMENT_X);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.Has(DISPLACEMENT_X));
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_SaveCurrentConfiguration, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_CHECK_GREATER_EQUAL(model_part.NumberOfNodes(), 0);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.Has(CURRENT_COORDINATES));
    }

    MapperUtilities::SaveCurrentConfiguration(model_part);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK(r_node.Has(CURRENT_COORDINATES));
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.X(), r_node[CURRENT_COORDINATES][0]);
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Y(), r_node[CURRENT_COORDINATES][1]);
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Z(), r_node[CURRENT_COORDINATES][2]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_RestoreCurrentConfiguration, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_CHECK_GREATER_EQUAL(model_part.NumberOfNodes(), 0);

    KRATOS_CHECK_EXCEPTION_IS_THROWN(MapperUtilities::RestoreCurrentConfiguration(model_part), "Nodes do not have CURRENT_COORDINATES for restoring the current configuration!");

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.Has(CURRENT_COORDINATES));
        r_node.X() += 0.1;
        r_node.Y() -= 0.125;
        r_node.Z() += 0.33;
    }

    MapperUtilities::SaveCurrentConfiguration(model_part);

    // X = X0
    VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes());

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK(r_node.Has(CURRENT_COORDINATES));
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.X(), r_node.X0());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Y(), r_node.Y0());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Z(), r_node.Z0());
    }

    MapperUtilities::RestoreCurrentConfiguration(model_part);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.Has(CURRENT_COORDINATES));
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.X(), (r_node.X0()+0.1));
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Y(), (r_node.Y0()-0.125));
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Z(), (r_node.Z0()+0.33));
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_PointsAreCollinear, KratosMappingApplicationSerialTestSuite)
{
    Point p1(0,0,0);
    Point p2(1,0,0);
    Point p3(2,0,0);
    Point p4(2,1,0);

    KRATOS_CHECK(MapperUtilities::PointsAreCollinear(p1,p2,p3));
    KRATOS_CHECK(MapperUtilities::PointsAreCollinear(p2,p3,p1));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointsAreCollinear(p1,p2,p4));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointsAreCollinear(p1,p3,p4));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointsAreCollinear(p2,p3,p4));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointsAreCollinear(p2,p3,p4));
}

}  // namespace Kratos::Testing