//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez
//

#include "compute_distance_sensitivities_process.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
// Constructor for ComputeDistanceSensitivities Process
ComputeDistanceSensitivitiesProcess::ComputeDistanceSensitivitiesProcess(ModelPart& rModelPart, ModelPart& rSkinModelPart,
                    Parameters ThisParameters
                ):
    Process(),
    mrModelPart(rModelPart),
    mrSkinModelPart(rSkinModelPart)
{
}

void ComputeDistanceSensitivitiesProcess::InitializeNodalVariables(ModelPart::NodeType::Pointer pNode) {

    KRATOS_TRY;

    pNode->SetValue(SHAPE_SENSITIVITY, ZeroVector(3));

    KRATOS_CATCH("")
}

void ComputeDistanceSensitivitiesProcess::Execute()
{
    KRATOS_TRY;

    const int dim = 2;
    const int num_nodes = 3;

    std::size_t aux_nodal_id = 1;
    std::size_t temp_cond_id = 1;
    std::unordered_map<std::set<IndexType>, IndexType, KeyHasherRange<std::set<IndexType>>, KeyComparorRange<std::set<IndexType>> >  edge_nodes_to_int_node_id_map;

    for (auto& r_elem : mrModelPart.Elements()) {
        // auto elemental_distances = r_elem.GetValue(ELEMENTAL_DISTANCES);
        array_1d<double, 3> elemental_distances;
        auto& r_geometry = r_elem.GetGeometry();
        for (std::size_t i = 0; i < r_geometry.size(); i++) {
            elemental_distances[i] = r_geometry[i].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }

        // KRATOS_WATCH(elemental_distances)
        const bool is_split = PotentialFlowUtilities::CheckIfElementIsCutByDistance<dim, num_nodes>(elemental_distances);

        if (is_split) {
            Triangle2D3ModifiedShapeFunctions modified_shape_fun(r_elem.pGetGeometry(), Vector(elemental_distances));
            DivideGeometry::Pointer p_split_utility = modified_shape_fun.pGetSplittingUtil();
            p_split_utility->GenerateIntersectionsSkin();

            // Get the split geometries from the splitting pattern
            const unsigned int n_pos_interface_geom = (p_split_utility->mPositiveInterfaces).size();
            KRATOS_ERROR_IF(n_pos_interface_geom > 1) << "The number of interfaces is greater than 1 in elem #" << r_elem.Id() << std::endl;
            std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_interface_geometries;
            split_interface_geometries.reserve(n_pos_interface_geom);
            split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mPositiveInterfaces).begin(), (p_split_utility->mPositiveInterfaces).end());

            for (unsigned int i_int_geom = 0; i_int_geom < split_interface_geometries.size(); ++i_int_geom){
                DivideGeometry::IndexedPointGeometryPointerType p_int_sub_geom = split_interface_geometries[i_int_geom];
                GeometryData::KratosGeometryType p_int_sub_geom_type = p_int_sub_geom->GetGeometryType();
                const unsigned int sub_int_geom_n_nodes = p_int_sub_geom->PointsNumber();

                // Fill the new condition nodes array
                Condition::NodesArrayType sub_int_geom_nodes_array;
                for (unsigned int i_node = 0; i_node < sub_int_geom_n_nodes; ++i_node){

                    std::set<IndexType> edge_node_set;
                    DivideGeometry::IndexedPointType& sub_int_geom_node = p_int_sub_geom->operator[](i_node);
                    auto& int_pt_coords = sub_int_geom_node.Coordinates();
                    const unsigned int local_id = sub_int_geom_node.Id();
                    const unsigned int intersected_edge_id = local_id - num_nodes;

                    // Get the intersected edge node_i and node_j
                    const unsigned int node_i = (p_split_utility->GetEdgeIdsI())[intersected_edge_id];
                    const unsigned int node_j = (p_split_utility->GetEdgeIdsJ())[intersected_edge_id];
                    edge_node_set.insert(r_geometry[node_i].Id());
                    edge_node_set.insert(r_geometry[node_j].Id());

                    if (edge_nodes_to_int_node_id_map.find(edge_node_set) == edge_nodes_to_int_node_id_map.end()) {
                        auto new_p_node = mrSkinModelPart.CreateNewNode(aux_nodal_id, int_pt_coords[0], int_pt_coords[1], int_pt_coords[2]);
                        edge_nodes_to_int_node_id_map[edge_node_set] = aux_nodal_id++;
                        InitializeNodalVariables(new_p_node);
                        sub_int_geom_nodes_array.push_back(new_p_node);
                    } else {
                        auto existing_p_node = mrSkinModelPart.pGetNode(edge_nodes_to_int_node_id_map[edge_node_set]);
                        sub_int_geom_nodes_array.push_back(existing_p_node);
                    }
                }

                // Set the new condition geometry
                Geometry< Node<3> >::Pointer p_new_geom = SetNewConditionGeometry(
                    p_int_sub_geom_type,
                    sub_int_geom_nodes_array);

                // Set the new condition properties
                Properties::Pointer p_cond_prop = mrModelPart.pGetProperties(0);

                // Create the new condition
                Condition::Pointer p_new_cond = Kratos::make_intrusive<Condition>(temp_cond_id++, p_new_geom, p_cond_prop);
                mrSkinModelPart.AddCondition(p_new_cond);

                ComputeLocalDistanceSensitivity(r_geometry, p_new_cond->GetGeometry());
            }
        }
    }

    KRATOS_CATCH("");
}

ModelPart::GeometryType::Pointer ComputeDistanceSensitivitiesProcess::SetNewConditionGeometry(
    const GeometryData::KratosGeometryType &rOriginGeometryType,
    const Condition::NodesArrayType &rNewNodesArray)
{
    switch(rOriginGeometryType){
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            return Kratos::make_shared<Line2D2< Node<3> > >(rNewNodesArray);
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return Kratos::make_shared<Triangle3D3< Node<3> > >(rNewNodesArray);
        default:
            KRATOS_ERROR << "Implement the visualization for the intersection geometry type " << rOriginGeometryType;
    }
}

void ComputeDistanceSensitivitiesProcess::ComputeLocalDistanceSensitivity(ModelPart::GeometryType& rVolumeGeometry, ModelPart::GeometryType& rSkinGeometry) {

    KRATOS_TRY;

    const double epsilon = 1e-1;
    const std::size_t Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const std::size_t num_nodes = 3;//rVolumeGeometry.size();

    std::vector<array_1d<double, 3>> int_pts_vector;
    for (std::size_t i_dim = 0; i_dim < Dim; i_dim++) {
        int_pts_vector.push_back(rSkinGeometry[i_dim].Coordinates());
    }
    if (Dim == 2) {
        array_1d<double, 3> z_coord_2d = rSkinGeometry[0].Coordinates();
        z_coord_2d[2] = 1.0;
        int_pts_vector.push_back(z_coord_2d);
    }
    array_1d<double, num_nodes> initial_distances;
    Plane3D initial_intersection_plane = Plane3D(Point{int_pts_vector[0]}, Point{int_pts_vector[1]}, Point{int_pts_vector[2]});
    double min_distance = std::numeric_limits<double>::max();
    for (std::size_t i_vol_node = 0; i_vol_node < rVolumeGeometry.size(); i_vol_node++) {
        initial_distances[i_vol_node] = initial_intersection_plane.CalculateSignedDistance(rVolumeGeometry[i_vol_node]);
        min_distance = std::min(std::abs(initial_distances[i_vol_node]), min_distance);
    }

    const double delta = epsilon*min_distance;
    for (std::size_t i_int_node = 0; i_int_node < rSkinGeometry.size(); i_int_node++) {
        array_1d<double, 3> shape_sensitivity = ZeroVector(3);
        for (std::size_t i_dim = 0; i_dim < Dim; i_dim++) {
            array_1d<double, num_nodes> perturbed_distances;
            auto& this_perturbed_point = int_pts_vector[i_int_node];

            this_perturbed_point[i_dim] += delta;
            Plane3D perturbed_plane = Plane3D(Point{int_pts_vector[0]}, Point{int_pts_vector[1]}, Point{int_pts_vector[2]});
            for (std::size_t i_vol_node = 0; i_vol_node < rVolumeGeometry.size(); i_vol_node++) {
                perturbed_distances[i_vol_node] = perturbed_plane.CalculateSignedDistance(rVolumeGeometry[i_vol_node]);
            }
            this_perturbed_point[i_dim] -= delta;
            auto diff_vector = perturbed_distances - initial_distances;
            for (auto diff : diff_vector) {
                shape_sensitivity[i_dim] += diff / (delta*diff_vector.size());
            }
        }
        auto current_shape_sensitivity = rSkinGeometry[i_int_node].GetValue(SHAPE_SENSITIVITY);
        rSkinGeometry[i_int_node].SetValue(SHAPE_SENSITIVITY, current_shape_sensitivity+shape_sensitivity);
    }

    KRATOS_CATCH("");

}
}// Namespace Kratos
