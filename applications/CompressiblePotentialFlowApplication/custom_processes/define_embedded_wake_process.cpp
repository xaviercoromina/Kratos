//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez
//


#include "define_embedded_wake_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"


namespace Kratos
{
// Constructor for DefineEmbeddedWakeProcess Process
DefineEmbeddedWakeProcess::DefineEmbeddedWakeProcess(ModelPart& rModelPart,
                    ModelPart& rWakeModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrWakeModelPart(rWakeModelPart)
{}

void DefineEmbeddedWakeProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]>2) << "DOMAIN_SIZE is greater than 2. DefineEmbeddedWakeProcess is only implemented for 2D cases!" << std::endl;

    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();
    MarkKuttaWakeElements();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        it_elem->Set(TO_SPLIT, false);
        if (it_elem->Is(BOUNDARY)) {
            it_elem->Set(TO_SPLIT);
        }
        auto geometry_elemental_distances = it_elem->GetValue(GEOMETRY_ELEMENTAL_DISTANCES);
        it_elem->SetValue(ELEMENTAL_DISTANCES, geometry_elemental_distances);
    }

    // auto& wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];
    // auto& wake_angle = mrModelPart.GetProcessInfo()[ROTATION_ANGLE];
    // double max_inactive_x = wake_origin[0];
    // for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
    //     auto it_elem = mrModelPart.ElementsBegin() + i;
    //     for(unsigned int i_node = 0; i_node<3; i_node++){
    //         geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    //         distances(i_node) = geometry_distances[i_node];
    //     }
    //     const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);
    // }

    KRATOS_CATCH("");
}


void DefineEmbeddedWakeProcess::ComputeDistanceToWake(){

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        auto geometry_elemental_distances = it_elem->GetValue(ELEMENTAL_DISTANCES);
        it_elem->SetValue(GEOMETRY_ELEMENTAL_DISTANCES, geometry_elemental_distances);
        it_elem->SetValue(WAKE, false);
        it_elem->SetValue(KUTTA, false);
        if (it_elem->Is(TO_SPLIT)) {
            it_elem->Set(BOUNDARY);
        }
    }

    CalculateDiscontinuousDistanceToSkinProcess<2> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->SetValue(WAKE_DISTANCE, 0.0);
        it_node->SetValue(KUTTA, 0);
        it_node->SetValue(WAKE, 0);
    }
}

void DefineEmbeddedWakeProcess::MarkWakeElements(){

    ModelPart& deactivated_model_part = mrModelPart.CreateSubModelPart("deactivated_model_part");
    std::vector<std::size_t> deactivated_elements_id_list;
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;

        BoundedVector<double, 3> nodal_distances_to_wake = it_elem->GetValue(ELEMENTAL_DISTANCES);
        it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(nodal_distances_to_wake);

        BoundedVector<double,3> geometry_distances;
        Vector distances(3);
        for(unsigned int i_node = 0; i_node<3; i_node++){
            geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
            distances(i_node) = geometry_distances[i_node];
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);
        if (is_embedded){
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_distances));
            // Computing Normal
            std::vector<Vector> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
            double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
            auto unit_normal = cut_normal[0]/norm_normal;
            it_elem->SetValue(VELOCITY_LOWER,unit_normal);
        }
        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element && it_elem->Is(ACTIVE)) {
            auto& r_geometry = it_elem->GetGeometry();
            if (is_embedded){
                #pragma omp critical
                {
                    deactivated_elements_id_list.push_back(it_elem->Id());
                }
                // it_elem->Set(ACTIVE, false);
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(KUTTA, 1);
                    r_geometry[i].UnSetLock();
                }
            }
            else{
                it_elem->SetValue(WAKE, true);
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE, 1);
                    r_geometry[i].UnSetLock();
                }
            }
            for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                r_geometry[i].UnSetLock();
            }
        } else {
            it_elem->SetValue(WAKE, false);
        }
    }
    deactivated_model_part.AddElements(deactivated_elements_id_list);
}

ModifiedShapeFunctions::Pointer DefineEmbeddedWakeProcess::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){
    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");

    std::vector<std::size_t> trailing_edge_node_list;
    std::vector<std::size_t> structure_element_list;
    std::vector<std::size_t> set_kutta_element_list;


    if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
        mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");


    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
        auto& r_geometry = it_elem->GetGeometry();
        std::size_t wake_nodes_counter = 0;
        std::size_t lower_nodes_counter = 0;
        for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
            if(r_geometry[i_node].GetValue(WAKE)) {
                wake_nodes_counter++;
            }
            if(r_geometry[i_node].GetValue(LOWER_SURFACE)) {
                lower_nodes_counter++;

            }
        }

        if(wake_nodes_counter == r_geometry.size()-1){
            // /*########## BELOW NEWGATIVE APPROACH  ##########*/

            // it_elem->Set(STRUCTURE, true);
            // it_elem->SetValue(WAKE, true);
            // BoundedVector<double, 3> distances_to_wake = it_elem->GetValue(ELEMENTAL_DISTANCES);

            // for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
            //     // if (r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE)<0.0) {
            //     if (distances_to_wake[i_node] > 0.0 && r_geometry[i_node].GetValue(KUTTA)) {
            //         r_geometry[i_node].SetValue(TRAILING_EDGE, true);
            //         r_geometry[i_node].SetValue(WING_TIP, true);
            //         trailing_edge_node_list.push_back(r_geometry[i_node].Id());

            //     }

            //     // if (distances_to_wake[i_node] < 0.0 && r_geometry[i_node].GetValue(KUTTA)) {
            //     //     // r_geometry[i_node].SetValue(TRAILING_EDGE, true);
            //     //     r_geometry[i_node].SetValue(WING_TIP, true);
            //     //     // trailing_edge_node_list.push_back(r_geometry[i_node].Id());

            //     // }
            // }
            // double max_x = -1e30;
            // std::size_t max_node_id = -1;
            // for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
            //     if (r_geometry[i_node].X() > max_x ) {
            //         max_x = r_geometry[i_node].X();
            //         max_node_id = r_geometry[i_node].Id();

            //     }
            // }

            // for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
            //     if (r_geometry[i_node].Id() == max_node_id) {
            //         r_geometry[i_node].SetValue(WING_TIP, false);
            //     }
            // }
            // /*###############################################*/


            /*########## BELOW POSITIVE APPROACH  ##########*/
            it_elem->SetValue(KUTTA, true);
            for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
                auto r_node_elem_candidates = r_geometry[i_node].GetValue(NEIGHBOUR_ELEMENTS);
                for (auto r_elem : r_node_elem_candidates) {
                    if (r_elem.GetValue(WAKE)) {
                        std::size_t kutta_nodes_counter = 0;
                        auto& r_neighbour_geometry = r_elem.GetGeometry();
                        for (unsigned int i_node_neigh= 0; i_node_neigh < r_neighbour_geometry.size(); i_node_neigh++) {

                            if(r_neighbour_geometry[i_node_neigh].GetValue(WAKE_DISTANCE)>0.0) {
                                kutta_nodes_counter++;
                            }
                        }

                        if (kutta_nodes_counter <2) {
                            set_kutta_element_list.push_back(r_elem.Id());
                        }
                        for (unsigned int i_node_neigh= 0; i_node_neigh < r_neighbour_geometry.size(); i_node_neigh++) {
                            // if(r_neighbour_geometry[i_node_neigh].GetValue(KUTTA) && r_neighbour_geometry[i_node_neigh].GetValue(WAKE_DISTANCE)>0.0) {
                            if(r_neighbour_geometry[i_node_neigh].GetValue(KUTTA) && r_neighbour_geometry[i_node_neigh].GetValue(WAKE_DISTANCE)>0.0 && kutta_nodes_counter>1) {
                                    structure_element_list.push_back(r_elem.Id());

                                    r_neighbour_geometry[i_node_neigh].SetValue(TRAILING_EDGE, true);
                                    r_neighbour_geometry[i_node_neigh].SetValue(WING_TIP, true);
                                    r_neighbour_geometry[i_node_neigh].SetValue(AIRFOIL, true);
                                    trailing_edge_node_list.push_back(r_neighbour_geometry[i_node_neigh].Id());
                            } else if ( r_neighbour_geometry[i_node_neigh].GetValue(WAKE_DISTANCE)<0.0 && kutta_nodes_counter>1){
                                // KRATOS_WATCH(r_elem.Id())
                                // KRATOS_WATCH(r_neighbour_geometry[i_node_neigh].Id())
                                // KRATOS_WATCH(kutta_nodes_counter)
                                // KRATOS_WATCH(r_neighbour_geometry[i_node_neigh].GetValue(WAKE_DISTANCE))
                                // KRATOS_WATCH(r_neighbour_geometry[i_node_neigh].GetValue(KUTTA))
                                // r_neighbour_geometry[i_node_neigh].SetValue(WING_TIP, true);
                            }
                        }
                    }
                }
            }
        }
    }
    KRATOS_WATCH(trailing_edge_node_list)
    KRATOS_WATCH(structure_element_list)
    KRATOS_ERROR_IF(structure_element_list.size() < 0) << "Structure element was not found" << std::endl;
    for (auto elem_id : structure_element_list){
        mrModelPart.GetElement(elem_id).Set(STRUCTURE);
    }
    for (auto elem_id : set_kutta_element_list){
        auto& kutta_elem = mrModelPart.GetElement(elem_id);
    //     bool is_wing_tip = false;
        auto& r_geometry = kutta_elem.GetGeometry();

    //     for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {

    //         if (r_geometry[i_node].GetValue(WING_TIP))
    //         {
    //             is_wing_tip = true;
    //         }
    //     }
    //     KRATOS_WATCH(elem_id)
    //     if (is_wing_tip){
        unsigned int kutta_nodes = 0;
        for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {
//
            if (r_geometry[i_node].GetValue(KUTTA))
            {
                kutta_nodes++;
            }
        }
        if (kutta_nodes > 1) {
            kutta_elem.SetValue(KUTTA, true);
            kutta_elem.SetValue(WAKE, false);

        }
    //         kutta_elem.SetValue(WAKE, true);
    //         // structure_element_id = kutta_elem.Id();
    //         kutta_elem.Set(STRUCTURE);
    //         // KRATOS_WATCH(kutta_elem.Id())
    //         // for (unsigned int i_node= 0; i_node < r_geometry.size(); i_node++) {

    //         //     if (r_geometry[i_node].GetValue(KUTTA) && r_geometry[i_node].GetValue(WAKE_DISTANCE)<0.0)
    //         //     {
    //         //         r_geometry[i_node].SetValue(WING_TIP, true);
    //         //     }
    //         // }
        // }
    }



    mrModelPart.RemoveSubModelPart("deactivated_model_part");
    KRATOS_ERROR_IF_NOT(trailing_edge_node_list.size()) << "No trailing edge nodes were found" << std::endl;

    mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);
}

void DefineEmbeddedWakeProcess::MarkKuttaWakeElements(){

}
}// Namespace Kratos

