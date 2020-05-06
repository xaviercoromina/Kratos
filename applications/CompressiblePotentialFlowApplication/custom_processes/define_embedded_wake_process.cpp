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

    KRATOS_CATCH("");
}


void DefineEmbeddedWakeProcess::ComputeDistanceToWake(){

    // #pragma omp parallel for
    // for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
    //     auto it_elem = mrModelPart.ElementsBegin() + i;
    //     auto geometry_elemental_distances = it_elem->GetValue(ELEMENTAL_DISTANCES);
    //     it_elem->SetValue(GEOMETRY_ELEMENTAL_DISTANCES, geometry_elemental_distances);
    //     it_elem->SetValue(WAKE, false);
    //     it_elem->SetValue(KUTTA, false);
    //     if (it_elem->Is(TO_SPLIT)) {
    //         it_elem->Set(BOUNDARY);
    //     }
    // }

    CalculateDiscontinuousDistanceToSkinProcess<2> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();

    // #pragma omp parallel for
    // for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
    //     auto it_node = mrModelPart.NodesBegin() + i;
    //     it_node->SetValue(WAKE_DISTANCE, 0.0);
    //     it_node->SetValue(KUTTA, 0);
    //     it_node->SetValue(WAKE, 0);
    // }
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
        for(unsigned int i_node = 0; i_node<3; i_node++){
            geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);
        // if (is_embedded){
        //     ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_distances));
        //     // Computing Normal
        //     std::vector<Vector> cut_normal;
        //     pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
        //     double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
        //     auto unit_normal = cut_normal[0]/norm_normal;
        //     it_elem->SetValue(VELOCITY_LOWER,unit_normal);
        // }
        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element) {
            if (is_embedded){
                #pragma omp critical
                {
                    deactivated_elements_id_list.push_back(it_elem->Id());
                }
                // it_elem->Set(ACTIVE, false);
                // it_elem->Set(KUTTA, false);
                auto& r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WING_TIP, true);
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].UnSetLock();
                }
            }
            else{
                it_elem->SetValue(WAKE, true);
                auto& r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].SetValue(WAKE, true);
                    r_geometry[i].UnSetLock();
                }
            }
        }
    }
    deactivated_model_part.AddElements(deactivated_elements_id_list);
}

void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){

    // double max_distance = 0.0;
    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    Element::Pointer p_max_elem;

    auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];

    // Find furthest deactivated element to the wake origin
    // Find deactivated element with Dim nodes shared with a wake element
    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
        std::size_t wake_nodes_counter = 0;
        auto& r_geometry = it_elem->GetGeometry();
        for (unsigned int i_node= 0; i_node < it_elem->GetGeometry().size(); i_node++) {
            if (r_geometry[i_node].GetValue(WAKE)) {
                wake_nodes_counter++;
            }
        }
        if (wake_nodes_counter == 2) {
            p_max_elem = mrModelPart.pGetElement(it_elem->Id());
        }
    }
    auto angle_in_deg = -1*mrModelPart.GetProcessInfo()[ROTATION_ANGLE];
    BoundedVector<double, 3> wake_direction;
    wake_direction[0] = cos(angle_in_deg*Globals::Pi/180);
    wake_direction[1] = sin(angle_in_deg*Globals::Pi/180);
    wake_direction[2] = 0.0;
    KRATOS_WATCH(wake_origin)
    KRATOS_WATCH(wake_direction)

    // Mark nodes of the furthest deactivated element and store its neighbour elements
    std::size_t max_node_i = -1;
    double max_x = -1e10;
    auto& r_max_element_geometry=p_max_elem->GetGeometry();
    p_max_elem->SetValue(KUTTA, true);
    KRATOS_WATCH(p_max_elem->Id())
    for (unsigned int i_node= 0; i_node < r_max_element_geometry.size(); i_node++) {
        // r_max_element_geometry[i_node].SetValue(AIRFOIL,true);
        // BoundedVector<double, 3> nodal_position_vector = r_max_element_geometry[i_node].Coordinates() - wake_origin;
        BoundedVector<double, 3> nodal_position_vector = r_max_element_geometry[i_node].Coordinates();
        KRATOS_WATCH(nodal_position_vector)
        KRATOS_WATCH(inner_prod(nodal_position_vector, wake_direction))
        double projection = std::abs(inner_prod(nodal_position_vector, wake_direction));
        KRATOS_WATCH(projection)
        // if (r_max_element_geometry[i_node].X()>max_x) {
        if (projection>max_x) {
            max_x = projection;
            max_node_i = i_node;
        }

        const GlobalPointersVector<Element>& r_node_elem_candidates = r_max_element_geometry[i_node].GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
            mKuttaWakeElementCandidates.push_back(r_node_elem_candidates(j));
        }
    }

    KRATOS_ERROR_IF(max_node_i < 0) << "No trailing edge nodes were found" << std::endl;
    r_max_element_geometry[max_node_i].SetValue(AIRFOIL, true);
    mrModelPart.RemoveSubModelPart("deactivated_model_part");

    std::vector<std::size_t> trailing_edge_node_list;
    trailing_edge_node_list.push_back(r_max_element_geometry[max_node_i].Id());

    if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
        mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");

    std::sort(trailing_edge_node_list.begin(),
              trailing_edge_node_list.end());
    mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);
}

void DefineEmbeddedWakeProcess::MarkKuttaWakeElements(){

    //TO-DO avoid using nodal distances
    double te_wake_distance;
    for (int i = 0; i < static_cast<int>(mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").Nodes().size()); i++) {
        auto it_node = mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").NodesBegin() + i;
        te_wake_distance=it_node->GetValue(WAKE_DISTANCE);

    }

    KRATOS_WATCH(te_wake_distance)
    std::vector<std::size_t> potentially_kutta_elements;
    Element::Pointer p_structure_element;

    // Find elements that touch the furthest deactivated element and that are part of the wake.
    for (std::size_t i = 0; i < mKuttaWakeElementCandidates.size(); i++)
    {
        auto& r_geometry = mKuttaWakeElementCandidates[i].GetGeometry();
        if (mKuttaWakeElementCandidates[i].GetValue(WAKE) && mKuttaWakeElementCandidates[i].Is(ACTIVE)) {
            std::size_t matching_sign_counter = 0;
            for (std::size_t i_node= 0; i_node < r_geometry.size(); i_node++) {
                //TO-DO CHECK IF product is 0
                if(te_wake_distance*r_geometry[i_node].GetValue(WAKE_DISTANCE)>0.0){
                    matching_sign_counter++;
                }
            }
            if (matching_sign_counter == 2) {
                mKuttaWakeElementCandidates[i].Set(STRUCTURE);
                p_structure_element = mrModelPart.pGetElement(mKuttaWakeElementCandidates[i].Id());
            }
            else {
                potentially_kutta_elements.push_back(mKuttaWakeElementCandidates[i].Id());
            }
        }
    }
    std::size_t wing_tip_structure_counter = 0;
    for (std::size_t i = 0; i < p_structure_element->GetGeometry().size(); i++) {
        if (p_structure_element->GetGeometry()[i].GetValue(WING_TIP)) {
            wing_tip_structure_counter++;
        }
    }

    for (std::size_t i = 0; i < potentially_kutta_elements.size(); i++) {
        if (wing_tip_structure_counter<2) {
            mrModelPart.GetElement(potentially_kutta_elements[i]).SetValue(WAKE, false);
            mrModelPart.GetElement(potentially_kutta_elements[i]).SetValue(KUTTA, true);
        }
    }

}
}// Namespace Kratos
