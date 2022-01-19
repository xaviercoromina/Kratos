// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Constants for the computation of the stress
        const double E = p_elem_prop->GetValue(YOUNG_MODULUS);
        const double NU = p_elem_prop->GetValue(POISSON_RATIO);



        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0000000000, 0.0100000000, 0.0100000000);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.0000000000, 0.0000000000, 0.0100000000);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0000000000, 0.0100000000, 0.0000000000);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0100000000, 0.0100000000, 0.0100000000);
        auto p_node_5 = r_model_part.CreateNewNode(5, 0.0100000000, 0.0000000000, 0.0100000000);
        auto p_node_6 = r_model_part.CreateNewNode(6, 0.0100000000, 0.0100000000, 0.0000000000);
        auto p_node_7 = r_model_part.CreateNewNode(7, 0.0000000000, 0.0000000000, 0.0000000000);
        auto p_node_8 = r_model_part.CreateNewNode(8, 0.0100000000, 0.0000000000, 0.0000000000);

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }

        std::vector<ModelPart::IndexType> element_nodes {4,1,3,6,5,2,7,8};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement3D8N", 1, element_nodes, p_elem_prop);

        p_element->Initialize(r_model_part.GetProcessInfo());

        // Define a matrix A and impose that the displacement on each node is u = A*x0
        Matrix A0 = ZeroMatrix(3,3);
        A0(0,0) = 2e-2;
        A0(0,1) = 5e-2;
        A0(1,0) = 5e-2;
        A0(1,1) = -1e-2;
        A0(2,2) = 1.0;

        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        double factor = 1.0;
        for(unsigned int step=0; step<5; step++)
        {
            factor+=1.0;
            Matrix A = factor*A0;
            for (auto& r_node : r_model_part.Nodes()){
                noalias(r_node.GetSolutionStepValue(DISPLACEMENT)) = prod(A, r_node.GetInitialPosition());
                r_node.Coordinates() = r_node.GetInitialPosition() + r_node.GetSolutionStepValue(DISPLACEMENT); //here i update the node coordinates
            }

            p_element->InitializeSolutionStep(const_procinfo_ref);
            p_element->InitializeNonLinearIteration(const_procinfo_ref);
            p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
            p_element->FinalizeNonLinearIteration(const_procinfo_ref);
            p_element->FinalizeSolutionStep(const_procinfo_ref);




        }
    }
}
}
