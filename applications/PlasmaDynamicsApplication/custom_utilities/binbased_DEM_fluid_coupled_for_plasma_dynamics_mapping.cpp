#include "binbased_DEM_fluid_coupled_for_plasma_dynamics_mapping.h"
#include "plasma_dynamics_application.h"
#include "plasma_dynamics_application_variables.h"


namespace Kratos
{

/// Interpolate fluid data onto the DEM model part.
/**
  * @param r_fluid_model_part: the origin model part from which to project
  * @param r_dem_model_part: the destination model part of which we want to interpolate its nodal values
  * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
  * @see binbased_nodes_in_element_locator
*/
// data_to_project to DEM mesh = alpha * new_data + (1 - alpha) * old_data

template <std::size_t TDim, typename TBaseTypeOfIonParticle>
void BinBasedDEMFluidCoupledPDMapping<TDim, TBaseTypeOfIonParticle>::InterpolateFromFluidMesh(
        ModelPart& r_fluid_model_part,
        ModelPart& r_dem_model_part,
        Parameters& parameters,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha)
{
    KRATOS_TRY

    // InterpolateFromFluidMesh 
    BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateFromFluidMesh(
        ModelPart& r_fluid_model_part,
        ModelPart& r_dem_model_part,
        Parameters& parameters,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha);

    GetElectricFieldValue(ModelPart& r_fluid_model_part);

    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    #pragma omp parallel for firstprivate(results, shape_function_values_at_point)
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        Node<3>::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            const bool element_located = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(),
                                                                              shape_function_values_at_point,
                                                                              p_element,
                                                                              results.begin(),
                                                                              max_results);
            // interpolating the variables

            if (element_located){

                p_particle->Set(INSIDE, true);

                for (unsigned int j = 0; j != mDEMCouplingVariables.size(); ++j){
                    ProjectPlasmaDynamics(p_element,
                            shape_function_values_at_point,
                            p_particle,
                            mDEMCouplingVariables[j],
                            alpha);                            
                }
            }

            else {
                p_particle->Set(INSIDE, false);
            }
        }
    }

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfIonParticle>
void BinBasedDEMFluidCoupledPDMapping<TDim, TBaseTypeOfIonParticle>::ProjectPlasmaDynamics(Element::Pointer p_elem,
             const Vector& N,
             Node<3>::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha)
{
    if (*r_destination_variable == ELECTRIC_FIELD_PROJECTED_TO_PARTICLE){
        BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::Interpolate(
            p_elem, N, p_node, ELECTRIC_FIELD, ELECTRIC_FIELD_PROJECTED_TO_PARTICLE, alpha);
    }

}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfIonParticle>
void BinBasedDEMFluidCoupledPDMapping<TDim, TBaseTypeOfIonParticle>::GetElectricFieldValue(ModelPart& r_fluid_model_part)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
		array_1d<double, 3 >& grad_phi = i_node->FastGetSolutionStepValue(NODAL_PHI_GRADIENT);
        array_1d<double, 3 >& electric_field = - grad_phi;

        i_node->SetSolutionStepValue(ELECTRIC_FIELD, electric_field);

    }
}
//***************************************************************************************************************
//***************************************************************************************************************

// Explicit instantiations
template class BinBasedDEMFluidCoupledPDMapping<2, IonParticle>;
template class BinBasedDEMFluidCoupledPDMapping<3, IonParticle>;

}//namespace Kratos