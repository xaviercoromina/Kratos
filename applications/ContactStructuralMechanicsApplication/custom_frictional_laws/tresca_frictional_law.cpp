// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_frictional_laws/tresca_frictional_law.h"

namespace Kratos
{
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
double TrescaFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetThresholdValue(
    const NodeType& rNode,
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_properties = rCondition.GetProperties();
    if (r_properties.Has(TRESCA_FRICTION_THRESHOLD)) {
        return r_properties.GetValue(TRESCA_FRICTION_THRESHOLD);
    } else if (rNode.Has(TRESCA_FRICTION_THRESHOLD)) {
        return rNode.GetValue(TRESCA_FRICTION_THRESHOLD);
    } else if (rCurrentProcessInfo.Has(TRESCA_FRICTION_THRESHOLD)) {
        return rCurrentProcessInfo.GetValue(TRESCA_FRICTION_THRESHOLD);
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
typename TrescaFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::DerivativesArray TrescaFrictionalLaw<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDerivativesThresholdArray(
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const DerivativeDataType& rDerivativeData,
    const MortarConditionMatrices& rMortarConditionMatrices
    )
{
    DerivativesArray derivatives_threshold_values;

    // Threshold is constant, derivative is always zero
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        for (std::size_t i_der = 0; i_der < TDim * (2 * TNumNodes + TNumNodesMaster); ++i_der) {
            derivatives_threshold_values[i_der][i_node] = 0.0;
        }
    }

    return derivatives_threshold_values;
}

/***********************************************************************************/
/***********************************************************************************/

template class TrescaFrictionalLaw<2, 2, false, 2>;
template class TrescaFrictionalLaw<3, 3, false, 3>;
template class TrescaFrictionalLaw<3, 4, false, 4>;
template class TrescaFrictionalLaw<3, 3, false, 4>;
template class TrescaFrictionalLaw<3, 4, false, 3>;
template class TrescaFrictionalLaw<2, 2, true,  2>;
template class TrescaFrictionalLaw<3, 3, true,  3>;
template class TrescaFrictionalLaw<3, 4, true,  4>;
template class TrescaFrictionalLaw<3, 3, true,  4>;
template class TrescaFrictionalLaw<3, 4, true,  3>;

}  // namespace Kratos.

