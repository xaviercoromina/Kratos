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
#include "custom_frictional_laws/frictional_law_with_derivative.h"

namespace Kratos
{
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
array_1d<double, TNumNodes> FrictionalLawWithDerivative<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetThresholdArray(
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double, TNumNodes> threshold_values;

    const auto& r_geometry = rCondition.GetGeometry();
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        threshold_values[i_node] = this->GetThresholdValue(r_node, rCondition, rCurrentProcessInfo);
    }

    return threshold_values;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
typename FrictionalLawWithDerivative<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::DerivativesArray FrictionalLawWithDerivative<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDerivativesThresholdArray(
    const PairedCondition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const DerivativeDataType& rDerivativeData,
    const MortarConditionMatrices& rMortarConditionMatrices
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetDerivativesThresholdArray, check your frictional law declaration" << std::endl;

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

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AddKratosComponent(std::string const& Name, FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster> const& ThisComponent)
{
    KratosComponents<FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>>::Add(Name, ThisComponent);
}

/***********************************************************************************/
/***********************************************************************************/

template class KratosComponents<FrictionalLawWithDerivative<2, 2, false, 2>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 3, false, 3>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 4, false, 4>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 3, false, 4>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 4, false, 3>>;
template class KratosComponents<FrictionalLawWithDerivative<2, 2, true,  2>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 3, true,  3>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 4, true,  4>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 3, true,  4>>;
template class KratosComponents<FrictionalLawWithDerivative<3, 4, true,  3>>;

/***********************************************************************************/
/***********************************************************************************/

template class FrictionalLawWithDerivative<2, 2, false, 2>;
template class FrictionalLawWithDerivative<3, 3, false, 3>;
template class FrictionalLawWithDerivative<3, 4, false, 4>;
template class FrictionalLawWithDerivative<3, 3, false, 4>;
template class FrictionalLawWithDerivative<3, 4, false, 3>;
template class FrictionalLawWithDerivative<2, 2, true,  2>;
template class FrictionalLawWithDerivative<3, 3, true,  3>;
template class FrictionalLawWithDerivative<3, 4, true,  4>;
template class FrictionalLawWithDerivative<3, 3, true,  4>;
template class FrictionalLawWithDerivative<3, 4, true,  3>;

}  // namespace Kratos.

