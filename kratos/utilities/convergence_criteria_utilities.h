//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONVERGENCE_CRITERIA_UTILITIES)
#define KRATOS_CONVERGENCE_CRITERIA_UTILITIES

// System includes
#include <vector>

// External includes

// Project includes
#include "containers/variable_data.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // Convergence variable list type
    typedef std::vector<std::tuple<const VariableData*, double, double>> ConvergenceVariableListType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @namespace ConvergenceCriteriaUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities for the converence criterias
 * @author Vicente Mataix Ferrandiz
 */
namespace ConvergenceCriteriaUtilities
{
    /**
     * @brief This generates the variable data vector
     * @param rList The list of variables
     * @return The variable data vector
     */
    std::vector<const VariableData*> KRATOS_API(KRATOS_CORE) GenerateVariableDataVector(const ConvergenceVariableListType& rList);

    /**
     * @brief This generates the ratio tolerance vector
     * @param rList The list of variables
     * @return The ratio tolerance vector
     */
    std::vector<double> KRATOS_API(KRATOS_CORE) GenerateRatioToleranceVector(const ConvergenceVariableListType& rList);

    /**
     * @brief This generates the absolute tolerance vector
     * @param rList The list of variables
     * @return The absolute tolerance vector
     */
    std::vector<double> KRATOS_API(KRATOS_CORE) GenerateAbsToleranceVector(const ConvergenceVariableListType& rList);

    /**
     * @brief This generates the local key map
     * @param rList The list of variables
     * @return The local key map
     */
    std::unordered_map<std::size_t, std::size_t> KRATOS_API(KRATOS_CORE) GenerateLocalKeyMap(const ConvergenceVariableListType& rList);

    /**
     * @brief This method generates the list of variables from Parameters
     * @param ThisParameters Input parameters
     * @return List of variables considered as input
     */
    ConvergenceVariableListType KRATOS_API(KRATOS_CORE) GenerateConvergenceVariableListFromParameters(Kratos::Parameters ThisParameters);

    /**
     * @brief Method to output the convergence status
     * @details This method prints the convergence status to the screen for each one of the checked variables
     * @param rConvergenceNorms Tuple containing the absolute and relative convergence values
     * @param VariableSize The number of variables
     * @param rVariableDataVector The variable data vector
     * @param rRatioToleranceVector The ratio tolerance vector
     * @param rAbsToleranceVector The absolute tolerance vector
     * @param rLocalKeyMap The local key map
     * @param EchoLevel The level of verbosity
     */
    void KRATOS_API(KRATOS_CORE) OutputConvergenceStatus(
        const std::tuple<std::vector<double>, std::vector<double>>& rConvergenceNorms,
        const int VariableSize,
        const std::vector<const VariableData*>& rVariableDataVector,
        const std::vector<double>& rRatioToleranceVector,
        const std::vector<double>& rAbsToleranceVector,
        std::unordered_map<std::size_t, std::size_t>& rLocalKeyMap,
        const std::size_t EchoLevel = 0
        );

    /**
     * @brief Method to check convergence
     * @details This method checks the convergence of the provided norms with the user-defined tolerances
     * @param rConvergenceNorms Tuple containing the absolute and relative convergence values
     * @param VariableSize The number of variables
     * @param rVariableDataVector The variable data vector
     * @param rRatioToleranceVector The ratio tolerance vector
     * @param rAbsToleranceVector The absolute tolerance vector
     * @param rLocalKeyMap The local key map
     * @param EchoLevel The level of verbosity
     * @return true Convergence is satisfied
     * @return false Convergence is not satisfied
     */
    bool KRATOS_API(KRATOS_CORE) CheckConvergence(
        const std::tuple<std::vector<double>, std::vector<double>>& rConvergenceNorms,
        const int VariableSize,
        const std::vector<const VariableData*>& rVariableDataVector,
        const std::vector<double>& rRatioToleranceVector,
        const std::vector<double>& rAbsToleranceVector,
        std::unordered_map<std::size_t, std::size_t>& rLocalKeyMap,
        const std::size_t EchoLevel = 0
        );

}; // namespace ConvergenceCriteriaUtilities
}  // namespace Kratos
#endif /* KRATOS_CONVERGENCE_CRITERIA_UTILITIES defined */
