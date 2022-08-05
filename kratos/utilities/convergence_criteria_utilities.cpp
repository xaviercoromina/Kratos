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

// System includes

// External includes

// Project includes
#include "utilities/convergence_criteria_utilities.h"

namespace Kratos
{
namespace ConvergenceCriteriaUtilities
{

std::vector<const VariableData*> GenerateVariableDataVector(const ConvergenceVariableListType& rList)
{
    std::size_t i = 0;
    std::vector<const VariableData*> aux_vect(rList.size());
    for (const auto& r_tup : rList) {
        aux_vect[i++] = std::get<0>(r_tup);
    }
    return aux_vect;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<double> GenerateRatioToleranceVector(const ConvergenceVariableListType& rList)
{
    std::size_t i = 0;
    std::vector<double> aux_vect(rList.size());
    for (const auto& r_tup : rList) {
        aux_vect[i++] = std::get<1>(r_tup);
    }
    return aux_vect;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<double> GenerateAbsToleranceVector(const ConvergenceVariableListType& rList)
{
    std::size_t i = 0;
    std::vector<double> aux_vect(rList.size());
    for (const auto& r_tup : rList) {
        aux_vect[i++] = std::get<2>(r_tup);
    }
    return aux_vect;
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<std::size_t, std::size_t> GenerateLocalKeyMap(const ConvergenceVariableListType& rList)
{
    std::size_t local_key = 0;
    std::unordered_map<std::size_t, std::size_t> aux_map;
    for (const auto& r_tup : rList) {
        const auto* p_var_data = std::get<0>(r_tup);
        if (aux_map.find(p_var_data->Key()) != aux_map.end()) {
            KRATOS_ERROR << "Convergence variable " << p_var_data->Name() << " is repeated. Check the input convergence variable list." << std::endl;
        } else {
            KRATOS_ERROR_IF(p_var_data->IsComponent()) << "Trying to check convergence with the " << p_var_data->Name() << " component variable. Use the corresponding vector one." << std::endl;
            aux_map[p_var_data->Key()] = local_key++;
        }
    }
    return aux_map;
}

/***********************************************************************************/
/***********************************************************************************/

ConvergenceVariableListType GenerateConvergenceVariableListFromParameters(Kratos::Parameters ThisParameters)
{
    // Iterate over variables
    ConvergenceVariableListType aux_list;
    if (!ThisParameters.Has("convergence_variables_list")) return aux_list;
    Kratos::Parameters convergence_variables_list = ThisParameters["convergence_variables_list"];
    for (auto param : convergence_variables_list) {
        if (param.Has("variable")) {
            const std::string& r_variable_name = param["variable"].GetString();

            // Variable pointer
            const VariableData* p_variable = KratosComponents<Variable<double>>::Has(r_variable_name) ? dynamic_cast<const VariableData*>(&KratosComponents<Variable<double>>::Get(r_variable_name)) : dynamic_cast<const VariableData*>(&KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));

            // Tolerances
            const double rel_tol = param.Has("relative_tolerance") ? param["relative_tolerance"].GetDouble() : 1.0e-4;
            const double abs_tol = param.Has("absolute_tolerance") ? param["absolute_tolerance"].GetDouble() : 1.0e-9;

            // Push back list
            aux_list.push_back(std::make_tuple(p_variable, rel_tol, abs_tol));
        }
    }

    return aux_list;
}

/***********************************************************************************/
/***********************************************************************************/

void KRATOS_API(KRATOS_CORE) OutputConvergenceStatus(
    const std::tuple<std::vector<double>, std::vector<double>>& rConvergenceNorms,
    const int VariableSize,
    const std::vector<const VariableData*>& rVariableDataVector,
    const std::vector<double>& rRatioToleranceVector,
    const std::vector<double>& rAbsToleranceVector,
    std::unordered_map<std::size_t, std::size_t>& rLocalKeyMap,
    const std::size_t EchoLevel
    )
{
    if (EchoLevel > 0) {
        const auto& r_var_ratio = std::get<0>(rConvergenceNorms);
        const auto& r_var_abs = std::get<1>(rConvergenceNorms);

        std::ostringstream stringbuf;
        stringbuf << "CONVERGENCE CHECK:\n";

        const int max_length_var_name = (*std::max_element(rVariableDataVector.begin(), rVariableDataVector.end(), [](const VariableData* p_var_data_1, const VariableData* p_var_data_2){
            return p_var_data_1->Name().length() < p_var_data_2->Name().length();
        }))->Name().length();

        for(int i = 0; i < VariableSize; i++) {
            const auto r_var_data = rVariableDataVector[i];
            const int key_map = rLocalKeyMap[r_var_data->Key()];
            const std::string space_str(max_length_var_name-r_var_data->Name().length(), ' ');
            stringbuf << " " << r_var_data->Name() << space_str <<" : ratio = " << r_var_ratio[key_map] << "; exp.ratio = " << rRatioToleranceVector[key_map] << " abs = " << r_var_abs[key_map] << " exp.abs = " << rAbsToleranceVector[key_map] << "\n";
        }
        KRATOS_INFO("") << stringbuf.str();
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool KRATOS_API(KRATOS_CORE) CheckConvergence(
    const std::tuple<std::vector<double>, std::vector<double>>& rConvergenceNorms,
    const int VariableSize,
    const std::vector<const VariableData*>& rVariableDataVector,
    const std::vector<double>& rRatioToleranceVector,
    const std::vector<double>& rAbsToleranceVector,
    std::unordered_map<std::size_t, std::size_t>& rLocalKeyMap,
    const std::size_t EchoLevel
    )
{
    bool is_converged = true;
    const auto& r_var_ratio = std::get<0>(rConvergenceNorms);
    const auto& r_var_abs = std::get<1>(rConvergenceNorms);

    for (int i = 0; i < VariableSize; ++i) {
        const auto r_var_data = rVariableDataVector[i];
        const std::size_t key_map = rLocalKeyMap[r_var_data->Key()];
        is_converged &= r_var_ratio[key_map] <= rRatioToleranceVector[key_map] || r_var_abs[key_map] <= rAbsToleranceVector[key_map];
    }

    // Note that this check ensures that all the convergence variables fulfil either the relative or the absolute criterion
    if (is_converged) {
        KRATOS_INFO_IF("", EchoLevel > 0) << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
        return true;
    } else {
        return false;
    }
}

} // namespace DelaunatorUtilities
} // namespace Kratos

#undef REAL
