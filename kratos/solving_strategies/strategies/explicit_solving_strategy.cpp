//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;

template <class TSparseSpace, class TDenseSpace>
void ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>::AssignSettings(const Parameters ThisParameters)
{
    const bool rebuild_level = ThisParameters["rebuild_level"].GetInt();
    const bool move_mesh_flag = ThisParameters["move_mesh_flag"].GetBool();
    SetMoveMeshFlag(move_mesh_flag);
    SetRebuildLevel(rebuild_level);

    // Setting up the default builder and solver
    const std::string& r_name = ThisParameters.Has("explicit_builder_settings") ? ThisParameters["explicit_builder_settings"].Has("name") ? ThisParameters["explicit_builder_settings"]["name"].GetString() : "explicit_builder" : "explicit_builder";
    if (KratosComponents<ExplicitBuilderType>::Has( r_name )) {
        // Defining the builder and solver
        mpExplicitBuilder = KratosComponents<ExplicitBuilderType>::Get(r_name).Create(ThisParameters["explicit_builder_settings"]);
    } else {
        KRATOS_ERROR << "Trying to construct explicit builder with name= " << r_name << std::endl <<
                        "Which does not exist. The list of available options (for currently loaded applications) are: " << std::endl <<
                        KratosComponents<ExplicitBuilderType>() << std::endl;
    }
}

//NOTE: here we must create persisting objects for the strategies
static ExplicitSolvingStrategyType msExplicitSolvingStrategy;

template<>
std::vector<Internals::RegisteredPrototypeBase<ExplicitSolvingStrategyType>> ExplicitSolvingStrategyType::msPrototypes{
    Internals::RegisteredPrototype<ExplicitSolvingStrategyType, ExplicitSolvingStrategyType>(ExplicitSolvingStrategyType::Name(), msExplicitSolvingStrategy)};

///@}

} /* namespace Kratos.*/
