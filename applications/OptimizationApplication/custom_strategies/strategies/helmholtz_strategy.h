//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/condition_number_utility.h"
#include "utilities/variable_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class HelmholtzStrategy: public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type definitions
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(HelmholtzStrategy);

    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    using TBuilderAndSolverType = typename BaseType::TBuilderAndSolverType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using SchemeType = Scheme<TSparseSpace, TDenseSpace>;

    ///@}
    ///@name Life Cycle
    ///@{

    /* This class is not required. */
    HelmholtzStrategy(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pNewLinearSolver,
        const bool ReformDofSetAtEachStep = false,
        const bool ComputeReactions = false,
        const int EchoLevel = 0,
        const double PoissonRatio = 0.3)
        : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;
        mComputeReactions = ComputeReactions;
        mEchoLevel = EchoLevel;

        bool calculate_norm_dx_flag = false;

        // TODO: This needs to be passed from python level to support MPI
        mpSheme = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());

        // TODO: This needs to be passed from python level to support MPI
        mpBuliderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));

        // TODO: This needs to be passed from python level to support MPI
        mpStrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
                                                                BaseType::GetModelPart(),
                                                                mpSheme,
                                                                mpBuliderAndSolver,
                                                                mComputeReactions,
                                                                mReformDofSetAtEachStep,
                                                                calculate_norm_dx_flag));

        mpStrategy->SetEchoLevel(mEchoLevel);

        KRATOS_CATCH("")
    }

    HelmholtzStrategy(const HelmholtzStrategy &Other) = delete;

    ~HelmholtzStrategy() override = default;

    ///@}
    ///@name Public operations
    ///@{

    double Solve() override
    {
        KRATOS_TRY;

        VariableUtils().UpdateCurrentToInitialConfiguration(
            BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

        // Solve for the mesh movement
        mpStrategy->Solve();

        // Clearing the system if needed
        if (mReformDofSetAtEachStep) {
            mpStrategy->Clear();
        }

        return 0.0;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{

    TSystemMatrixType &GetSystemMatrix() override
    {
        return mpStrategy->GetSystemMatrix();
    }

    TSystemVectorType &GetSystemVector() override
    {
        return mpStrategy->GetSystemVector();
    }

    TSystemVectorType &GetSolutionVector() override
    {
        return mpStrategy->GetSolutionVector();
    }

    typename BaseType::Pointer GetStrategy()
    {
        return mpStrategy;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    typename SchemeType::Pointer mpSheme;

    typename BaseType::Pointer mpStrategy;

    typename TBuilderAndSolverType::Pointer mpBuliderAndSolver;

    int mEchoLevel;

    bool mReformDofSetAtEachStep;

    bool mComputeReactions;

    ///@}

};

} // namespace Kratos
