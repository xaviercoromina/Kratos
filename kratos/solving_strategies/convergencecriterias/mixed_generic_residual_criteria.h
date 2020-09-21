//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Riccardo Rossi
//                   Carlos Roig
//                   Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//

#ifndef KRATOS_MIXED_GENERIC_RESIDUAL_CRITERIA_H
#define	KRATOS_MIXED_GENERIC_RESIDUAL_CRITERIA_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "convergence_criteria.h"
#include "utilities/constraint_utilities.h"
#include "utilities/convergence_criteria_utilities.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class MixedGenericResidualCriteria
 * @brief Convergence criteria for mixed vector-scalar problems.
 * @ingroup KratosCore
 * @details This class implements a convergence control based on the check of the residual of a nodal vector variable and a nodal scalar variable. The error is evaluated separately for each of them, and relative and absolute tolerances for both must be specified.
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
template< class TSparseSpace, class TDenseSpace >
class MixedGenericResidualCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MixedGenericResidualCriteria);

    /// The definition of the base ConvergenceCriteria
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    /// The definition of the current class
    typedef MixedGenericResidualCriteria< TSparseSpace, TDenseSpace > ClassType;

    /// The data type
    typedef typename BaseType::TDataType TDataType;

    /// The dofs array type
    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// The sparse matrix type
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    /// The dense vector type
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /// The convergence variable list type
    typedef std::vector<std::tuple<const VariableData*, TDataType, TDataType>> ConvergenceVariableListType;

    /// Definition of the IndexType
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    explicit MixedGenericResidualCriteria()
        : BaseType(),
          mVariableSize(0)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit MixedGenericResidualCriteria(Kratos::Parameters ThisParameters)
        : MixedGenericResidualCriteria(ConvergenceCriteriaUtilities::GenerateConvergenceVariableListFromParameters(ThisParameters))
    {
        this->mActualizeRHSIsNeeded = true;
    }

    /**
     * @brief Construct a new Mixed Generic Criteria object
     * Construct the mixed generic convergence criteria from a convergence variables list.
     * The convergence variable list contains for each variable the variable itself as well as the corresponding relative and absolute tolerances.
     * @param rConvergenceVariablesList List containing tuples with the convergence variables to be checked. The tuples are set as <Variable, relative tolerance, absolute tolerance>
     */
    MixedGenericResidualCriteria(const ConvergenceVariableListType& rConvergenceVariablesList)
        : BaseType()
        , mVariableSize([&] (const ConvergenceVariableListType& rList) -> int {return rList.size();} (rConvergenceVariablesList))
        , mVariableDataVector(ConvergenceCriteriaUtilities::GenerateVariableDataVector(rConvergenceVariablesList))
        , mRatioToleranceVector(ConvergenceCriteriaUtilities::GenerateRatioToleranceVector(rConvergenceVariablesList))
        , mAbsToleranceVector(ConvergenceCriteriaUtilities::GenerateAbsToleranceVector(rConvergenceVariablesList))
        , mLocalKeyMap(ConvergenceCriteriaUtilities::GenerateLocalKeyMap(rConvergenceVariablesList))
    {
        mInitialResidualNormVector = std::vector<TDataType>(mVariableSize, 0.0);
        mCurrentResidualNormVector = std::vector<TDataType>(mVariableSize, 0.0);
        mReferenceDispNormVector   = std::vector<TDataType>(mVariableSize, 0.0);

        this->mActualizeRHSIsNeeded = true;
    }

    /// Destructor.
    ~MixedGenericResidualCriteria() override
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Check if we are solving for something
        if (TSparseSpace::Size(rb) != 0) {
            // Calculate the convergence ratio and absolute norms
            const auto convergence_norms = CalculateConvergenceNorms(rModelPart, rDofSet, rb);

            // Output convergence status
            OutputConvergenceStatus(convergence_norms);

            // Check convergence
            return ConvergenceCriteriaUtilities::CheckConvergence(convergence_norms, mVariableSize, mVariableDataVector, mRatioToleranceVector, mAbsToleranceVector, mLocalKeyMap, this->GetEchoLevel());
        } else {
            // Case in which all the DOFs are constrained!
            return true;
        }
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
        KRATOS_ERROR_IF(rModelPart.IsDistributed() && rModelPart.NumberOfMasterSlaveConstraints() > 0) << "This Criteria does not yet support constraints in MPI!" << std::endl;
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);

        // Filling mActiveDofs when MPC exist
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            ConstraintUtilities::ComputeActiveDofs(rModelPart, mActiveDofs, rDofSet);
        }

        SizeType size_residual;
        CalculateResidualNorm(rModelPart, mInitialResidualNormVector, size_residual, rDofSet, rb);
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "mixed_generic_residual_criteria";
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MixedGenericResidualCriteria";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
protected:
    ///@name Protected Static Member Variables
    ///@{


    ///@}
    ///@name Protected Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Get the Variable Size object
     * Get the number of variables to be checked
     * @return const int Number of variables to check
     */
    int GetVariableSize() const
    {
        return mVariableSize;
    }

    /**
     * @brief Get the Variable Data Vector object
     * Get the member vector that stores pointers to the variables to check
     * @return std::vector<VariableData*> Vector containing pointers to the variables to check
     */
    std::vector<const VariableData*> GetVariableDataVector() const
    {
        return mVariableDataVector;
    }

    /**
     * @brief Get the Ratio Tolerance Vector object
     * Get the member vector containing the ratio tolerances for each variable to check
     * @return std::vector<TDataType> Vector containing the ratio tolerances
     */
    std::vector<TDataType> GetRatioToleranceVector() const
    {
        return mRatioToleranceVector;
    }

    /**
     * @brief Get the Abs Tolerance Vector object
     * Get the member vector containing the absolute tolerances for each variable to check
     * @return std::vector<TDataType> Vector containing the absolute tolerances
     */
    std::vector<TDataType> GetAbsToleranceVector() const
    {
        return mAbsToleranceVector;
    }

    /**
     * @brief Get the Local Key Map object
     * Returns a reference to the variable key local map
     * @return std::unordered_map<IndexType, IndexType>& Reference to the local key map
     */
    std::unordered_map<IndexType, IndexType>& GetLocalKeyMap()
    {
        return mLocalKeyMap;
    }

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void CalculateResidualNorm(
        ModelPart& rModelPart,
        std::vector<TDataType>& rResidualSolutionNorm,
        SizeType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        std::vector<TDataType> residual_solution_norm = std::vector<TDataType>(mVariableSize, 0.0);
        SizeType dof_num = 0;

        // Auxiliar values
        TDataType residual_dof_value = 0.0;
        const auto it_dof_begin = rDofSet.begin();
        const int number_of_dof = static_cast<int>(rDofSet.size());

        // Loop over Dofs
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:dof_num)
            for (int i = 0; i < number_of_dof; ++i) {
                auto it_dof = it_dof_begin + i;

                const IndexType dof_id = it_dof->EquationId();

                if (mActiveDofs[dof_id] == 1) {
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                    const auto& r_current_variable = it_dof->GetVariable();
                    const IndexType var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                    #pragma omp atomic
                    residual_solution_norm[var_local_key] += std::pow(residual_dof_value, 2);
                    ++dof_num;
                }
            }
        } else {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:dof_num)
            for (int i = 0; i < number_of_dof; ++i) {
                auto it_dof = it_dof_begin + i;

                if (!it_dof->IsFixed()) {
                    const IndexType dof_id = it_dof->EquationId();
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                    const auto& r_current_variable = it_dof->GetVariable();
                    const IndexType var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                    #pragma omp atomic
                    residual_solution_norm[var_local_key] += std::pow(residual_dof_value, 2);
                    ++dof_num;
                }
            }
        }

        rDofNum = dof_num;
        #pragma omp parallel for
        for (int i = 0; i < number_of_dof; ++i) {
            rResidualSolutionNorm[i] = std::sqrt(residual_solution_norm[i]);
        }
    }

    /**
     * @brief Calculate the convergence norms
     * This method calculates the convergence norms for all the variables to be checked
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual)
     * @return std::tuple<std::vector<TDataType>, std::vector<TDataType>> Tuple containing the absolute and relative convergence values
     */
    std::tuple<std::vector<TDataType>, std::vector<TDataType>> CalculateConvergenceNorms(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        std::vector<int> dofs_count(mVariableSize, 0);
        std::vector<TDataType> solution_norms_vector(mVariableSize, 0.0);
        std::vector<TDataType> increase_norms_vector(mVariableSize, 0.0);

        // Accumulate the norm values
        GetNormValues(rModelPart, rDofSet, rb, dofs_count, solution_norms_vector, increase_norms_vector);

        // Synchronize the norm values
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        auto global_solution_norms_vector = r_data_comm.SumAll(solution_norms_vector);
        auto global_increase_norms_vector = r_data_comm.SumAll(increase_norms_vector);
        auto global_dofs_count = r_data_comm.SumAll(dofs_count);

        // Check division by zero in global solution norms
        const double zero_tol = 1.0e-12;
        for(int i = 0; i < mVariableSize; ++i) {
            if (global_solution_norms_vector[i] < zero_tol) {
                global_solution_norms_vector[i] = 1.0;
            }
        }

        // Calculate the norm values
        std::vector<TDataType> var_ratio(mVariableSize, 0.0);
        std::vector<TDataType> var_abs(mVariableSize, 0.0);
        for(int i = 0; i < mVariableSize; ++i) {
            var_ratio[i] = std::sqrt(global_increase_norms_vector[i] / global_solution_norms_vector[i]);
            var_abs[i] = std::sqrt(global_increase_norms_vector[i]) / static_cast<TDataType>(global_dofs_count[i]);
        }

        // Output the ratio and absolute norms as a tuple
        return std::make_tuple(var_ratio, var_abs);
    }

    /**
     * @brief Method to output the convergence status
     * This method prints the convergence status to the screen for each one of the checked variables
     * @param rConvergenceNorms Tuple containing the absolute and relative convergence values
     */
    virtual void OutputConvergenceStatus(const std::tuple<std::vector<TDataType>,std::vector<TDataType>>& rConvergenceNorms)
    {
        ConvergenceCriteriaUtilities::OutputConvergenceStatus(rConvergenceNorms, mVariableSize, mVariableDataVector, mRatioToleranceVector, mAbsToleranceVector, mLocalKeyMap, this->GetEchoLevel());
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Private Static Member Variables
    ///@{


    ///@}
    ///@name Private Member Variables
    ///@{

    const int mVariableSize;                                    /// The size of the number of variables

    const std::vector<const VariableData*> mVariableDataVector; /// The variables to be checked

    const std::vector<TDataType> mRatioToleranceVector;         /// The ratio threshold for the norm of the residual

    const std::vector<TDataType> mAbsToleranceVector;           /// The absolute value threshold for the norm of the residual

    std::vector<TDataType> mInitialResidualNormVector;          /// The reference norm of the residual

    std::vector<TDataType> mCurrentResidualNormVector;          /// The current norm of the residual

    std::vector<TDataType> mReferenceDispNormVector;            /// The norm at the beginning of the iterations

    std::unordered_map<IndexType, IndexType> mLocalKeyMap;          /// The map containing the local keys

    std::vector<int> mActiveDofs;                               /// This vector contains the dofs that are active

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Get the Norm Values
     * This function accumulates the solution and increment norm values in the provided arrays.
     * Note that these arrays are assumed to be already initialized to zero.
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual)
     * @param rDofsCount Array containing the number of DOFs per variable
     * @param rSolutionNormsVector Array containing the solution norms accumulated values for each variable checked
     * @param rIncreaseNormsVector Array containing the correction norms accumulated values for each variable checked
     */
    virtual void GetNormValues(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rb,
        std::vector<int>& rDofsCount,
        std::vector<TDataType>& rSolutionNormsVector,
        std::vector<TDataType>& rIncreaseNormsVector)
    {
        const int n_dofs = rDofSet.size();

        // Loop over Dofs
        #pragma omp parallel
        {
            // Local thread variables
            int dof_id;
            TDataType dof_rhs;
            TDataType dof_value;

            // Local reduction variables
            std::vector<TDataType> var_solution_norm_reduction(mVariableSize);
            std::vector<TDataType> var_correction_norm_reduction(mVariableSize);
            std::vector<int> dofs_counter_reduction(mVariableSize);
            for (int i = 0; i < mVariableSize; ++i) {
                var_solution_norm_reduction[i] = 0.0;
                var_correction_norm_reduction[i] = 0.0;
                dofs_counter_reduction[i] = 0;
            }

            const auto it_dof_begin = rDofSet.begin();

            #pragma omp for
            for (int i = 0; i < n_dofs; ++i) {
                auto it_dof = it_dof_begin + i;
                if (it_dof->IsFree()) {
                    dof_id = it_dof->EquationId();
                    dof_value = it_dof->GetSolutionStepValue(0);
                    dof_rhs = TSparseSpace::GetValue(rb, dof_id);

                    const auto& r_current_variable = it_dof->GetVariable();
                    const IndexType var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                    var_solution_norm_reduction[var_local_key] += dof_value * dof_value;
                    var_correction_norm_reduction[var_local_key] += dof_rhs * dof_rhs;
                    dofs_counter_reduction[var_local_key]++;
                }
            }

            #pragma omp critical
            {
                for (int i = 0; i < mVariableSize; ++i) {
                    rDofsCount[i] += dofs_counter_reduction[i];
                    rSolutionNormsVector[i] += var_solution_norm_reduction[i];
                    rIncreaseNormsVector[i] += var_correction_norm_reduction[i];
                }
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
};
///@} // Kratos classes

///@} // Application group
}

#endif // KRATOS_MIXED_GENERIC_RESIDUAL_CRITERIA_H
