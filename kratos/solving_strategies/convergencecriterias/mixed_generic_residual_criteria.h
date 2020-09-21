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

        std::vector<SizeType> size_residual;
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
        std::vector<SizeType>& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        std::vector<SizeType> dofs_count(mVariableSize, 0);
        std::vector<TDataType> residual_solution_norm = std::vector<TDataType>(mVariableSize, 0.0);

        // Auxiliar values
        TDataType residual_dof_value = 0.0;
        const auto it_dof_begin = rDofSet.begin();
        const int number_of_dof = static_cast<int>(rDofSet.size());

        // Loop over Dofs
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            #pragma omp parallel for firstprivate(residual_dof_value)
            for (int i = 0; i < number_of_dof; ++i) {
                auto it_dof = it_dof_begin + i;

                const IndexType dof_id = it_dof->EquationId();

                if (mActiveDofs[dof_id] == 1) {
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                    const auto& r_current_variable = it_dof->GetVariable();
                    const IndexType var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                    #pragma omp atomic
                    residual_solution_norm[var_local_key] += std::pow(residual_dof_value, 2);
                    #pragma omp atomic
                    dofs_count[var_local_key] += 1;
                }
            }
        } else {
            #pragma omp parallel for firstprivate(residual_dof_value)
            for (int i = 0; i < number_of_dof; ++i) {
                auto it_dof = it_dof_begin + i;

                if (!it_dof->IsFixed()) {
                    const IndexType dof_id = it_dof->EquationId();
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                    const auto& r_current_variable = it_dof->GetVariable();
                    const IndexType var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                    #pragma omp atomic
                    residual_solution_norm[var_local_key] += std::pow(residual_dof_value, 2);
                    #pragma omp atomic
                    dofs_count[var_local_key] += 1;
                }
            }
        }

        rDofNum = dofs_count;
        #pragma omp parallel for
        for (int i = 0; i < mVariableSize; ++i) {
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
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Calculate the norm values
        std::vector<TDataType> var_ratio(mVariableSize, 0.0);
        std::vector<TDataType> var_abs(mVariableSize, 0.0);

        std::vector<SizeType> size_residual(mVariableSize, 0);
        CalculateResidualNorm(rModelPart, mCurrentResidualNormVector, size_residual, rDofSet, rb);

        #pragma omp parallel for
        for (int i = 0; i < mVariableSize; ++i) {
            if(mInitialResidualNormVector[i] < std::numeric_limits<TDataType>::epsilon()) {
                var_ratio[i] = 0.0;
            } else {
                var_ratio[i] = mCurrentResidualNormVector[i]/mInitialResidualNormVector[i];
            }

            const TDataType float_size_residual = static_cast<TDataType>(size_residual[i]);
            var_abs[i] = (mCurrentResidualNormVector[i]/float_size_residual);
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

    std::unordered_map<IndexType, IndexType> mLocalKeyMap;      /// The map containing the local keys

    std::vector<int> mActiveDofs;                               /// This vector contains the dofs that are active

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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
