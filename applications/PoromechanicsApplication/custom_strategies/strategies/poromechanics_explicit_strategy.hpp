
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_POROMECHANICS_EXPLICIT_STRATEGY)
#define KRATOS_POROMECHANICS_EXPLICIT_STRATEGY

/* System includes */
// #include <fstream>

// Project includes
#include "custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp"
#include "utilities/parallel_utilities.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class PoromechanicsExplicitStrategy
    : public MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsExplicitStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// DoF types definition
    typedef typename Node<3>::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    using MotherType::mInitializeWasPerformed;
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mpScheme;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsExplicitStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        Parameters& rParameters,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            Parameters default_parameters( R"(
                {
                    "rebuild_level": 0,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "initial_radius": 1.0e-12,
                    "characteristic_length": 0.05,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list" : []
                }  )" );

            // Validate agains defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);

            mpParameters = &rParameters;

            // Set Load SubModelParts and Variable names
            if(rParameters["loads_sub_model_part_list"].size() > 0)
            {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                if( mSubModelPartList.size() != mVariableNames.size() )
                    KRATOS_THROW_ERROR( std::logic_error, "For each SubModelPart there must be a corresponding nodal Variable", "" )

                for(unsigned int i = 0; i < mVariableNames.size(); i++)
                {
                    mSubModelPartList[i] = &( model_part.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }

            BaseType::SetRebuildLevel(rParameters["rebuild_level"].GetInt());

            mNumberOfStepsToConverge = 0;
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoromechanicsExplicitStrategy() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false) {

            ModelPart& r_model_part = BaseType::GetModelPart();

            TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (!mpScheme->SchemeIsInitialized())mpScheme->Initialize(r_model_part);

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (!mpScheme->ElementsAreInitialized())mpScheme->InitializeElements(r_model_part);

            // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
            if (!mpScheme->ConditionsAreInitialized())mpScheme->InitializeConditions(r_model_part);

            // Set Nodal Mass to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
            const array_1d<double, 3> nodal_mass_array = ZeroVector(3);
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS_ARRAY, nodal_mass_array, r_nodes);
            // TODO: Set Nodal AntiCompressibility to zero for mass-balance equation (C=1/Q, with Q being the compressibility coeff.)

            // Iterate over the elements
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

            Vector dummy_vector;
            #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                // Getting nodal mass and inertia from element
                // this function needs to be implemented in the respective
                // element to provide nodal masses
                auto it_elem = it_elem_begin + i;
                it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initial operations ... things that are constant over the Solution Step
        mpScheme->InitializeSolutionStep(r_model_part, matrix_a_dummy, rDx, rb);

        if (BaseType::mRebuildLevel > 0) {
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();

            // Set Nodal Mass and Damping to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
            const array_1d<double, 3> nodal_mass_array = ZeroVector(3);
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS_ARRAY, nodal_mass_array, r_nodes);

            Vector dummy_vector;
            #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                // Getting nodal mass and inertia from element
                // this function needs to be implemented in the respective
                // element to provide nodal masses
                auto it_elem = it_elem_begin + i;
                it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
            }
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();

        // Some dummy sets and matrices
        DofsArrayType dof_set_dummy;
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initialize the non linear iteration
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        mpScheme->Predict(r_model_part, dof_set_dummy, rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            BaseType::MoveMesh();

        // Explicitly integrates the equation of motion.
        mpScheme->Update(r_model_part, dof_set_dummy, rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            BaseType::MoveMesh();

        // Finalize the non linear iteration
        mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Calculate reactions if required
        if (mCalculateReactionsFlag) {
            this->CalculateReactions(mpScheme, r_model_part);
        }

        // CONVERGENCE CHECK
        return this->CheckConvergence(r_model_part);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();
        // Finalisation of the solution step,
        // operations to be done after achieving convergence, for example the
        // Final Residual Vector (rb) has to be saved in there
        // to avoid error accumulation
        mpScheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);

        // Cleaning memory after the solution
        mpScheme->Clean();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters* mpParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart
    int mNumberOfStepsToConverge;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief This method computes the reactions of the problem
     * @param pScheme The pointer to the integration scheme used
     * @param rModelPart The model part which defines the problem
     * @param rA The LHS of the system (empty)
     * @param rDx The solution of the system (empty)
     * @param rb The RHS of the system (empty)
     */
    virtual void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        // We iterate over the nodes
        auto& r_nodes = rModelPart.Nodes();

        // Auxiliar values
        array_1d<double, 3> force_residual = ZeroVector(3);
        double flux_residual = 0.0;

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int dim = r_current_process_info[DOMAIN_SIZE];

        // Getting
        const auto it_node_begin = r_nodes.begin();

        // Iterating nodes
        #pragma omp parallel for firstprivate(force_residual,flux_residual), schedule(guided,512)
        for(int i=0; i<static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = it_node_begin + i;

            noalias(force_residual) = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            flux_residual = it_node->FastGetSolutionStepValue(FLUX_RESIDUAL);

            if( it_node->IsFixed(DISPLACEMENT_X) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_X);
                r_reaction = -force_residual[0];
            }
            if( it_node->IsFixed(DISPLACEMENT_Y) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Y);
                r_reaction = -force_residual[1];
            }
            if( it_node->IsFixed(WATER_PRESSURE) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_WATER_PRESSURE);
                r_reaction = -flux_residual;
            }
            if(dim==3) {
                if( it_node->IsFixed(DISPLACEMENT_Z) == true ) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Z);
                    r_reaction = -force_residual[2];
                }
            }
        }
    }

    virtual bool CheckConvergence(ModelPart& rModelPart)
    {
        // Initialize variables
        mNumberOfStepsToConverge += 1;

        bool is_converged = false;
        bool is_converged_rx = false;
        bool is_converged_ry = false;
        bool is_converged_rz = false;
        bool is_converged_rwp = false;

        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        const auto it_node_begin = rModelPart.NodesBegin();

        // Calculate maximum reactions
        unsigned int NumThreads = ParallelUtilities::GetNumThreads();
        std::vector<double> maximum_reaction_partition(NumThreads);
        std::vector<double> maximum_reaction_water_pressure_partition(NumThreads);
        
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            maximum_reaction_partition[k] = it_node_begin->FastGetSolutionStepValue(REACTION_X);
            maximum_reaction_water_pressure_partition[k] = it_node_begin->FastGetSolutionStepValue(REACTION_WATER_PRESSURE);

            double reaction_x,reaction_y,reaction_z,reaction_water_pressure;

            #pragma omp for
            for(int i = 0; i < NNodes; i++)
            {
                auto itCurrentNode = it_node_begin + i;

                reaction_x = std::abs(itCurrentNode->FastGetSolutionStepValue(REACTION_X));
                reaction_y = std::abs(itCurrentNode->FastGetSolutionStepValue(REACTION_Y));
                reaction_z = std::abs(itCurrentNode->FastGetSolutionStepValue(REACTION_Z));
                reaction_water_pressure = std::abs(itCurrentNode->FastGetSolutionStepValue(REACTION_WATER_PRESSURE));

                if( reaction_x > maximum_reaction_partition[k] ){
                    maximum_reaction_partition[k] = reaction_x;
                }
                if( reaction_y > maximum_reaction_partition[k] ){
                    maximum_reaction_partition[k] = reaction_y;
                }
                if( reaction_z > maximum_reaction_partition[k] ){
                    maximum_reaction_partition[k] = reaction_z;
                }
                if( reaction_water_pressure > maximum_reaction_water_pressure_partition[k] ){
                    maximum_reaction_water_pressure_partition[k] = reaction_water_pressure;
                }
            }
        }

        double maximum_reaction = maximum_reaction_partition[0];
        double maximum_reaction_water_pressure = maximum_reaction_water_pressure_partition[0];

        for(unsigned int i=1; i < NumThreads; i++)
        {
            if(maximum_reaction_partition[i] > maximum_reaction){
                maximum_reaction = maximum_reaction_partition[i];
            }
            if(maximum_reaction_water_pressure_partition[i] > maximum_reaction_water_pressure){
                maximum_reaction_water_pressure = maximum_reaction_water_pressure_partition[i];
            }
        }

        // Calculate total reactions
        double total_reaction_x = 0.0;
        double total_reaction_y = 0.0;
        double total_reaction_z = 0.0;
        double total_reaction_water_pressure = 0.0;
        #pragma omp parallel for reduction(+:total_reaction_x,total_reaction_y,total_reaction_z,total_reaction_water_pressure)
        for (int i = 0; i < NNodes; ++i) {
            auto itCurrentNode = it_node_begin + i;
            const double& r_reaction_x = itCurrentNode->FastGetSolutionStepValue(REACTION_X);
            const double& r_reaction_y = itCurrentNode->FastGetSolutionStepValue(REACTION_Y);
            const double& r_reaction_z = itCurrentNode->FastGetSolutionStepValue(REACTION_Z);
            const double& r_reaction_water_pressure = itCurrentNode->FastGetSolutionStepValue(REACTION_WATER_PRESSURE);

            total_reaction_x += r_reaction_x;
            total_reaction_y += r_reaction_y;
            total_reaction_z += r_reaction_z;
            total_reaction_water_pressure += r_reaction_water_pressure;
        }

        // Calculate relative total reactions
        double relative_total_reaction_x = 0.0;
        double relative_total_reaction_y = 0.0;
        double relative_total_reaction_z = 0.0;
        double relative_total_reaction_water_pressure = 0.0;
        const double abs_total_reaction_x = std::abs(total_reaction_x);
        const double abs_total_reaction_y = std::abs(total_reaction_y);
        const double abs_total_reaction_z = std::abs(total_reaction_z);
        const double abs_total_reaction_water_pressure = std::abs(total_reaction_water_pressure);
        if(maximum_reaction > 1.0e-12) {
            relative_total_reaction_x = abs_total_reaction_x/maximum_reaction;
            relative_total_reaction_y = abs_total_reaction_y/maximum_reaction;
            relative_total_reaction_z = abs_total_reaction_z/maximum_reaction;
        }
        if(maximum_reaction_water_pressure > 1.0e-12) {
            relative_total_reaction_water_pressure = abs_total_reaction_water_pressure/maximum_reaction_water_pressure;
        }

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        if(relative_total_reaction_x <= r_current_process_info[ERROR_RATIO]){
            is_converged_rx = true;
        }
        if(relative_total_reaction_y <= r_current_process_info[ERROR_RATIO]){
            is_converged_ry = true;
        }
        if(relative_total_reaction_z <= r_current_process_info[ERROR_RATIO]){
            is_converged_rz = true;
        }
        if(relative_total_reaction_water_pressure <= r_current_process_info[ERROR_RATIO]){
            is_converged_rwp = true;
        }
        if(is_converged_rx==true && is_converged_ry==true && is_converged_rz==true && is_converged_rwp==true) {
            is_converged = true;

            KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "Reaction convergence is achieved after: " << mNumberOfStepsToConverge << " steps." << std::endl;
            KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "relative_total_reaction_x: " << relative_total_reaction_x << std::endl;
            KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "relative_total_reaction_y: " << relative_total_reaction_y << std::endl;
            KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "relative_total_reaction_z: " << relative_total_reaction_z << std::endl;
            KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "relative_total_reaction_water_pressure: " << relative_total_reaction_water_pressure << std::endl;
            // KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "total_reaction_x: " << total_reaction_x << std::endl;
            // KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "total_reaction_y: " << total_reaction_y << std::endl;
            // KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "total_reaction_z: " << total_reaction_z << std::endl;
            // KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "total_reaction_water_pressure: " << total_reaction_water_pressure << std::endl;
            mNumberOfStepsToConverge = 0;
        }

        return is_converged;
    }

}; // Class PoromechanicsExplicitStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_EXPLICIT_STRATEGY  defined
