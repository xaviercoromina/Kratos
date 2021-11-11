//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Riccardo Rossi
//                  Raul Bravo
//
#if !defined(KRATOS_ROM_BUILDER_AND_SOLVER)
#define KRATOS_ROM_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/qr_utility.h"
#include "utilities/atomic_utilities.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ROMBuilderAndSolver : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise builder and solver
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the Element_Variables and Condition_Variables
     * which will be asked and set by another strategy object
     */

    //pointer definition

    KRATOS_CLASS_POINTER_DEFINITION(ROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ROMBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters)
        : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns": [],
            "number_of_rom_dofs": 0,
            "build_petrov_galerkin": false,
            "solve_petrov_galerkin": false,
            "solve_least_squares" : false,
            "nodal_residual_unknowns": [],
            "number_of_rom_residual_dofs": 0
        })");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // We set the other member variables
        mpLinearSystemSolver = pNewLinearSystemSolver;

        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();

        mNodalDofs = mNodalVariablesNames.size();
        mRomDofs = ThisParameters["number_of_rom_dofs"].GetInt();
        if (ThisParameters["build_petrov_galerkin"].GetBool() == true){
            mBuildPetrovGalerkin = true;
        }
        else if (ThisParameters["solve_least_squares"].GetBool() == true){
            mSolveLeastSquares = true;
            mRomDofs_petrov = mRomDofs;
            }
        else if (ThisParameters["solve_petrov_galerkin"].GetBool() == true){
            mSolvePetrovGalerkin = true;
            mRomDofs_petrov = ThisParameters["number_of_rom_residual_dofs"].GetInt();
            }
        else {
            mRomDofs_petrov = mRomDofs;
        }

        // Setting up mapping: VARIABLE_KEY --> CORRECT_ROW_IN_BASIS
        for(int k=0; k<mNodalDofs; k++){
            if(KratosComponents<Variable<double>>::Has(mNodalVariablesNames[k]))
            {
                const auto& var = KratosComponents<Variable<double>>::Get(mNodalVariablesNames[k]);
                mMapPhi[var.Key()] = k;
            }
            else
                KRATOS_ERROR << "variable \""<< mNodalVariablesNames[k] << "\" not valid" << std::endl;

        }
    }

    /** Destructor.
     */
    ~ROMBuilderAndSolver() = default;

    virtual void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        //Gets the array of elements from the modeler
        auto &r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

        typedef std::unordered_set<NodeType::DofType::Pointer, DofPointerHasher> set_type;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set;
        dof_global_set.reserve(number_of_elements * 20);

        double number_of_hrom_elements=0.0;
        #pragma omp parallel firstprivate(dof_list, second_dof_list) reduction(+:number_of_hrom_elements)
        {
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We cleate the temporal set and we reserve some space on them
            set_type dofs_tmp_set;
            dofs_tmp_set.reserve(20000);
            // Gets the array of elements from the modeler
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i)
            {
                auto it_elem = r_elements_array.begin() + i;
                //detect whether the element has a Hyperreduced Weight (H-ROM simulation) or not (ROM simulation)
                if ((it_elem)->Has(HROM_WEIGHT))
                    number_of_hrom_elements++;
                else
                    it_elem->SetValue(HROM_WEIGHT, 1.0);
                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of conditions from the modeler
            ConditionsArrayType &r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());

            ModelPart::ConditionsContainerType selected_conditions_private;
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_conditions; ++i)
            {
                auto it_cond = r_conditions_array.begin() + i;
                // Gather the H-reduced conditions that are to be considered for assembling. Ignoring those for displaying results only
                if (it_cond->Has(HROM_WEIGHT)){
                    selected_conditions_private.push_back(*it_cond.base());
                    number_of_hrom_elements++;
                }
                else
                    it_cond->SetValue(HROM_WEIGHT, 1.0);
                // Gets list of Dof involved on every elem*pcurrent_rom_nodal_basis = nullptrent
                pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }
            #pragma omp critical
            {
                for (auto &cond : selected_conditions_private){
                    mSelectedConditions.push_back(&cond);
                }
            }

            // Gets the array of constraints from the modeler
            auto &r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_constraints; ++i)
            {
                auto it_const = r_constraints_array.begin() + i;

                // Gets list of Dof involved on every element
                it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }
        if (number_of_hrom_elements>0){
             mHromSimulation = true;
        }

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dof_global_set.size());
        for (auto it = dof_global_set.begin(); it != dof_global_set.end(); it++)
        {
            Doftemp.push_back(*it);
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "End of setup dof set\n"
                                                                          << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                                                                 << "Node : " << dof_iterator->Id() << std::endl
                                                                 << "Dof : " << (*dof_iterator) << std::endl
                                                                 << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        KRATOS_CATCH("");
    }

    /**
            organises the dofset in order to speed up the building phase
     */
    virtual void SetUpSystem(
        ModelPart &r_model_part
    ) override
    {
        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();
        int ndofs = static_cast<int>(BaseType::mDofSet.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < static_cast<int>(ndofs); i++){
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + i;
            dof_iterator->SetEquationId(i);
        }
    }


    // Vector ProjectToReducedBasis(
	// 	const TSystemVectorType& rX,
	// 	ModelPart::NodesContainerType& rNodes
	// )
    // {
    //     Vector rom_unknowns = ZeroVector(mRomDofs);
    //     for(const auto& node : rNodes)
    //     {
    //         unsigned int node_aux_id = node.GetValue(AUX_ID);
    //         const auto& nodal_rom_basis = node.GetValue(ROM_BASIS);
	// 			for (int i = 0; i < mRomDofs; ++i) {
	// 				for (int j = 0; j < mNodalDofs; ++j) {
	// 					rom_unknowns[i] += nodal_rom_basis(j, i)*rX(node_aux_id*mNodalDofs + j);
	// 				}
	// 			}
    //     }
    //     return rom_unknowns;
	// }

    void ProjectToFineBasis(
        const TSystemVectorType &rRomUnkowns,
        ModelPart &rModelPart,
        TSystemVectorType &Dx)
    {
        const auto dofs_begin = BaseType::mDofSet.begin();
        const auto dofs_number = BaseType::mDofSet.size();

        #pragma omp parallel firstprivate(dofs_begin, dofs_number)
        {
            const Matrix *pcurrent_rom_nodal_basis = nullptr;
            unsigned int old_dof_id;
            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(dofs_number); k++){
                auto dof = dofs_begin + k;
                if(pcurrent_rom_nodal_basis == nullptr){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                else if(dof->Id() != old_dof_id ){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                Dx[dof->EquationId()] = inner_prod(  row(  *pcurrent_rom_nodal_basis    , mMapPhi[dof->GetVariable().Key()]   )     , rRomUnkowns);
            }
        }
    }

    void GetPhiElemental(
        Matrix &PhiElemental,
        const Element::DofsVectorType &dofs,
        const Element::GeometryType &geom)
    {
        const Matrix *pcurrent_rom_nodal_basis = nullptr;
        int counter = 0;
        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
            auto variable_key = dofs[k]->GetVariable().Key();
            if(k==0)
                pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
            else if(dofs[k]->Id() != dofs[k-1]->Id()){
                counter++;
                pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
            }
            if (dofs[k]->IsFixed())
                noalias(row(PhiElemental, k)) = ZeroVector(PhiElemental.size2());
            else
                noalias(row(PhiElemental, k)) = row(*pcurrent_rom_nodal_basis, mMapPhi[variable_key]);
        }
    }


    void GetPhiElementalResidual(
        Matrix &PhiElementalResidual,
        const Element::DofsVectorType &dofs,
        const Element::GeometryType &geom)
    {
        const Matrix *pcurrent_rom_nodal_basis = nullptr;
        int counter = 0;
        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
            auto variable_key = dofs[k]->GetVariable().Key();
            if(k==0)
                pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS_ASSEMBLED_RESIDUALS));
            else if(dofs[k]->Id() != dofs[k-1]->Id()){
                counter++;
                pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS_ASSEMBLED_RESIDUALS));
            }
            if (dofs[k]->IsFixed())
                noalias(row(PhiElementalResidual, k)) = ZeroVector(PhiElementalResidual.size2());
            else
                noalias(row(PhiElementalResidual, k)) = row(*pcurrent_rom_nodal_basis, mMapPhi[variable_key]);
        }
    }

    void GetPhiGlobal(Matrix &PhiGlobal,
    ModelPart &rModelPart)
    {
        const auto dofs_begin = BaseType::mDofSet.begin();
        const auto dofs_number = BaseType::mDofSet.size();

        #pragma omp parallel firstprivate(dofs_begin, dofs_number)
        {
            const Matrix *pcurrent_rom_nodal_basis = nullptr;
            unsigned int old_dof_id;
            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(dofs_number); k++){
                auto dof = dofs_begin + k;
                if(pcurrent_rom_nodal_basis == nullptr){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                else if(dof->Id() != old_dof_id ){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                noalias(row(PhiGlobal,dof->EquationId())) = row(*pcurrent_rom_nodal_basis,mMapPhi[dof->GetVariable().Key()]);
            }
        }
    }

    void GetPsiGlobal(Matrix &PsiGlobal,
    ModelPart &rModelPart)
    {
        const auto dofs_begin = BaseType::mDofSet.begin();
        const auto dofs_number = BaseType::mDofSet.size();

        #pragma omp parallel firstprivate(dofs_begin, dofs_number)
        {
            const Matrix *pcurrent_rom_nodal_basis = nullptr;
            unsigned int old_dof_id;
            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(dofs_number); k++){
                auto dof = dofs_begin + k;
                if(pcurrent_rom_nodal_basis == nullptr){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS_ASSEMBLED_RESIDUALS));
                    old_dof_id = dof->Id();
                }
                else if(dof->Id() != old_dof_id ){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS_ASSEMBLED_RESIDUALS));
                    old_dof_id = dof->Id();
                }
                noalias(row(PsiGlobal,dof->EquationId())) = row(*pcurrent_rom_nodal_basis,mMapPhi[dof->GetVariable().Key()]);
            }
        }
    }

//     inline void AssembleRowContribution(TSystemMatrixType& A, const Matrix& Alocal, const unsigned int i, const unsigned int i_local, Element::EquationIdVectorType& EquationId)
//     {
//         double* values_vector = A.value_data().begin();
//         std::size_t* index1_vector = A.index1_data().begin();
//         std::size_t* index2_vector = A.index2_data().begin();

//         size_t left_limit = index1_vector[i];
// //    size_t right_limit = index1_vector[i+1];

//         //find the first entry
//         size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
//         size_t last_found = EquationId[0];

//         double& r_a = values_vector[last_pos];
//         const double& v_a = Alocal(i_local,0);
//         AtomicAdd(r_a,  v_a);

//         //now find all of the other entries
//         size_t pos = 0;
//         for (unsigned int j=1; j<EquationId.size(); j++) {
//             unsigned int id_to_find = EquationId[j];
//             if(id_to_find > last_found) {
//                 pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
//             } else if(id_to_find < last_found) {
//                 pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
//             } else {
//                 pos = last_pos;
//             }

//             double& r = values_vector[pos];
//             const double& v = Alocal(i_local,j);
//             AtomicAdd(r,  v);

//             last_found = id_to_find;
//             last_pos = pos;
//         }
//     }

//     inline unsigned int ForwardFind(const unsigned int id_to_find,
//                                     const unsigned int start,
//                                     const size_t* index_vector)
//     {
//         unsigned int pos = start;
//         while(id_to_find != index_vector[pos]) pos++;
//         return pos;
//     }

//     inline unsigned int BackwardFind(const unsigned int id_to_find,
//                                      const unsigned int start,
//                                      const size_t* index_vector)
//     {
//         unsigned int pos = start;
//         while(id_to_find != index_vector[pos]) pos--;
//         return pos;
//     }

    // void Assemble(
    //     TSystemMatrixType& A,
    //     TSystemVectorType& b,
    //     const LocalSystemMatrixType& LHS_Contribution,
    //     const LocalSystemVectorType& RHS_Contribution,
    //     Element::EquationIdVectorType& EquationId
    // )
    // {
    //     unsigned int local_size = LHS_Contribution.size1();

    //     for (unsigned int i_local = 0; i_local < local_size; i_local++) {
    //         unsigned int i_global = EquationId[i_local];

    //         double& r_a = b[i_global];
    //         const double& v_a = RHS_Contribution(i_local);
    //         AtomicAdd(r_a, v_a);

    //         AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);
    //     }
    // }




    /*@{ */

    /**
            Function to perform the building and solving phase at the same time.
            It is ideally the fastest and safer function to use when it is possible to solve
            just after building
     */
    virtual void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        //define a dense matrix to hold the reduced problem
        Matrix Arom;
        Vector brom;
        if (mBuildPetrovGalerkin== true){
            Arom = ZeroMatrix(Dx.size(), mRomDofs);
            brom = ZeroVector(Dx.size());
        }
        else {
            Arom = ZeroMatrix(mRomDofs_petrov, mRomDofs);
            brom = ZeroVector(mRomDofs_petrov);
        }
        
        TSystemVectorType x(Dx.size());

        double project_to_reduced_start = OpenMPUtils::GetCurrentTime();
        Vector xrom = ZeroVector(mRomDofs);
        //this->ProjectToReducedBasis(x, rModelPart.Nodes(),xrom);
        const double project_to_reduced_end = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to reduced basis time: " << project_to_reduced_end - project_to_reduced_start << std::endl;

        //build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());
        const auto el_begin = rModelPart.ElementsBegin();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        auto help_cond_begin = rModelPart.ConditionsBegin();
        auto help_nconditions = static_cast<int>(rModelPart.Conditions().size());

        if ( mHromSimulation == true){
            // Only selected conditions are considered for the calculation on an H-ROM simualtion.
            help_cond_begin = mSelectedConditions.begin();
            help_nconditions = static_cast<int>(mSelectedConditions.size());
        }

        // Getting the array of the conditions
        const auto cond_begin = help_cond_begin;
        const auto nconditions = help_nconditions;


        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();
        if (mBuildPetrovGalerkin== true){
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix tempA;
                Vector tempb;
                tempA = ZeroMatrix(Dx.size(),mRomDofs);
                tempb = ZeroVector(Dx.size());
                Matrix aux;
                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nelements); k++)
                {
                    auto it_el = el_begin + k;
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    
                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
                            if(dofs[k]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[k]->EquationId()))+=row(aux,k);
                            tempb[dofs[k]->EquationId()]+=RHS_Contribution(k);
                        }
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nconditions); k++)
                {
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);

                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
                            if(dofs[k]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[k]->EquationId()))+=row(aux,k);
                            tempb[dofs[k]->EquationId()]+=RHS_Contribution(k);
                        }
                    }
                }

                #pragma omp critical
                {
                    noalias(Arom) +=tempA;
                    noalias(brom) +=tempb;
                }
            } 
            
        }
        else if (mSolveLeastSquares== true){
            Matrix globalAphi;
            Vector globalb;
            globalAphi = ZeroMatrix(Dx.size(),mRomDofs);
            globalb = ZeroVector(Dx.size());
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix tempA;
                Vector tempb;
                tempA = ZeroMatrix(Dx.size(),mRomDofs);
                tempb = ZeroVector(Dx.size());
                Matrix aux;
                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nelements); k++)
                {
                    auto it_el = el_begin + k;
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    
                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
                            if(dofs[k]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[k]->EquationId()))+=row(aux,k);
                            tempb[dofs[k]->EquationId()]+=RHS_Contribution(k);
                        }
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nconditions); k++)
                {
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);

                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
                            if(dofs[k]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[k]->EquationId()))+=row(aux,k);
                            tempb[dofs[k]->EquationId()]+=RHS_Contribution(k);
                        }
                    }
                }

                #pragma omp critical
                {
                    noalias(globalAphi) +=tempA;
                    noalias(globalb) +=tempb;
                }
            }
            noalias(Arom) = prod(trans(globalAphi),globalAphi);
            noalias(brom) = prod(trans(globalAphi),globalb);
        }
        else if (mBuildPetrovGalerking== true){
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix tempA;
                Vector tempb;
                tempA = ZeroMatrix(mRomDofs,mRomDofs);
                tempb = ZeroVector(mRomDofs);
                
                Matrix aux;
                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nelements); k++)
                {
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element
                    //is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    
                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        // for(int l = 0; l < static_cast<int>(dofs.size()); ++l){
                        //     if(dofs[l]->IsFixed()){  //When dof is fixed set to zero the corresponging elemental LHS row (==not adding contribution)
                        //         noalias(row(aux,l)) = ZeroVector(aux.size2());
                        //     }
                        // }
                        double h_rom_weight = it_el->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(aux),aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(aux),RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nconditions); k++){
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);

                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int l = 0; l < static_cast<int>(dofs.size()); ++l){
                            if(dofs[l]->IsFixed()){  //When dof is fixed set to zero the corresponging elemental LHS row (==not adding contribution)
                                noalias(row(aux,l)) = ZeroVector(aux.size2());
                            }
                        }
                        double h_rom_weight = it->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(aux),aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(aux),RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp critical
                {
                    noalias(Arom) +=tempA;
                    noalias(brom) +=tempb;
                }

            }
        }
        else if (mSolvePetrovGalerking == true){
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix PhiElementalResidual;
                Matrix tempA;
                Vector tempb;
                tempA = ZeroMatrix(mRomDofs_petrov,mRomDofs);
                tempb = ZeroVector(mRomDofs_petrov);
            
                Matrix aux;

                #pragma omp for nowait
                for (int k = 0; k < nelements; k++)
                {
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element
                    //is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);

                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        if(PhiElementalResidual.size1() != dofs.size() || PhiElementalResidual.size2() != mRomDofs_petrov)
                            PhiElementalResidual.resize(dofs.size(), mRomDofs_petrov,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        GetPhiElementalResidual(PhiElementalResidual,dofs,geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        double h_rom_weight = it_el->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(PhiElementalResidual), aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(PhiElementalResidual), RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < nconditions; k++){
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        if(PhiElementalResidual.size1() != dofs.size() || PhiElementalResidual.size2() != mRomDofs_petrov)
                            PhiElementalResidual.resize(dofs.size(), mRomDofs_petrov,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        GetPhiElementalResidual(PhiElementalResidual,dofs,geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        double h_rom_weight = it->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(PhiElementalResidual), aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(PhiElementalResidual), RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp critical
                {
                    noalias(Arom) +=tempA;
                    noalias(brom) +=tempb;
                }
            }
        }
        else if (mSolvePetrovGalerkin == true){
            Matrix globalAphi;
            Vector globalb;
            globalAphi = ZeroMatrix(Dx.size(),mRomDofs);
            globalb = ZeroVector(Dx.size());
            // #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix tempA;
                Vector tempb;
                tempA = ZeroMatrix(Dx.size(),mRomDofs);
                tempb = ZeroVector(Dx.size());
                Matrix aux;
                // #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nelements); k++)
                {
                    auto it_el = el_begin + k;
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    
                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int l = 0; l < static_cast<int>(dofs.size()); l++){
                            if(dofs[l]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[l]->EquationId()))+=row(aux,l);
                            tempb[dofs[l]->EquationId()]+=RHS_Contribution(l);
                        }
                    }
                }

                // #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(nconditions); k++)
                {
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);

                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        for(int l = 0; l < static_cast<int>(dofs.size()); l++){
                            if(dofs[l]->IsFixed()==false)  //When dof is fixed set to zero the corresponging LHS row (==not adding contribution)
                                noalias(row(tempA,dofs[l]->EquationId()))+=row(aux,l);
                            tempb[dofs[l]->EquationId()]+=RHS_Contribution(l);
                        }
                    }
                }

                // #pragma omp critical
                {
                    noalias(globalAphi) +=tempA;
                    noalias(globalb) +=tempb;
                }
            }
            Matrix globalPsi;
            globalPsi = ZeroMatrix(Dx.size(),mRomDofs_petrov);
            GetPsiGlobal(globalPsi,rModelPart);
            noalias(Arom) = prod(trans(globalPsi),globalAphi);
            noalias(brom) = prod(trans(globalPsi),globalb);
        }
        else {
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                Matrix PhiElemental;
                Matrix tempA = ZeroMatrix(mRomDofs,mRomDofs);
                Vector tempb = ZeroVector(mRomDofs);
                Matrix aux;

                #pragma omp for nowait
                for (int k = 0; k < nelements; k++)
                {
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element
                    //is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE)){
                        element_is_active = (it_el)->Is(ACTIVE);
                    }


                    //KRATOS_WATCH("entered elements loop")

                    if (element_is_active){
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        const auto &geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        double h_rom_weight = it_el->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < nconditions; k++){
                    auto it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the condition
                    //is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE)){
                        condition_is_active = (it)->Is(ACTIVE);
                    }

                    //KRATOS_WATCH("entered conditions loop")

                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        const auto &geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                            aux.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(aux) = prod(LHS_Contribution, PhiElemental);
                        double h_rom_weight = it->GetValue(HROM_WEIGHT);
                        noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                        noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                    }
                }

                #pragma omp critical
                {
                    noalias(Arom) +=tempA;
                    noalias(brom) +=tempb;
                }

            }
        }


        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;


        //solve for the rom unkowns dunk = Arom^-1 * brom
        Vector dxrom(xrom.size());
        double start_solve = OpenMPUtils::GetCurrentTime();
        if (mBuildPetrovGalerkin==true){
            QR<double, row_major> qr_util;
            qr_util.compute(Dx.size(), mRomDofs, &(Arom)(0,0));
            qr_util.solve(&(brom)(0), &(dxrom)(0));
        }
        else if (mSolvePetrovGalerkin == true){
            //Resolver con QR
            QR<double, row_major> qr_util;
            qr_util.compute(mRomDofs_petrov, mRomDofs, &(Arom)(0,0));
            qr_util.solve(&(brom)(0), &(dxrom)(0));
        }
        else if (mSolveLeastSquares == true){
            MathUtils<double>::Solve(Arom, dxrom, brom);
        }
        else {
            MathUtils<double>::Solve(Arom, dxrom, brom);
        }
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Solve reduced system time: " << stop_solve - start_solve << std::endl;

        // //update database
        // noalias(xrom) += dxrom;

        // project reduced solution back to full order model
        double project_to_fine_start = OpenMPUtils::GetCurrentTime();
        ProjectToFineBasis(dxrom, rModelPart, Dx);
        const double project_to_fine_end = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to fine basis time: " << project_to_fine_end - project_to_fine_start << std::endl;
    }

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType &pA,
        TSystemVectorPointerType &pDx,
        TSystemVectorPointerType &pb,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }

        // TSystemMatrixType& A = *pA; //To solve least squares problem A@Phi@Dq=b
        TSystemVectorType &Dx = *pDx;
        TSystemVectorType &b = *pb;

        //resizing the system vectors and matrix
        // if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        // {
        //     A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
        //     ConstructMatrixStructure(pScheme, A, rModelPart);
        // }
        // else
        // {
        //     if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
        //     {
        //         KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
        //         A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
        //         ConstructMatrixStructure(pScheme, A, rModelPart);
        //     }
        // }

        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        // TSparseSpace::SetToZero(Dx);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);
        // TSparseSpace::SetToZero(b);

        KRATOS_CATCH("")
    }

    // virtual void ConstructMatrixStructure(
    //     typename TSchemeType::Pointer pScheme,
    //     TSystemMatrixType& A,
    //     ModelPart& rModelPart)
    // {
    //     //filling with zero the matrix (creating the structure)
    //     Timer::Start("MatrixStructure");

    //     const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

    //     const std::size_t equation_size = BaseType::mEquationSystemSize;

    //     std::vector< LockObject > lock_array(equation_size);

    //     std::vector<std::unordered_set<std::size_t> > indices(equation_size);

    //     block_for_each(indices, [](std::unordered_set<std::size_t>& rIndices){
    //         rIndices.reserve(40);
    //     });

    //     Element::EquationIdVectorType ids(3, 0);

    //     block_for_each(rModelPart.Elements(), ids, [&](Element& rElem, Element::EquationIdVectorType& rIdsTLS){
    //         pScheme->EquationId(rElem, rIdsTLS, CurrentProcessInfo);
    //         for (std::size_t i = 0; i < rIdsTLS.size(); i++) {
    //             lock_array[rIdsTLS[i]].lock();
    //             auto& row_indices = indices[rIdsTLS[i]];
    //             row_indices.insert(rIdsTLS.begin(), rIdsTLS.end());
    //             lock_array[rIdsTLS[i]].unlock();
    //         }
    //     });

    //     block_for_each(rModelPart.Conditions(), ids, [&](Condition& rCond, Element::EquationIdVectorType& rIdsTLS){
    //         pScheme->EquationId(rCond, rIdsTLS, CurrentProcessInfo);
    //         for (std::size_t i = 0; i < rIdsTLS.size(); i++) {
    //             lock_array[rIdsTLS[i]].lock();
    //             auto& row_indices = indices[rIdsTLS[i]];
    //             row_indices.insert(rIdsTLS.begin(), rIdsTLS.end());
    //             lock_array[rIdsTLS[i]].unlock();
    //         }
    //     });

    //     if (rModelPart.MasterSlaveConstraints().size() != 0) {
    //         struct TLS
    //         {
    //             Element::EquationIdVectorType master_ids = Element::EquationIdVectorType(3,0);
    //             Element::EquationIdVectorType slave_ids = Element::EquationIdVectorType(3,0);
    //         };
    //         TLS tls;

    //         block_for_each(rModelPart.MasterSlaveConstraints(), tls, [&](MasterSlaveConstraint& rConst, TLS& rTls){
    //             rConst.EquationIdVector(rTls.slave_ids, rTls.master_ids, CurrentProcessInfo);

    //             for (std::size_t i = 0; i < rTls.slave_ids.size(); i++) {
    //                 lock_array[rTls.slave_ids[i]].lock();
    //                 auto& row_indices = indices[rTls.slave_ids[i]];
    //                 row_indices.insert(rTls.slave_ids[i]);
    //                 lock_array[rTls.slave_ids[i]].unlock();
    //             }

    //             for (std::size_t i = 0; i < rTls.master_ids.size(); i++) {
    //                 lock_array[rTls.master_ids[i]].lock();
    //                 auto& row_indices = indices[rTls.master_ids[i]];
    //                 row_indices.insert(rTls.master_ids[i]);
    //                 lock_array[rTls.master_ids[i]].unlock();
    //             }
    //         });

    //     }

    //     //destroy locks
    //     lock_array = std::vector< LockObject >();

    //     //count the row sizes
    //     unsigned int nnz = 0;
    //     for (unsigned int i = 0; i < indices.size(); i++) {
    //         nnz += indices[i].size();
    //     }

    //     A = CompressedMatrixType(indices.size(), indices.size(), nnz);

    //     double* Avalues = A.value_data().begin();
    //     std::size_t* Arow_indices = A.index1_data().begin();
    //     std::size_t* Acol_indices = A.index2_data().begin();

    //     //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
    //     Arow_indices[0] = 0;
    //     for (int i = 0; i < static_cast<int>(A.size1()); i++) {
    //         Arow_indices[i+1] = Arow_indices[i] + indices[i].size();
    //     }

    //     IndexPartition<std::size_t>(A.size1()).for_each([&](std::size_t i){
    //         const unsigned int row_begin = Arow_indices[i];
    //         const unsigned int row_end = Arow_indices[i+1];
    //         unsigned int k = row_begin;
    //         for (auto it = indices[i].begin(); it != indices[i].end(); it++) {
    //             Acol_indices[k] = *it;
    //             Avalues[k] = 0.0;
    //             k++;
    //         }

    //         indices[i].clear(); //deallocating the memory

    //         std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

    //     });

    //     A.set_filled(indices.size()+1, nnz);

    //     Timer::Stop("MatrixStructure");
    // }

    /*@} */
    /**@name Operations */
    /*@{ */

    /*@} */
    /**@name Access */
    /*@{ */

    /*@} */
    /**@name Inquiry */
    /*@{ */

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /*@} */
    /**@name Friends */
    /*@{ */

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /** Pointer to the Model.
     */
    typename TLinearSolver::Pointer mpLinearSystemSolver;

    //DofsArrayType mDofSet;
    std::vector<DofPointerType> mDofList;

    bool mReshapeMatrixFlag = false;

    /// flag taking care if the dof set was initialized ot not
    bool mDofSetIsInitialized = false;

    /// flag taking in account if it is needed or not to calculate the reactions
    bool mCalculateReactionsFlag = false;

    /// number of degrees of freedom of the problem to be solve
    unsigned int mEquationSystemSize;
    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    int mEchoLevel = 0;

    TSystemVectorPointerType mpReactionsVector;

    std::vector<std::string> mNodalVariablesNames;
    int mNodalDofs;
    unsigned int mRomDofs;
    unsigned int mRomDofs_petrov;
    std::unordered_map<Kratos::VariableData::KeyType,int> mMapPhi;
    ModelPart::ConditionsContainerType mSelectedConditions;
    bool mSolvePetrovGalerkin = false,mBuildPetrovGalerkin = false,mBuildPetrovGalerking = false, mSolveLeastSquares = false, mSolvePetrovGalerking = false;
    bool mHromSimulation = false;

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /*@} */
    /**@name Private  Access */
    /*@{ */

    /*@} */
    /**@name Private Inquiry */
    /*@{ */

    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class ROMBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_ROM_BUILDER_AND_SOLVER  defined */
