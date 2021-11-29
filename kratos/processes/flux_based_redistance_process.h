//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Pablo Becker
//

#if !defined(KRATOS_POTENTIAL_FILLING_UTILITY_INCLUDED)
#define KRATOS_POTENTIAL_FILLING_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "elements/distance_calculation_flux_based_element_simplex.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "utilities/variable_utils.h"
#include "spatial_containers/spatial_containers.h" 
#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/merge_variable_lists_utility.h"

// #include "filling_application_variables.h"
// #include "altair_extension_application_variables.h"
// #include "custom_utilities/compute_filled_volume_utility.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and computes a field which resembles the filltime
	 * This field is to be used to compute the optimization problem of gate location
	 * 
	*/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class FluxBasedRedistanceProcess : public Process
{
public:
	///@name Type Definitions
	///@{

	typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
	typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;
	typedef typename Geometry<Node<3>>::Pointer GeometryPointer;
	typedef std::vector<GeometryPointer> GeometryPointerVector;

    typedef Node<3>::Pointer NodePointerType;
    typedef std::vector<NodePointerType> NodePointerTypeVector;
    typedef NodePointerTypeVector::iterator NodePointerIterator;
    typedef std::vector<double>::iterator DistanceIterator;
    //typedef Bucket<3, Node<3>, NodePointerTypeVector, NodePointerType, NodePointerIterator, DistanceIterator> BucketType;
    //typedef Tree<KDTreePartition<BucketType>> KdtreeType; //Kdtree
	typedef Bins<TDim, Node<3>, NodePointerTypeVector, NodePointerType, NodePointerIterator, DistanceIterator> StaticBins;


	///@}
	///@name Pointer Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition
	KRATOS_CLASS_POINTER_DEFINITION(FluxBasedRedistanceProcess);

	///@}
	///@name Life Cycle
	///@{

	FluxBasedRedistanceProcess(
		ModelPart &rBaseModelPart,
		typename TLinearSolver::Pointer pLinearSolver,
		Parameters &Settings) : mrBaseModelPart(rBaseModelPart)
	{
		KRATOS_TRY

		//by default we compute the pseudofilltime
		Parameters default_parameters(R"(
			{
				"echo_level"        : 0
			}  )");
		Settings.ValidateAndAssignDefaults(default_parameters);

		//checking model part is correct;
		Check(rBaseModelPart, Settings);

		// Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
		mPartIsInitialized = false;
		ReGenerateModelPart(rBaseModelPart);

		//compute density to get a Fo=1
		const double domain_length = CalculateDomainLength();
		// ModelPart &rSkinPart = rBaseModelPart.GetSubModelPart(Settings["skin_part"].GetString());
		// const double avg_domain_thickness = 1.0 / static_cast<double>(rSkinPart.NumberOfNodes()) * block_for_each<SumReduction<double>>(mpModelPart->Nodes(), [&](Node<3> &rNode) { return rNode.FastGetSolutionStepValue(THICKNESS); });
		// const double density = (avg_domain_thickness * avg_domain_thickness) / (domain_length * domain_length);
		const double density = 1.0 / (domain_length * domain_length);
		mpModelPart->pGetProcessInfo()->SetValue(DENSITY, density);

		//compute edge size (necessary to fix nodes in a radius)
		mMaxEdgeSize = FindMaxEdgeSize();

		//find neighbours
		// auto &r_comm = mpModelPart->GetCommunicator().GetDataCommunicator();
		// FindGlobalNodalNeighboursProcess global_nodal_neigh(r_comm, *mpModelPart);
		// global_nodal_neigh.Execute();

		// //generate search structure
        // for (ModelPart::NodesContainerType::iterator inode = mpModelPart->NodesBegin(); inode != mpModelPart->NodesEnd(); inode++)
        // {
        //         mpTotalNodes.push_back(*(inode.base()));
        // }
        // mpSearchStructure = new StaticBins(mpTotalNodes.begin(), mpTotalNodes.end());


		// Generate solver
		CreateSolvingStrategy(pLinearSolver, Settings);

		//define solution mode:
		// mSolvePseudoFilltime = false;
		// if (Settings["problem_to_solve"].GetString() == "pseudo_filltime")
		// 	mSolvePseudoFilltime = true;

		KRATOS_CATCH("")
	}

	/// Destructor.
	virtual ~FluxBasedRedistanceProcess()
	{
	}

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	//TODO Charlie: This has to gone if you derive from process.
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "FluxBasedRedistanceProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "FluxBasedRedistanceProcess";
	}

	/**
		 * @brief Computes the pseudo filltime/flowlength by solving two implicit systems of equations
		 * the first system returns a potential field used to compute the pseudovelocity
		 * and the second system uses the previous field to compute the pseudo filltime
		 */
	void Execute() override
	{
		KRATOS_TRY;

		const double max_edge_size = mMaxEdgeSize;

		//reading settings.
		// std::vector<array_1d<double, 3>> InjectionCoordinates;
		// SetVectorsOfInjectionPoints(rInjectionSettings, InjectionCoordinates);

		// //list of nodes that were fixed. At least 1 per inlet
		// GlobalPointersVector<Node<3>> InjectionNodes;

		//Step 0:  Fix distance in injection conditions.
		// VariableUtils().ApplyFixity(DISTANCE, false, mpModelPart->Nodes());			 //unfixing  all nodes
		// VariableUtils().SetHistoricalVariableToZero(DISTANCE, mpModelPart->Nodes()); //unfixing  all nodes
		// for (unsigned int i_injection = 0; i_injection < InjectionCoordinates.size(); i_injection++){
		// 	//find closest node(s)

		// 	double squaredDistance = 0.0;
		//     Node<3> injection_point(1, InjectionCoordinates[i_injection][0],InjectionCoordinates[i_injection][1],InjectionCoordinates[i_injection][2]);
		//     NodePointerType closest_node = mpSearchStructure->SearchNearestPoint(injection_point, squaredDistance);
		// 	KRATOS_ERROR_IF(sqrt(squaredDistance) > max_edge_size) << "Injection point: " << InjectionCoordinates[i_injection] << "lies outside the domain" << std::endl;

		// 	//fix injection node
		// 	closest_node->Fix(DISTANCE);
		// 	InjectionNodes.push_back(closest_node);
		// }

        const int nnodes = static_cast<int>(mpModelPart->NumberOfNodes());

        block_for_each(mpModelPart->Nodes(), [](Node<3>& rNode){
            double& d = rNode.FastGetSolutionStepValue(DISTANCE);
            double& fix_flag = rNode.FastGetSolutionStepValue(FLAG_VARIABLE);

            // Free the DISTANCE values
            fix_flag = 1.0;
            rNode.Free(DISTANCE);

            // Save the distances
            rNode.SetValue(DISTANCE, d);

            if(d == 0){
                d = 1.0e-15;
                fix_flag = -1.0;
                rNode.Fix(DISTANCE);
            } else {
                if(d > 0.0){
                    d = 1.0e15; // Set to a large number, to make sure that that the minimal distance is computed according to CaculateTetrahedraDistances
                } else {
                    d = -1.0e15;
                }
            }
        });

        block_for_each(mpModelPart->Elements(), [this](Element& rElem){
            array_1d<double,TDim+1> distances;
            auto& geom = rElem.GetGeometry();

            for(unsigned int i=0; i<TDim+1; i++){
                distances[i] = geom[i].GetValue(DISTANCE);
            }

            const array_1d<double,TDim+1> original_distances = distances;

            // The element is cut by the interface
            if(this->IsSplit(distances)){
                // Compute the unsigned distance using GeometryUtils
                GeometryUtils::CalculateExactDistancesToPlane(geom, distances);
                // if (mOptions.Is(CALCULATE_EXACT_DISTANCES_TO_PLANE)) {
                //     GeometryUtils::CalculateExactDistancesToPlane(geom, distances);
                // }
                // else {
                //     if(TDim==3){
                //         GeometryUtils::CalculateTetrahedraDistances(geom, distances);
                //     }
                //     else {
                //         GeometryUtils::CalculateTriangleDistances(geom, distances);
                //     }
                // }

                // Assign the sign using the original distance values
                for(unsigned int i = 0; i < TDim+1; ++i){
                    if(original_distances[i] < 0){
                        distances[i] = -distances[i];
                    }
                }

                for(unsigned int i = 0; i < TDim+1; ++i){
                    double &d = geom[i].FastGetSolutionStepValue(DISTANCE);
                    double &fix_flag = geom[i].FastGetSolutionStepValue(FLAG_VARIABLE);
                    geom[i].SetLock();
                    if(std::abs(d) > std::abs(distances[i])){
                        d = distances[i];
                    }
                    fix_flag = -1.0;
                    geom[i].Fix(DISTANCE);
                    geom[i].UnSetLock();
                }
            }
        });

        // SHALL WE SYNCHRONIZE SOMETHING IN HERE?¿?¿??¿ WE'VE CHANGED THE NODAL DISTANCE VALUES FROM THE ELEMENTS...
        this->SynchronizeFixity();
        this->SynchronizeDistance();

        // Compute the maximum and minimum distance for the fixed nodes
        double max_dist = 0.0;
        double min_dist = 0.0;
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = mpModelPart->NodesBegin() + i_node;
            if(it_node->IsFixed(DISTANCE)){
                const double& d = it_node->FastGetSolutionStepValue(DISTANCE);
                if(d > max_dist){
                    max_dist = d;
                }
                if(d < min_dist){
                    min_dist = d;
                }
            }
        }

        // Synchronize the maximum and minimum distance values
        const auto &r_communicator = mpModelPart->GetCommunicator().GetDataCommunicator();
        max_dist = r_communicator.MaxAll(max_dist);
        min_dist = r_communicator.MinAll(min_dist);

        // Assign the max dist to all of the non-fixed positive nodes
        // and the minimum one to the non-fixed negatives
        block_for_each(mpModelPart->Nodes(), [&min_dist, &max_dist](Node<3>& rNode){
            if(!rNode.IsFixed(DISTANCE)){
                double& d = rNode.FastGetSolutionStepValue(DISTANCE);
                if(d>0){
                    d = max_dist;
                } else {
                    d = min_dist;
                }
            }
        });

		// Step1 - solve a transient diffusion problem to get the potential field.
		//This step is the same, whether filltime or flowlength is required
		mpModelPart->pGetProcessInfo()->SetValue(FRACTIONAL_STEP, 1);
		mpSolvingStrategy->Solve();

		//FOR DEBUGGING results of step 1)
		// block_for_each(mpModelPart->Nodes(), [&](Node<3> &rNode) {
		// 	rNode.SetValue(DISTANCE, rNode.FastGetSolutionStepValue(DISTANCE));
		// });

		//step 2: compute velocity using the gradient of the potential field:
		//This step is the same, whether filltime or flowlength is required
		ComputeVelocities();

		//assinging and fixing the values in nodes that are close to the injection point(s)
		// for (unsigned int i_injection = 0; i_injection < InjectionCoordinates.size(); i_injection++){
		// 	VariableUtils().SetFlag(VISITED, false, mpModelPart->Nodes());
		// 	ExpandNodeFixity(InjectionNodes[i_injection], InjectionCoordinates[i_injection], 2.05 * mMaxEdgeSize);
		// }

		//seting the FS # according to the problem to be solved
		// if (mSolvePseudoFilltime){
		// 	mpModelPart->pGetProcessInfo()->SetValue(FRACTIONAL_STEP, 3);
		// }
		// else{ //flowlenght
		// 	mpModelPart->pGetProcessInfo()->SetValue(FRACTIONAL_STEP, 2);
		// }
		mpModelPart->pGetProcessInfo()->SetValue(FRACTIONAL_STEP, 2);

		//solving filltime/flowlength system
		mpSolvingStrategy->Solve();

		KRATOS_CATCH("")
	}



	void Clear()
	{
		mpModelPart->Nodes().clear();
		mpModelPart->Conditions().clear();
		mpModelPart->Elements().clear();
		// mpModelPart->GetProcessInfo().clear();
		mPartIsInitialized = false;

		mpSolvingStrategy->Clear();
	}

	///@}

protected:
	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	ModelPart &mrBaseModelPart;
	ModelPart *mpModelPart;
    //KdtreeType *mpSearchStructure;
	// StaticBins *mpSearchStructure;
	NodePointerTypeVector mpTotalNodes;

	double mMaxEdgeSize;
	bool mPartIsInitialized;
	// bool mSolvePseudoFilltime;
	typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

	///@}
	///@name Protected Operators
	///@{

	/// Constructor without linear solver for derived classes
	FluxBasedRedistanceProcess(
		ModelPart &rBaseModelPart) : mrBaseModelPart(rBaseModelPart)
	{
		mPartIsInitialized = false;
	}

	///@}
	///@name Protected Operations
	///@{

	virtual void ReGenerateModelPart(ModelPart &rBaseModelPart)
	{
		KRATOS_TRY

		Model &current_model = rBaseModelPart.GetModel();

		if (current_model.HasModelPart("RedistanceModelPart"))
			current_model.DeleteModelPart("RedistanceModelPart");
		mpModelPart = &(current_model.CreateModelPart("RedistanceModelPart"));

		// Ensure that the nodes have distance as a DOF
		VariableUtils().AddDof<Variable<double>>(DISTANCE, rBaseModelPart);

		//varibles
		MergeVariableListsUtility().Merge(*mpModelPart, rBaseModelPart);  

		// Generating the model part
	    Element::Pointer p_element = Kratos::make_intrusive<DistanceCalculationFluxBasedElement<TDim>>();
		ConnectivityPreserveModeler modeler;
		modeler.GenerateModelPart(rBaseModelPart, *mpModelPart, *p_element);

		mPartIsInitialized = true;

		KRATOS_CATCH("")
	}

	/**
		 * @brief Computes the nodal velocities 
		 * The direction of the velocity is defined as the gradient of the DISTANCE
		 * While the modulus of the velocity is proportional to its thickness.
		 * NOTE: It is not a real velocity since it does not take into account flowrate or flowing section
		 * It is only meant to provide a faster moving flow front in thicker regions.
		 */
	void ComputeVelocities()
	{
		KRATOS_TRY

		const ProcessInfo &rCurrentProcessInfo = mpModelPart->GetProcessInfo();

		//not using variable utils to do the two tasks in the same loop
		block_for_each(mpModelPart->Nodes(), [&](Node<3> &rNode) {
			rNode.FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
			rNode.FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
		});

		block_for_each(mpModelPart->Elements(), [&](ModelPart::ElementType &rElement) {
			rElement.AddExplicitContribution(rCurrentProcessInfo);
		});

		//not using variable utils to do the two tasks in the same loop
		block_for_each(mpModelPart->Nodes(), [&](Node<3> &rNode) {
			array_1d<double, 3> &vel = rNode.FastGetSolutionStepValue(VELOCITY);

			//dividing by nodal volume in case we want to the get gradient with original modulus
			vel /= rNode.FastGetSolutionStepValue(NODAL_VOLUME);

			//in case we want to scale with thickness using nodal values instead of elemental values.
			//getting |vel|=1
			//vel /= MathUtils<double>::Norm3(vel);
			//scaling using the thickness
			//const double thickness = rNode.FastGetSolutionStepValue(THICKNESS);
			//const double new_vel_modulus = thickness*thickness;
			//vel *= new_vel_modulus;
		});

		KRATOS_CATCH("")
	}

	//MISC FUNCTIONS

	void Check(ModelPart &rBaseModelPart, Parameters &Settings)
	{
		// Check that there is at least one element and node in the model
		const auto n_nodes = rBaseModelPart.NumberOfNodes();
		const auto n_elems = rBaseModelPart.NumberOfElements();

		KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
		KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

		VariableUtils().CheckVariableExists<Variable<double>>(DISTANCE, rBaseModelPart.Nodes());
		VariableUtils().CheckVariableExists<Variable<double>>(NODAL_VOLUME, rBaseModelPart.Nodes());
		// VariableUtils().CheckVariableExists<Variable<double>>(THICKNESS, rBaseModelPart.Nodes());
		VariableUtils().CheckVariableExists<Variable<array_1d<double, 3>>>(VELOCITY, rBaseModelPart.Nodes());

		if (TDim == 2){
			KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) << "In 2D the element type is expected to be a triangle" << std::endl;
		} else if (TDim == 3) {
			KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) << "In 3D the element type is expected to be a tetrahedra" << std::endl;
		}

		//checking the problem to be solved is correct
		// auto problem_to_solve = Settings["problem_to_solve"].GetString();
		// KRATOS_ERROR_IF_NOT(problem_to_solve == "pseudo_filltime" || problem_to_solve == "flowlength") << "'problem_to_solve' not recognized. Accepted options: [pseudo_filltime,flowlength]" << std::endl;
	}

	void CreateSolvingStrategy(typename TLinearSolver::Pointer pLinearSolver, Parameters &Settings)
	{
		typename SchemeType::Pointer pscheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>>();
		typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

		bool CalculateReactions = false;
		bool ReformDofAtEachIteration = false;
		bool CalculateNormDxFlag = false;

		BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pLinearSolver);
		mpSolvingStrategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>>(
			*mpModelPart,
			pscheme,
			pBuilderSolver,
			CalculateReactions,
			ReformDofAtEachIteration,
			CalculateNormDxFlag);

		mpSolvingStrategy->SetEchoLevel(Settings["echo_level"].GetInt());

		mpSolvingStrategy->Check();
	}

	double CalculateDomainLength()
	{
		typedef CombinedReduction<MinReduction<double>,
								  MinReduction<double>,
								  MinReduction<double>,
								  MaxReduction<double>,
								  MaxReduction<double>,
								  MaxReduction<double>>
			MultipleReduction;
		double x_min, y_min, z_min, x_max, y_max, z_max;

		std::tie(
			x_min,
			y_min,
			z_min,
			x_max,
			y_max,
			z_max) = block_for_each<MultipleReduction>(mpModelPart->Nodes(), [&](Node<3> &rNode) {
			const array_1d<double, 3> coord = rNode.Coordinates();
			return std::make_tuple(
				coord[0],
				coord[1],
				coord[2],
				coord[0],
				coord[1],
				coord[2]);
		});

		return std::sqrt((x_max - x_min) * (x_max - x_min) + (y_max - y_min) * (y_max - y_min) + (z_max - z_min) * (z_max - z_min));
	}

	///@}

	//////////////// Auxiliary Functions
	double FindMaxEdgeSize()
	{
		KRATOS_TRY;
		return block_for_each<MaxReduction<double>>(mpModelPart->Elements(), [&](ModelPart::ElementType &rElement) {
			Geometry<Node<3>> &rGeometry = rElement.GetGeometry();
			return rGeometry.MaxEdgeLength();
		});
		KRATOS_CATCH("")
	}

	// void ExpandNodeFixity(Node<3> &rNode, const array_1d<double, 3> &rInjectionPoint, const double &Radius)
	// {
	// 	const double this_node_dist = MathUtils<double>::Norm3(rNode.Coordinates() - rInjectionPoint);
	// 	const bool is_not_visited = rNode.IsNot(VISITED);
	// 	if (this_node_dist <= Radius && is_not_visited)	{
	// 		rNode.Set(VISITED, true);
	// 		double old_nodal_value = std::numeric_limits<double>::max(); // setting a high value
	// 		if (rNode.IsFixed(DISTANCE)){ //if it is already fixed, we must take the lowest value
	// 			old_nodal_value = rNode.FastGetSolutionStepValue(DISTANCE);
	// 		} else { //fixing it.
	// 			rNode.Fix(DISTANCE);
	// 		}

	// 		//assigning the right value according to the problem to be solved.
	// 		double this_inlet_nodal_value = this_node_dist; //for flowlenght, we take directly the distance to the injection point
	// 		if (mSolvePseudoFilltime) {
	// 			const double nodal_vel = MathUtils<double>::Norm3(rNode.FastGetSolutionStepValue(VELOCITY));
	// 			this_inlet_nodal_value /= nodal_vel; //for filltime, we divide by velocity to compute the time it takes to reach this node
	// 		}

	// 		//taking the lowest and saving
	// 		if (this_inlet_nodal_value < old_nodal_value) {
	// 			rNode.FastGetSolutionStepValue(DISTANCE) = this_inlet_nodal_value;
	// 		}

	// 		//looping the neighbouring nodes of the nodes
	// 		GlobalPointersVector<Node<3>> &r_neighbour_nodes = rNode.GetValue(NEIGHBOUR_NODES); //Calcular otra vez?
	// 		for (GlobalPointersVector<Node<3>>::iterator i_neighbour_node = r_neighbour_nodes.begin(); i_neighbour_node != r_neighbour_nodes.end(); i_neighbour_node++)	{
	// 			ExpandNodeFixity(*i_neighbour_node, rInjectionPoint, Radius);
	// 		}
	// 	}
	// }

    void SynchronizeDistance(){
        auto &r_communicator = mpModelPart->GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(mpModelPart->NumberOfNodes());

            // Set the distance absolute value
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mpModelPart->NodesBegin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = std::abs(it_node->FastGetSolutionStepValue(DISTANCE));
            }

            // Synchronize the unsigned value to minimum
            r_communicator.SynchronizeCurrentDataToMin(DISTANCE);

            // Set the distance sign again by retrieving it from the non-historical database
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mpModelPart->NodesBegin() + i_node;
                if(it_node->GetValue(DISTANCE) < 0.0){
                    it_node->FastGetSolutionStepValue(DISTANCE) = -it_node->FastGetSolutionStepValue(DISTANCE);
                }
            }
        }
    }

    void SynchronizeFixity(){
        auto &r_communicator = mpModelPart->GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(mpModelPart->NumberOfNodes());

            // Synchronize the fixity flag variable to minium
            // (-1.0 means fixed and 1.0 means free)
            r_communicator.SynchronizeCurrentDataToMin(FLAG_VARIABLE);

            // Set the fixity according to the synchronized flag
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mpModelPart->NodesBegin() + i_node;
                const double &r_fix_flag = it_node->FastGetSolutionStepValue(FLAG_VARIABLE);
                if (r_fix_flag == -1.0){
                    it_node->Fix(DISTANCE);
                }
            }
        }
    }

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

    bool IsSplit(const array_1d<double,TDim+1> &rDistances){
        unsigned int positives = 0, negatives = 0;

        for(unsigned int i = 0; i < TDim+1; ++i){
            if(rDistances[i] >= 0){
                ++positives;
            } else {
                ++negatives;
            }
        }

        if (positives > 0 && negatives > 0){
            return true;
        }

        return false;
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

	/// Assignment operator.
	FluxBasedRedistanceProcess &operator=(FluxBasedRedistanceProcess const &rOther);

	///@}

}; // Class FluxBasedRedistanceProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function

///@}

} // namespace Kratos.

#endif //