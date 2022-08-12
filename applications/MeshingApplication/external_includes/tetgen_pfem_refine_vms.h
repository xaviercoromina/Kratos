//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.8 $
//




#if !defined(KRATOS_TETGEN_PFEM_MODELER_VMS_H_INCLUDED )
#define  KRATOS_TETGEN_PFEM_MODELER_VMS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application_variables.h"
#include "processes/node_erase_process.h"

#include "spatial_containers/spatial_containers.h"
//#include "containers/bucket.h"
//#include "containers/kd_tree.h"
//#include "external_includes/trigen_refine.h"
#include "tetgen_pfem_refine.h"

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
	/** Detail class definition.
	*/
    class TetGenPfemModelerVms : public TetGenPfemModeler
	{
	public:
		///@name Type Definitions
		///@{

		typedef Node<3> PointType;
		typedef Node<3>::Pointer PointPointerType;
		//typedef PointerVector<PointType>           PointVector;
		typedef std::vector<PointType::Pointer>           PointVector;
		typedef PointVector::iterator PointIterator;
		typedef std::vector<double>               DistanceVector;
		typedef std::vector<double>::iterator     DistanceIterator;

		/// Pointer definition of TetGenPfemModeler
        KRATOS_CLASS_POINTER_DEFINITION(TetGenPfemModelerVms);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
        TetGenPfemModelerVms() : TetGenPfemModeler() {}

		/// Destructor.
        ~TetGenPfemModelerVms() override= default;


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{


		std::string GetCharCodeWithMeshGenerationOptions(const bool add_nodes_option) override {
			if(add_nodes_option) {
				return "rQJYq1.4/20nS";
			} else {
				return "rQJYnS";
			}
		}

		void MarkNodesToEraseIfNeeded(ModelPart& ThisModelPart, const double h_factor) override {
			unsigned int max_results = 50;
			PointVector res(max_results);
			DistanceVector res_distances(max_results);
			Node<3> work_point(0,0.0,0.0,0.0);
 			//if the remove_node switch is activated, we check if the nodes got too close
			PointVector list_of_nodes;

			list_of_nodes.reserve(ThisModelPart.Nodes().size());

			for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
			{
					(list_of_nodes).push_back(*(i_node.base()));
			}

			unsigned int bucket_size = 20;
			kd_tree  nodes_tree1(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);

			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
			{
				//radius means the distance, closer than which no node shall be allowd. if closer -> mark for erasing
				const double radius=h_factor*in->FastGetSolutionStepValue(NODAL_H);

				work_point[0]=in->X();
				work_point[1]=in->Y();
				work_point[2]=in->Z();

				unsigned int  n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(),res_distances.begin(), max_results);
				if (n_points_in_radius>1) {
					if (in->FastGetSolutionStepValue(IS_BOUNDARY)==0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0) {
						//look if we are already erasing any of the other nodes
						double erased_nodes = 0;
						for(auto i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
							erased_nodes += in->Is(TO_ERASE); // TODO: this looks like a bug, why sum several times the info from the same node?

						if( erased_nodes < 1) //we cancel the node if no other nodes are being erased
							in->Set(TO_ERASE,true);
					}
					else if ( (in)->FastGetSolutionStepValue(IS_STRUCTURE)!=1.0) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2 * h_factor)
					{
						//here we loop over the neighbouring nodes and if there are nodes
						//with IS_BOUNDARY=1 which are closer than 0.2*nodal_h from our we remove the node we are considering
						unsigned int k = 0;
						unsigned int counter = 0;
						for(auto i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
						{
							if ( (*i)->FastGetSolutionStepValue(IS_BOUNDARY,1)==1.0 && res_distances[k] < 0.2*radius && res_distances[k] > 0.0 )
							{
								counter += 1;
							}
							k++;
						}
						if(counter > 0)
							in->Set(TO_ERASE,true);
					}
				}
			}
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
		std::string Info() const override{return "";}

		/// Print information about this object.
		void PrintInfo(std::ostream& rOStream) const override{}

		/// Print object's data.
		void PrintData(std::ostream& rOStream) const override{}


		///@}
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables
		///@{


		///@}
		///@name Protected member Variables
		///@{


		///@}
		///@name Protected Operators
		///@{


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


		///@}

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
        TetGenPfemModelerVms& operator=(TetGenPfemModelerVms const& rOther);


		///@}

    }; // Class TetGenPfemModelerVms

	///@}

	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,
        TetGenPfemModelerVms& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream,
        const TetGenPfemModelerVms& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}


}  // namespace Kratos.

#endif // KRATOS_TETGEN_PFEM_MODELER_VMS_H_INCLUDED  defined


