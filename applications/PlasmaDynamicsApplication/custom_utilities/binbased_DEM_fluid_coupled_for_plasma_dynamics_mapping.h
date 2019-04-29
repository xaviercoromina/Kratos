//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Marc CHUNG TO SANG, mchungtosang@cimne.upc.edu$
//   Date:                $Date: 2019-04-10 18:56:00 $
//
//

#if !defined(KRATOS_BINBASED_DEM_FLUID_COUPLED_FOR_PLASMA_DYNAMICS_MAPPING)
#define  KRATOS_BINBASED_DEM_FLUID_COUPLED_FOR_PLASMA_DYNAMICS_MAPPING

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "mollification/density_function_polynomial.h"
#include "custom_functions.h"
#include "../../SwimmingDEMApplication/custom_utilities/fields/velocity_field.h"
#include "../../SwimmingDEMApplication/custom_utilities/fields/fluid_field_utility.h"

#include "../../SwimmingDEMApplication/custom_utilities/binbased_DEM_fluid_coupled_mapping.h"

//Database includes


// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

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

/// This class allows the interpolation between non-matching simplicial meshes in 2D and 3D with linear shape functions. it is designed for DEM-CFD coupling problems

/*
* @author  Marc CHUNG TO SANG, mchungtosang@cimne.upc.edu
*
* Derived class of BinBasedDEMFluidCoupledMapping for plasma dynamics application: 
*                  BinBasedDEMFluidCoupledPDMapping
*
* @author  Guillermo Casas Gonzalez <gcasas@cimne.upc.edu>
*
* For every node of the destination model part it is checked in which element of the origin model part it is
* contained and a linear interpolation is performed
*
* The data structure used by default is a bin,
*
* For a more general tool that allows the mapping between 2 and 3D non-matching meshes, please see /kratos/applications/MeshingApplication/custom_utilities/projection.h
*/




// Some function definitions
//***************************************************************************************************************
//***************************************************************************************************************


//***************************************************************************************************************
//***************************************************************************************************************

template <std::size_t TDim, typename TBaseTypeOfIonParticle>
class KRATOS_API(PLASMA_DYNAMICS_APPLICATION) BinBasedDEMFluidCoupledPDMapping : public BinBasedDEMFluidCoupledMapping
{
public:
///@name Type Definitions
///@{
typedef ModelPart::ElementsContainerType                      ElementsArrayType;
typedef ElementsArrayType::ContainerType                      ResultElementsContainerType;
typedef std::vector<ResultElementsContainerType>              VectorResultElementsContainerType;
typedef ModelPart::ElementsContainerType::iterator            ElementIteratorType;
typedef IonParticle<TBaseTypeOfIonParticle>  ParticleType;

typedef ModelPart::NodesContainerType                         NodesArrayType;
typedef NodesArrayType::ContainerType                         ResultNodesContainerType;
typedef std::vector<ResultNodesContainerType>                 VectorResultNodesContainerType;
typedef std::vector<Node<3>::Pointer>                         NodalPointersContainerType;
typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;

typedef std::size_t                                           ListIndexType;
typedef SpatialSearch::DistanceType                           DistanceType;
typedef SpatialSearch::VectorDistanceType                     VectorDistanceType;

/// Pointer definition of BinBasedDEMFluidCoupledPDMapping
typedef BinBasedDEMFluidCoupledPDMapping<TDim, TBaseTypeOfIonParticle> BinBasedDEMFluidCoupledPDMapping_TDim_TBaseTypeOfIonParticle;
KRATOS_CLASS_POINTER_DEFINITION(BinBasedDEMFluidCoupledPDMapping_TDim_TBaseTypeOfIonParticle);

///@}
///@name Life Cycle
///@{

/// Default constructor.
//----------------------------------------------------------------
//                       Key for coupling_type
//----------------------------------------------------------------
//        Averaged variables       |  Fluid Fraction
//   Fluid-to-DEM | DEM-to-fluid   |
//----------------------------------------------------------------
// 0:   Linear         Constant            Constant
// 1:   Linear         Linear              Constant
// 2:   Linear         Linear              Linear
// 3:   Linear         Filtered            Filtered
//----------------------------------------------------------------


BinBasedDEMFluidCoupledPDMapping(Parameters& rParameters)
                             : mMustCalculateMaxNodalArea(true),
                               mFluidDeltaTime(0.0),
                               mFluidLastCouplingFromDEMTime(0.0),
                               mMaxNodalAreaInv(0.0),
                               mNumberOfDEMSamplesSoFarInTheCurrentFluidStep(0)
{

}

/// Destructor.
virtual ~BinBasedDEMFluidCoupledPDMapping() {}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{
void InterpolateFromFluidMesh(
        ModelPart& r_fluid_model_part,
        ModelPart& r_dem_model_part,
        Parameters& parameters,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha) override;
void ProjectPlasmaDynamics(Element::Pointer p_elem,
             const Vector& N,
             Node<3>::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha) ;
void GetElectricFieldValue(ModelPart& r_fluid_model_part);


///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> string.
virtual std::string Info() const
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

///@}
///@name Friends
///@{

///@}

protected:



private:




//***************************************************************************************************************
//***************************************************************************************************************

//***************************************************************************************************************
//***************************************************************************************************************


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
BinBasedDEMFluidCoupledPDMapping& operator=(BinBasedDEMFluidCoupledPDMapping const& rOther);

///@}

}; // Class BinBasedDEMFluidCoupledPDMapping

/// output stream function
template<std::size_t TDim, typename TBaseTypeOfIonParticle>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinBasedDEMFluidCoupledPDMapping<TDim, TBaseTypeOfIonParticle>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_BINBASED_DEM_FLUID_COUPLED_FOR_PLASMA_DYNAMICS_MAPPING  defined
