//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_PORO_EXPLICIT_CDF_SCHEME_HPP_INCLUDED)
#define KRATOS_PORO_EXPLICIT_CDF_SCHEME_HPP_INCLUDED

/* External includes */

/* Project includes */
#include "custom_strategies/schemes/poro_explicit_cd_scheme.hpp"
#include "utilities/variable_utils.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos {

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

/**
 * @class PoroExplicitCDFScheme
 * @ingroup StructuralMechanicsApplciation
 * @brief An explicit forward euler scheme with a split of the inertial term
 * @author Ignasi de Pouplana
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class PoroExplicitCDFScheme
    : public PoroExplicitCDScheme<TSparseSpace, TDenseSpace> {

public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef Scheme<TSparseSpace, TDenseSpace> BaseofBaseType;
    typedef PoroExplicitCDScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Some definitions related with the base class
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition fo the node iterator
    typedef typename ModelPart::NodeIterator NodeIterator;

    /// The definition of the numerical limit
    static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

    using BaseType::mDeltaTime;
    using BaseType::mAlpha;
    using BaseType::mBeta;
    using BaseType::mTheta;
    using BaseType::mGCoefficient;

    /// Counted pointer of PoroExplicitCDFScheme
    KRATOS_CLASS_POINTER_DEFINITION(PoroExplicitCDFScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details The PoroExplicitCDFScheme method
     */
    PoroExplicitCDFScheme()
        : PoroExplicitCDScheme<TSparseSpace, TDenseSpace>()
    {

    }

    /** Destructor.
    */
    virtual ~PoroExplicitCDFScheme() {}

    ///@}
    ///@name Operators
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        mDelta = r_current_process_info[DELTA];
        mDeltab = r_current_process_info[DELTA_B];
        mGamma = r_current_process_info[GAMMA];
        mKappa0 = r_current_process_info[KAPPA_0];
        mKappa1 = r_current_process_info[KAPPA_1];
        mAlpha1 = r_current_process_info[RAYLEIGH_ALPHA_1];
        mBeta1 = r_current_process_info[RAYLEIGH_BETA_1];

        KRATOS_CATCH("")
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        ) override
    {
        array_1d<double, 3>& r_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3> displacement_aux;
        noalias(displacement_aux) = r_displacement;
        array_1d<double, 3>& r_displacement_old = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT_OLD);
        array_1d<double, 3>& r_displacement_older = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT_OLDER);
        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);

        double& r_current_water_pressure = itCurrentNode->FastGetSolutionStepValue(WATER_PRESSURE);
        double& r_current_dt_water_pressure = itCurrentNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);      

        const array_1d<double, 3>& r_external_force = itCurrentNode->FastGetSolutionStepValue(EXTERNAL_FORCE);
        const array_1d<double, 3>& r_external_force_old = itCurrentNode->FastGetSolutionStepValue(EXTERNAL_FORCE,1);
        // array_1d<double, 3>& r_external_force_older = itCurrentNode->FastGetSolutionStepValue(EXTERNAL_FORCE_OLDER);
        const array_1d<double, 3>& r_internal_force = itCurrentNode->FastGetSolutionStepValue(INTERNAL_FORCE);
        const array_1d<double, 3>& r_internal_force_old = itCurrentNode->FastGetSolutionStepValue(INTERNAL_FORCE,1);
        array_1d<double, 3>& r_internal_force_older = itCurrentNode->FastGetSolutionStepValue(INTERNAL_FORCE_OLDER);

        std::array<bool, 3> fix_displacements = {false, false, false};
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        // CDF-1d
        // for (IndexType j = 0; j < DomainSize; j++) {
        //     if (fix_displacements[j] == false) {
        //             r_displacement[j] = ( (2.0*(1.0+mDelta)-mDeltaTime*(mAlpha+mDelta*mAlpha1))*nodal_mass*r_displacement[j]
        //                                   + (mDeltaTime*(mAlpha+2.0*mDelta*mAlpha1)-(1.0+mDelta))*nodal_mass*r_displacement_old[j]
        //                                   - mDeltaTime*mDelta*mAlpha1*nodal_mass*r_displacement_older[j]
        //                                   - mDeltaTime*(mBeta+mDelta*mBeta1+mDeltaTime*(mGamma+mDelta*mKappa0))*r_internal_force[j]
        //                                   + mDeltaTime*(mBeta+2.0*mDelta*mBeta1-mDeltaTime*(1.0-mGamma+mDelta*mKappa1))*r_internal_force_old[j]
        //                                   - mDeltaTime*mDelta*mBeta1*r_internal_force_older[j]
        //                                   + mDeltaTime*mDeltaTime*((mGamma+mDelta*mKappa0)*r_external_force[j]+(1.0-mGamma+mDelta*mKappa1)*r_external_force_old[j])
        //                                 ) / ( nodal_mass*(1.0+mDelta) );
        //     }
        // }

        // CDF-1d*
        for (IndexType j = 0; j < DomainSize; j++) {
            if (fix_displacements[j] == false) {
                    r_displacement[j] = ( (2.0+mDelta-mDeltaTime*(mAlpha+mDelta*mAlpha1))*nodal_mass*r_displacement[j]
                                          + (mDeltaTime*(mAlpha+2.0*mDelta*mAlpha1)-(1.0+2.0*mDelta))*nodal_mass*r_displacement_old[j]
                                          + (mDelta-mDeltaTime*mDelta*mAlpha1)*nodal_mass*r_displacement_older[j]
                                          - mDeltaTime*(mBeta+mDelta*mBeta1+mDeltaTime*(mGamma-mDelta*mKappa0))*r_internal_force[j]
                                          + mDeltaTime*(mBeta+2.0*mDelta*mBeta1-mDeltaTime*(1.0-mGamma-mDelta*mKappa1))*r_internal_force_old[j]
                                          - mDeltaTime*mDelta*mBeta1*r_internal_force_older[j]
                                          + mDeltaTime*mDeltaTime*((mGamma-mDelta*mKappa0)*r_external_force[j]+(1.0-mGamma-mDelta*mKappa1)*r_external_force_old[j])
                                        ) / ( nodal_mass );
            }
        }

        // CDF-1
        // for (IndexType j = 0; j < DomainSize; j++) {
        //     if (fix_displacements[j] == false) {
        //             r_displacement[j] = ( (2.0*(1.0-mDelta)-mAlpha*mDeltaTime)*nodal_mass*r_displacement[j]
        //                                   + (mDelta-1.0+mAlpha*mDeltaTime)*nodal_mass*r_displacement_old[j]
        //                                   - mDeltaTime*(mBeta+mDeltaTime*(mGamma-mDelta*mKappa0))*r_internal_force[j]
        //                                   + mDeltaTime*(mBeta-mDeltaTime*(1.0-mGamma-mDelta*mKappa1))*r_internal_force_old[j]
        //                                   + mDeltaTime*mDeltaTime*((mGamma-mDelta*mKappa0)*r_external_force[j]+(1.0-mGamma-mDelta*mKappa1)*r_external_force_old[j])
        //                                 ) / ( nodal_mass*(1.0-mDelta) );
        //     }
        // }

        // CDF-12bd
        // for (IndexType j = 0; j < DomainSize; j++) {
        //     if (fix_displacements[j] == false) {
        //             r_displacement[j] = ( (2.0+mDelta+3.5*mDeltab-mAlpha*mDeltaTime)*nodal_mass*r_displacement[j]
        //                                   + (mAlpha*mDeltaTime-1.0+mDelta-4.0*mDeltab)*nodal_mass*r_displacement_old[j]
        //                                   + (1.5*mDeltab-mDelta)*nodal_mass*r_displacement_older[j]
        //                                   - mDeltaTime*(mBeta+0.5*mDeltaTime*(0.75+mDelta-mDeltab))*r_internal_force[j]
        //                                   + mDeltaTime*(mBeta-0.5*mDeltaTime*(1.0+mDelta-mDeltab))*r_internal_force_old[j]
        //                                   - 0.125*mDeltaTime*mDeltaTime*r_internal_force_older[j]
        //                                   + mDeltaTime*mDeltaTime*(0.5*(0.75+mDelta-mDeltab)*r_external_force[j]+0.5*(1.0+mDelta-mDeltab)*r_external_force_old[j]+0.125*r_external_force_older[j])
        //                                 ) / ( nodal_mass*(1.0+mDelta+mDeltab) );
        //     }
        // }

        // Solution of the darcy_equation
        if( itCurrentNode->IsFixed(WATER_PRESSURE) == false ) {
            // TODO: this is on standby
            r_current_water_pressure = 0.0;
            r_current_dt_water_pressure = 0.0;
        }

        noalias(r_displacement_older) = r_displacement_old;
        noalias(r_displacement_old) = displacement_aux;
        // noalias(r_external_force_older) = r_external_force_old;
        noalias(r_internal_force_older) = r_internal_force_old;
        const array_1d<double, 3>& r_velocity_old = itCurrentNode->FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double, 3>& r_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& r_acceleration = itCurrentNode->FastGetSolutionStepValue(ACCELERATION);

        noalias(r_velocity) = (1.0/mDeltaTime) * (r_displacement - r_displacement_old);
        noalias(r_acceleration) = (1.0/mDeltaTime) * (r_velocity - r_velocity_old);
    }


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    double mDelta;
    double mDeltab;
    double mGamma;
    double mKappa0;
    double mKappa1;
    double mAlpha1;
    double mBeta1;

    ///@}
    ///@name Protected Structs
    ///@{


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

    ///@}

}; /* Class PoroExplicitCDFScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_PORO_EXPLICIT_CDF_SCHEME_HPP_INCLUDED  defined */
