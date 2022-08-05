// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED )
#define  FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED

// System includes

// External includes

// Project includes
#include "custom_frictional_laws/frictional_law.h"

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

/**
 * @class FrictionalLawWithDerivative
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the base class for frictional laws with derivative
 * @details This class does nothing, define derived frictional laws in order to make use of it
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative
    : public FrictionalLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// Define the base class
    typedef FrictionalLaw BaseType;

    /// Node definition
    typedef Node<3> NodeType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// Definition of the derivative data
    typedef DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster> DerivativeDataType;

    /// The definition of the mortar operators
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, true, TNumNodesMaster> MortarConditionMatrices;

    /// Definition of the derivatives array
    typedef array_1d<array_1d<double, TNumNodes>, TDim * (2 * TNumNodes + TNumNodesMaster)> DerivativesArray;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of FrictionalLawWithDerivative
    KRATOS_CLASS_POINTER_DEFINITION( FrictionalLawWithDerivative );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FrictionalLawWithDerivative()
    {
    }

    ///Copy constructor  (not really required)
    FrictionalLawWithDerivative(const FrictionalLawWithDerivative& rhs)
    {
    }

    /// Destructor.
    ~FrictionalLawWithDerivative()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the threshold array considered for computing friction
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     * @return The threshold derivative array considered for computing friction
     */
    virtual array_1d<double, TNumNodes> GetThresholdArray(
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method computes the threshold derivative array considered for computing friction
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     * @param rDerivativeData The reference to the derivative database
     * @param rMortarConditionMatrices The container of the mortar operators
     * @return The threshold derivative array considered for computing friction
     */
    virtual DerivativesArray GetDerivativesThresholdArray(
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const DerivativeDataType& rDerivativeData,
        const MortarConditionMatrices& rMortarConditionMatrices
        );

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
        return "FrictionalLawWithDerivative";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << std::endl;
    }

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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    }

    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class FrictionalLawWithDerivative

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
inline std::istream & operator >>(std::istream& rIStream,
                                  FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>& rThis);

/// output stream function
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<2, 2, false, 2>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 3, false, 3>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 4, false, 4>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 3, false, 4>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 4, false, 3>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<2, 2, true,  2>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 3, true,  3>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 4, true,  4>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 3, true,  4>;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLawWithDerivative<3, 4, true,  3>;

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AddKratosComponent(std::string const& Name, FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster> const& ThisComponent);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2N
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2N
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2N(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<2, 2, false, 2> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 3, false, 3> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 4, false, 4> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4N
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4N
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4N(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 3, false, 4> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3N
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3N
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3N(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 4, false, 3> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2NNV
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2NNV
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_2D2NNV(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<2, 2, true, 2> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3NNV
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3NNV
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3NNV(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 3, true, 3> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4NNV
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4NNV
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4NNV(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 4, true, 4> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4NNV
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4NNV
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D3N4NNV(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 3, true, 4> >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3NNV
#undef KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3NNV
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW_WITH_DERIVATIVES_3D4N3NNV(name, reference) \
    KratosComponents<FrictionalLawWithDerivative<3, 4, true, 3> >::Add(name, reference); \
    Serializer::Register(name, reference);

}  // namespace Kratos.

#endif // FRICTIONAL_LAW_WITH_DERIVATIVE_H_DEFINED  defined
