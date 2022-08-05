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

#if !defined(FRICTIONAL_LAW_H_DEFINED )
#define  FRICTIONAL_LAW_H_DEFINED

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "custom_conditions/paired_condition.h"
#include "includes/process_info.h"
#include "includes/mortar_classes.h"

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
 * @class FrictionalLaw
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class defines the base class for frictional laws
 * @details This class does nothing, define derived frictional laws in order to make use of it
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FrictionalLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// Node definition
    typedef Node<3> NodeType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of FrictionalLaw
    KRATOS_CLASS_POINTER_DEFINITION( FrictionalLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FrictionalLaw()
    {
    }

    ///Copy constructor  (not really required)
    FrictionalLaw(const FrictionalLaw& rhs)
    {
    }

    /// Destructor.
    virtual ~FrictionalLaw()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the friction coefficient
     * @param rNode The node where the friction coefficient is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    virtual double GetFrictionCoefficient(
        const NodeType& rNode,
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method computes the threshold value considered for computing friction
     * @param rNode The node where the threshold value is obtained
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    virtual double GetThresholdValue(
        const NodeType& rNode,
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rCondition The condition where the friction is computed
     * @param rCurrentProcessInfo The current instance of the process info
     */
    virtual int Check(
        const PairedCondition& rCondition,
        const ProcessInfo& rCurrentProcessInfo
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
    virtual std::string Info() const
    {
        return "FrictionalLaw";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class FrictionalLaw

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  FrictionalLaw& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const FrictionalLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) KratosComponents<FrictionalLaw >;
KRATOS_API_EXTERN template class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) KratosComponents< Variable<FrictionalLaw::Pointer> >;

void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AddKratosComponent(std::string const& Name, FrictionalLaw const& ThisComponent);
void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AddKratosComponent(std::string const& Name, Variable<FrictionalLaw::Pointer> const& ThisComponent);

#ifdef KRATOS_REGISTER_FRICTIONAL_LAW
#undef KRATOS_REGISTER_FRICTIONAL_LAW
#endif
#define KRATOS_REGISTER_FRICTIONAL_LAW(name, reference) \
    KratosComponents<FrictionalLaw >::Add(name, reference); \
    Serializer::Register(name, reference);

}  // namespace Kratos.

#endif // FRICTIONAL_LAW_H_DEFINED  defined
