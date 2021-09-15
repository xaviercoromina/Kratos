//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED)
#define  KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/modified_mises_nonlocal_damage_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ModifiedMisesNonlocalDamage3DLawMix : public ModifiedMisesNonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMisesNonlocalDamage3DLawMix);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ModifiedMisesNonlocalDamage3DLawMix();
    
    /// Second Constructor
    ModifiedMisesNonlocalDamage3DLawMix(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    ModifiedMisesNonlocalDamage3DLawMix (const ModifiedMisesNonlocalDamage3DLawMix& rOther);

    /// Destructor
    ~ModifiedMisesNonlocalDamage3DLawMix() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */

    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient ) override;



private:
    
    /// Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class ModifiedMisesNonlocalDamage3DLawMix
}  // namespace Kratos.
#endif // KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED  defined 
