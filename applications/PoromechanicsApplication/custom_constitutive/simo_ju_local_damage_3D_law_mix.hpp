//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_SIMO_JU_LOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED)
#define  KRATOS_SIMO_JU_LOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/simo_ju_local_damage_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SimoJuLocalDamage3DLawMix : public SimoJuLocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SimoJuLocalDamage3DLawMix);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SimoJuLocalDamage3DLawMix();
    
    /// Second Constructor
    SimoJuLocalDamage3DLawMix(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    SimoJuLocalDamage3DLawMix (const SimoJuLocalDamage3DLawMix& rOther);

    /// Destructor
    ~SimoJuLocalDamage3DLawMix() override;

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

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

}; // Class SimoJuLocalDamage3DLawMix
}  // namespace Kratos.
#endif // KRATOS_SIMO_JU_LOCAL_DAMAGE_3D_LAW_MIX_H_INCLUDED  defined 
