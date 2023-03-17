//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/simo_ju_local_damage_3D_law_mix_ortho.hpp"

namespace Kratos
{

//Default Constructor
SimoJuLocalDamage3DLawMixOrtho::SimoJuLocalDamage3DLawMixOrtho() : SimoJuLocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuLocalDamage3DLawMixOrtho::SimoJuLocalDamage3DLawMixOrtho(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : SimoJuLocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuLocalDamage3DLawMixOrtho::SimoJuLocalDamage3DLawMixOrtho(const SimoJuLocalDamage3DLawMixOrtho& rOther) : SimoJuLocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuLocalDamage3DLawMixOrtho::~SimoJuLocalDamage3DLawMixOrtho() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer SimoJuLocalDamage3DLawMixOrtho::Clone() const
{
    SimoJuLocalDamage3DLawMixOrtho::Pointer p_clone(new SimoJuLocalDamage3DLawMixOrtho(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void SimoJuLocalDamage3DLawMixOrtho::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();

    // TODO: DCB_test. We calculated C matrix in Matlab/Excel
    // const double E_1      = 1.38e11;
    // const double E_2      = 8.96e9;
    // const double E_3      = E_2;
    // const double G_12     = 7.1e9;
    // const double G_23     = 3.446e9;
    // const double nu_12    = 0.3;
    // const double nu_23    = nu_12;
    // const double nu_13    = nu_12;
    // const double nu_21    = E_2/E_1*nu_12;
    // const double nu_32    = E_3/E_2*nu_23;
    // const double nu_31    = E_3/E_1*nu_13;
    // const double G_31_inv = (1.0+nu_31)/E_1+(1.0+nu_13)/E_3;
    // const double G_31     = 1.0/G_31_inv;

    // 3D linear elastic constitutive matrix

    Matrix TempLinearElasticMatrix(6,6);
    noalias(TempLinearElasticMatrix) = ZeroMatrix(6,6);
    TempLinearElasticMatrix ( 0 , 0 ) = 140343119915.104;
    TempLinearElasticMatrix ( 0 , 1 ) = 3905199858.50725;
    TempLinearElasticMatrix ( 0 , 2 ) = 3905199858.50725;
    TempLinearElasticMatrix ( 1 , 0 ) = 3905199858.50725;
    TempLinearElasticMatrix ( 1 , 1 ) = 9954820276.99927;
    TempLinearElasticMatrix ( 1 , 2 ) = 3062512584.69157;
    TempLinearElasticMatrix ( 2 , 0 ) = 3905199858.50725;
    TempLinearElasticMatrix ( 2 , 1 ) = 3062512584.69157;
    TempLinearElasticMatrix ( 2 , 2 ) = 9954820276.99926;
    TempLinearElasticMatrix ( 3 , 3 ) = 3446000000.00000; // Using formulas we obtained this (same as in paper)
    TempLinearElasticMatrix ( 4 , 4 ) = 7100000000.00000; // paper data
    // TempLinearElasticMatrix ( 4 , 4 ) = 6558374380.36490; // Using formulas we obtained this...
    TempLinearElasticMatrix ( 5 , 5 ) = 7100000000.00000; // just like D(4,4)
    // TempLinearElasticMatrix ( 5 , 5 ) = 6558374380.36490; // Using formulas we obtained this... (the solution is not stiff enough)

    noalias(rLinearElasticMatrix) = TempLinearElasticMatrix;
}

} // Namespace Kratos
