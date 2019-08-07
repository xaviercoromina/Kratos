//
// Author: Marc Chung To Sang mchungtosang@cimne.upc.edu
//          

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "electron_particle.h"

namespace Kratos {

ElectronParticle::~ElectronParticle(){}

ElectronParticle& ElectronParticle::operator=(ElectronParticle const& rOther) {
    
    SphericParticle::operator=(rOther);
        
    mElectronCharge = rOther.mElectronCharge;

    return *this;
}

void ElectronParticle::Initialize(const ProcessInfo& r_process_info) {   
    SphericParticle::Initialize(r_process_info);
    double added_mass_coefficient = 1.0;
    SetMass(added_mass_coefficient * GetDensity() * CalculateVolume());
}

void ElectronParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {   
    SphericParticle::Calculate(rVariable, Output, r_process_info);
    // ELECTRIC_FIELD_PROJECTED_TO_PARTICLE is a common variable to all types of charged particle used in the plasma dynamics application
    // should be read as if we consider a charged particle then we define its ID_type 
    // this ID_type is then used to create a SubModelPart of the spheres_model_part containing all the particles of this type
    if (rVariable == ELECTRIC_FIELD_PROJECTED_TO_PARTICLE){

        Output = 2.0;

    };
}

void ElectronParticle::CalculateCoulombForce(array_1d<double, 3>& Coulomb_force)
{
    array_1d<double, 3 >& ElectricFieldPoisson = this->GetGeometry()[0].FastGetSolutionStepValue(ELECTRIC_FIELD_PROJECTED_TO_PARTICLE);
    //KRATOS_INFO("DEM: DEM: ElectricFieldPoisson")<< ElectricFieldPoisson << std::endl;
    Coulomb_force[0] =  mElectronCharge * (IonParticle::mExternalElectricField[0]+ ElectricFieldPoisson[0] );
    Coulomb_force[1] =  mElectronCharge * (IonParticle::mExternalElectricField[1]+ ElectricFieldPoisson[1] );
    Coulomb_force[2] =  mElectronCharge * (IonParticle::mExternalElectricField[2]+ ElectricFieldPoisson[2] );
    //KRATOS_INFO("DEM: DEM: Coulomb Force")<< Coulomb_force << std::endl;
}


void ElectronParticle::CalculateLaplaceForce(array_1d<double, 3>& Laplace_force)
{
    array_1d<double, 3 >& VelocityPreviousStep = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    //double Density =  this->GetDensity();
    DEM_SET_TO_CROSS_OF_FIRST_TWO_3(VelocityPreviousStep, IonParticle::mExternalMagneticField, Laplace_force)
    Laplace_force *= mElectronCharge;
    //KRATOS_INFO("DEM: DEM: Velocity")<< VelocityPreviousStep << std::endl;
    //KRATOS_INFO("DEM: DEM: Magnetic Field")<< IonParticle::mExternalMagneticField << std::endl;
    //KRATOS_INFO("DEM: DEM: mElectronCharge")<< mElectronCharge << std::endl;
    //KRATOS_INFO("DEM: DEM: Laplace Force")<< Laplace_force << std::endl;
    
    // Laplace_force[0] = coeff * (VelocityPreviousStep[1] * mExternalMagneticField[2] - VelocityPreviousStep[2] * mExternalMagneticField[1]);
    // Laplace_force[1] = coeff * (VelocityPreviousStep[2] * mExternalMagneticField[0] - VelocityPreviousStep[0] * mExternalMagneticField[2]);
    // Laplace_force[2] = coeff * (VelocityPreviousStep[0] * mExternalMagneticField[1] - VelocityPreviousStep[1] * mExternalMagneticField[0]);
}

void ElectronParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) 
{
    SphericParticle::MemberDeclarationFirstStep(r_process_info);
}

    
} // namespace Kratos
