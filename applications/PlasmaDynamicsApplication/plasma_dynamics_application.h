
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Chung To Sang, Ignasi de Pouplana, Guillermo Casas
//


#if !defined(KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "custom_elements/ion_particle.h"
#include "custom_elements/electron_particle.h"
#include "plasma_dynamics_application_variables.h"

namespace Kratos
{

class KRATOS_API(PLASMA_DYNAMICS_APPLICATION) KratosPlasmaDynamicsApplication : public KratosApplication
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(KratosPlasmaDynamicsApplication);

    // Default constructor
    KratosPlasmaDynamicsApplication();

    // Destructor
    virtual ~KratosPlasmaDynamicsApplication() {}


    void Register() override;

    // Turn back information as a string
    std::string Info() const override
    {
        return "KratosPlasmaDynamicsApplication";
    }

    // Print information about this object
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    // Print object's data
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

private:

// Member Variables
const IonParticle mIonParticle3D;
const ElectronParticle mElectronParticle3D;   

// Assignment operator.
KratosPlasmaDynamicsApplication& operator=(KratosPlasmaDynamicsApplication const& rOther);

// Copy constructor.
KratosPlasmaDynamicsApplication(KratosPlasmaDynamicsApplication const& rOther);

}; // Class KratosPlasmaDynamicsApplication
}  // namespace Kratos.

#endif // KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED  defined


