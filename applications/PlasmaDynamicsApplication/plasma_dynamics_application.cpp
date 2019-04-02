
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


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "geometries/sphere_3d_1.h"

#include "custom_constitutive/DEM_electromagnetic_CL.h"

// Application includes
#include "plasma_dynamics_application.h"

namespace Kratos
{

KratosPlasmaDynamicsApplication::KratosPlasmaDynamicsApplication()
    : KratosApplication("PlasmaDynamicsApplication"),
      mIonParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
      mElectronParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1))))
      {}

void KratosPlasmaDynamicsApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("PlasmaDyn") << std::endl;
    KRATOS_INFO("PlasmaDyn") << "     KRATOS |  _  ) |  _ \\  |  _ \\ __ _  ___| | __      "<< std::endl;
    KRATOS_INFO("PlasmaDyn") << "            | |_| ) | | | | | |_) / _` |/ __| |/ /      "<< std::endl;
    KRATOS_INFO("PlasmaDyn") << "            |  ___) | |_| | |  __/ (_| | (__|   <       "<< std::endl;
    KRATOS_INFO("PlasmaDyn") << "            |_|     |____/  |_|   \\__,_|\\___|_|\\_\\      " << std::endl;
    KRATOS_INFO("PlasmaDyn") << std::endl;
    KRATOS_INFO("PlasmaDyn") << "Importing PlasmaDynamicsApplication... ";
  
    KRATOS_INFO("") << "Initializing KratosPlasmaDynamicsApplication... " << std::endl;

    //Register Elements
    KRATOS_REGISTER_ELEMENT("IonParticle3D", mIonParticle3D)
    KRATOS_REGISTER_ELEMENT("ElectronParticle3D", mElectronParticle3D)

    // Register Constitutive Laws
    Serializer::Register("DEM_electromagnetic", DEM_electromagnetic());
    KRATOS_INFO("") << " done." << std::endl;
}

}// namespace Kratos.
