// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes 

// External includes 

// Project includes 
#include "custom_elements/conv_diff_1d.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"

namespace Kratos
{

ConvDiff1D::ConvDiff1D(
  IndexType NewId, 
  GeometryType::Pointer pGeometry
  ) : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

ConvDiff1D::ConvDiff1D(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    ) : Element(NewId, pGeometry, pProperties)
{

}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer ConvDiff1D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiff1D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer ConvDiff1D::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiff1D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

ConvDiff1D::~ConvDiff1D()
{
}

/***********************************************************************************/
/***********************************************************************************/

void ConvDiff1D::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ConvDiff1D::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, 
    VectorType& rRightHandSideVector, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
      

    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ConvDiff1D::CalculateRightHandSide(
    VectorType& rRightHandSideVector, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ConvDiff1D::EquationIdVector(
    EquationIdVectorType& rResult, 
    const ProcessInfo& rCurrentProcessInfo
    ) const 
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_variable = p_settings->GetUnknownVariable();
    if (rResult.size() != 2)
        rResult.resize(2, false);

    for (unsigned int i = 0; i < 2; i++) {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_variable).EquationId();
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ConvDiff1D::GetDofList(
  DofsVectorType& rElementalDofList, 
  const ProcessInfo& rCurrentProcessInfo
  ) const
{
    KRATOS_TRY;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_variable = p_settings->GetUnknownVariable();
    if (rElementalDofList.size() != 2)
        rElementalDofList.resize(2);

    for (unsigned int i = 0; i < 2; i++) {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_variable);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

int ConvDiff1D::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    return check;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos


