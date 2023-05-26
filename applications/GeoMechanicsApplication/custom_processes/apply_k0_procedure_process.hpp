// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

// System includes
#include <cmath>
#include <iostream>
#include <string>

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyK0ProcedureProcess : public Process
{
  public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyK0ProcedureProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyK0ProcedureProcess(ModelPart&  model_part,
                           Parameters& rParameters) : Process(Flags()), mrModelPart(model_part)
    {
        KRATOS_TRY

        mModelPartName      =  rParameters["model_part_name"].GetString();

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyK0ProcedureProcess() override{}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfElements() > 0) {
             // K0 procedure for the model part:
             block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
                 CalculateK0Stresses(rElement);
              });
         }
 
        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
  protected:
    /// Member Variables

    ModelPart& mrModelPart;
    std::string mModelPartName;

  private:

      void CalculateK0Stresses(Element& rElement)
      {
          // Get K0 material parameters of this element ( probably there is something more efficient )
          const Element::PropertiesType& rProp = rElement.GetProperties();
          ConstitutiveLaw::Pointer pConstitutiveLaw = rProp.GetValue(CONSTITUTIVE_LAW);
          const int& K0MainDirection = rProp[K0_MAIN_DIRECTION];
          const double& K0ValueXX = rProp[K0_VALUE_XX];
          const double& K0ValueYY = rProp[K0_VALUE_YY];
          const double& K0ValueZZ = rProp[K0_VALUE_ZZ];
         
          //Loop over integration points
          const Element::GeometryType& rGeom = rElement.GetGeometry();
          const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(rElement.GetIntegrationMethod());
          
          // Get element stress tensor 
          ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();
          std::vector<ConstitutiveLaw::StressVectorType> rStressVector;
          rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVector, rCurrentProcessInfo);

          for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {

              // Apply K0 procedure
              if (K0MainDirection == 0) {
                  rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_YY] = K0ValueYY * rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_XX];
                  rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_ZZ] = K0ValueZZ * rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_XX];
              }
              else if (K0MainDirection == 1) {
                  rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_XX] = K0ValueXX * rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_YY];
                  rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_ZZ] = K0ValueZZ * rStressVector[GPoint][INDEX_2D_PLANE_STRAIN_YY];
              }
              else {
                  KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in ApplyK0ProcedureProcess: " << K0MainDirection << std::endl;
              }
          }
          // Set element integration point stress tensors
          rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVector, rCurrentProcessInfo);

      }

}; //Class

} /* namespace Kratos.*/
