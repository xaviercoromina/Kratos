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
    ~ApplyK0ProcedureProcess() override = default;

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
   
  private:
      ModelPart& mrModelPart;
      std::string mModelPartName;

      void CalculateK0Stresses(Element& rElement)
      {
          // Get K0 material parameters of this element ( probably there is something more efficient )
          const Element::PropertiesType& rProp = rElement.GetProperties();
          ConstitutiveLaw::Pointer pConstitutiveLaw = rProp.GetValue(CONSTITUTIVE_LAW);
          const int K0MainDirection = rProp[K0_MAIN_DIRECTION];
          if (K0MainDirection < 0 || K0MainDirection > 1) {
              KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in ApplyK0ProcedureProcess: " << K0MainDirection << std::endl;
          }
          KRATOS_INFO("apply_k0_procedure_process") << " K0_NC " << rProp.Has(K0_NC) << " index voor Phi  " << rProp.Has(NUMBER_OF_UMAT_PHI_PARAMETER) << std::endl;

          //Check for alternative K0 specifications
          double PoissonUR = 0.;
          if (rProp.Has(POISSON_UNLOADING_RELOADING)) PoissonUR = rProp[POISSON_UNLOADING_RELOADING];
          array_1d<double, 3> K0Vector;
          if (rProp.Has(K0_NC)) {
              std::fill(K0Vector.begin(), K0Vector.end(), rProp[K0_NC]);
           }
          else if (rProp.Has(NUMBER_OF_UMAT_PHI_PARAMETER) && rProp.Has(NUMBER_OF_UMAT_PARAMETERS) && rProp.Has(UMAT_PARAMETERS)) {
              if (rProp[NUMBER_OF_UMAT_PHI_PARAMETER] < 1 || rProp[NUMBER_OF_UMAT_PHI_PARAMETER] > rProp[NUMBER_OF_UMAT_PARAMETERS]) {
                  KRATOS_ERROR << "undefined NUMBER_OF_UMAT_PHI_PARAMETER in ApplyK0ProcedureProcess: " << rProp[NUMBER_OF_UMAT_PHI_PARAMETER] << std::endl;
              }
              // is more checking is possible and should that happen here?
              double Phi = rProp[UMAT_PARAMETERS][rProp[NUMBER_OF_UMAT_PHI_PARAMETER] - 1];
              if (Phi < 0. || Phi > 90.) {
                  KRATOS_ERROR << "friction angle Phi out of range in ApplyK0ProcedureProcess: " << Phi << std::endl;
              }
              // convert to Radians Use the conversion function of Anne i.s.o. the nonsens to Radians below.
              Phi *= Globals::Pi / 180.;
              std::fill(K0Vector.begin(), K0Vector.end(), 1.0 - sin(Phi));
          }
          else if (rProp.Has(K0_VALUE_XX) && rProp.Has(K0_VALUE_YY) && rProp.Has(K0_VALUE_ZZ)) {
              K0Vector[0] = rProp[K0_VALUE_XX];
              K0Vector[1] = rProp[K0_VALUE_YY];
              K0Vector[2] = rProp[K0_VALUE_ZZ];
          }
          else {
              KRATOS_ERROR << "Insufficient material dat for K0 procedure process: " << std::endl;
          }

          //Loop over integration points
          const Element::GeometryType& rGeom = rElement.GetGeometry();
          const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(rElement.GetIntegrationMethod());
          
          // Get element stress tensor 
          const ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();
          std::vector<ConstitutiveLaw::StressVectorType> rStressVectors;
          rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);

          for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {

              // Determine OCR dependent K0 values
              if ((rProp.Has(K0_NC) || rProp.Has(NUMBER_OF_UMAT_PHI_PARAMETER)) && rProp.Has(OCR)) {
                  //Modify for presence of OCR (or POP?) field values
                  double K0Value = rProp[K0_NC] * rProp[OCR] - (PoissonUR / (1.0 - PoissonUR)) * (rProp[OCR] - 1.0);
                  std::fill(K0Vector.begin(), K0Vector.end(), K0Value);
              }

              // Apply K0 procedure
             for (int IDir = 0; IDir <= 2; ++IDir) {
                  if (IDir != K0MainDirection) {
                      rStressVectors[GPoint][IDir] = K0Vector[IDir] * rStressVectors[GPoint][K0MainDirection];
                  }
              }
             // Erase shear stresses
             std::fill(rStressVectors[GPoint].begin()+3, rStressVectors[GPoint].end(), 0.0);
           }
          // Set element integration point stress tensors
          rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);

      }

}; //Class

} /* namespace Kratos.*/
