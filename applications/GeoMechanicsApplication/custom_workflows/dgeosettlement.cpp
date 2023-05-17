// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "dgeosettlement.h"
#include "input_output/logger.h"
#include "custom_utilities/input_utilities.h"


namespace Kratos
{

KratosGeoSettlement::KratosGeoSettlement()
{
    KRATOS_INFO("KratosGeoSettlement") << "Setting up Kratos" << std::endl;

    if (!mKernel.IsImported("GeoMechanicsApplication"))
    {
        KRATOS_INFO("KratosGeoSettlement") << "Importing GeoMechanicsApplication" << std::endl;
        mpApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpApp);
    }
}

int KratosGeoSettlement::RunStage(const std::string&          rWorkingDirectory,
                                  const std::string&          rProjectParametersFileName,
                                  std::function<void(char*)>  logCallback,
                                  std::function<void(double)> reportProgress,
                                  std::function<void(char*)>  reportTextualProgress,
                                  std::function<bool()>       shouldCancel)
{
    KRATOS_INFO("KratosGeoSettlement") << "About to run a stage..." << std::endl;

    const auto project_parameters_file_path = rWorkingDirectory + "/" + rProjectParametersFileName;
    const auto project_parameters = makeProjectParametersFrom(project_parameters_file_path);

    KRATOS_INFO("KratosGeoSettlement") << "Parsed project parameters file " << project_parameters_file_path << std::endl;

    const auto model_part_name = project_parameters["solver_settings"]["model_part_name"].GetString();
    ModelPart& model_part = mModel.CreateModelPart(model_part_name);
    model_part.SetBufferSize(2);

    KRATOS_INFO("KratosGeoSettlement") << "Created a model part" << std::endl;

    return 1;
}

}