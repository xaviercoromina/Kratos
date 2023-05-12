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

#include "input_output/logger.h"
#include "dgeosettlement.h"


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

int KratosGeoSettlement::Run(const std::string&          rWorkingDirectory,
                             const std::string&          rParameterName,
                             std::function<void(char*)>  logCallback,
                             std::function<void(double)> reportProgress,
                             std::function<void(char*)>  reportTextualProgress,
                             std::function<bool()>       shouldCancel)
{
    return 1;
}

}