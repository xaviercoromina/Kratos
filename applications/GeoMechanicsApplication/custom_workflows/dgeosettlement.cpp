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


namespace Kratos
{

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