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

#pragma once

#include <string>
#include <functional>

#include "includes/kratos_export_api.h"


namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoSettlement
{
public:
    int Run(const std::string&          rWorkingDirectory,
            const std::string&          rParameterName,
            std::function<void(char*)>  logCallback,
            std::function<void(double)> reportProgress,
            std::function<void(char*)>  reportTextualProgress,
            std::function<bool()>       shouldCancel);
};

}
