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

#include "input_utilities.h"


namespace Kratos
{

Parameters makeProjectParametersFrom(const std::string& rProjectFilePath)
{
    std::ifstream t(rProjectFilePath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    Parameters projFile{buffer.str()};
    return projFile;
}

}
