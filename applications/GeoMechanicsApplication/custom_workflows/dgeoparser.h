// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include <string>
#include <vector>
#include <memory>

namespace Kratos
{
    // Forward declarations
    class Model;
    class ModelPart;
    class Parameters;
    class Process;

    class KratosGeoParser
    {
    public:
        KratosGeoParser() = delete;
        static void parseMaterial(Model& model, const std::string& rFilepath);
        static void parseMesh(ModelPart& model_part, const std::string& rFilepath);
        static Parameters openProjectParamsFile(const std::string& rFilepath);
    };
}