// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra,
//                   Anne van de Graaf
//
#pragma once

#include "processes/process.h"


namespace Kratos
{

class ModelPart;
class Parameters;


class ApplyVectorConstraintsTableProcess : public Process
{
public:
    ApplyVectorConstraintsTableProcess(ModelPart&        rModelPart,
                                       const Parameters& rSettings);

    ~ApplyVectorConstraintsTableProcess() override;

    using ProcessUniquePointer = Kratos::unique_ptr<Process>;

    void ExecuteInitialize() override;
    void ExecuteInitializeSolutionStep() override;

private:
    std::vector<Parameters> CreateParametersForActiveComponents(const Parameters& rSettings) const;
    std::vector<char> ActiveComponents(const Parameters& rSettings) const;
    Parameters CreateParametersForComponent(const Parameters& rSettings, char component) const;
    static std::size_t ComponentToIndex(char component);
    ProcessUniquePointer MakeProcessFor(const Parameters& rParameters) const;

    ModelPart& mrModelPart;
    std::vector<ProcessUniquePointer> mProcesses;
};

}
