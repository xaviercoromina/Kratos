//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "entity_specific_properties_process.h"

namespace Kratos {

const Parameters EntitySpecificPropertiesProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
            "container_type" : "elements",
            "echo_level"     : 0
        })");

    return default_parameters;
}

EntitySpecificPropertiesProcess::EntitySpecificPropertiesProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();

    const auto& container_type = rParameters["container_type"].GetString();
    if (container_type == "conditions") {
        mIsConditions = true;
    } else if (container_type == "elements") {
        mIsElements = true;
    } else if (container_type == "all") {
        mIsConditions = true;
        mIsElements = true;
    } else {
        KRATOS_ERROR << "EntitySpecificPropertiesProcess: Unsupported "
                        "container type requested for "
                     << mModelPartName << ". Followings are the supported container types:"
                     << "\n\tconditions"
                     << "\n\telements"
                     << "\n\tall";
    }

    KRATOS_CATCH("");
}

int EntitySpecificPropertiesProcess::Check()
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

void EntitySpecificPropertiesProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (!mIsPropertiesReplaced) {
        if (mIsConditions) {
            CreateEntitySpecificProperties(r_model_part.Conditions());
        }

        if (mIsElements) {
            CreateEntitySpecificProperties(r_model_part.Elements());
        }

        mIsPropertiesReplaced = true;

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Created entity specific properties for"
            << (mIsConditions ? " conditions" : "")
            << (mIsElements ? " elements" : "") << " in " << mModelPartName << ".\n";
    }

    KRATOS_CATCH("");
}

std::string EntitySpecificPropertiesProcess::Info() const
{
    return std::string("EntitySpecificPropertiesProcess");
}

void EntitySpecificPropertiesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void EntitySpecificPropertiesProcess::PrintData(std::ostream& rOStream) const
{
}

template<class ContainerType>
void EntitySpecificPropertiesProcess::CreateEntitySpecificProperties(ContainerType& rContainer)
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    // creation of properties is done in serial
    int properties_id = r_model_part.NumberOfProperties() + 1;
    for (auto& r_entity : rContainer) {
        auto p_properties = r_model_part.CreateNewProperties(properties_id++);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        r_entity.SetProperties(p_properties);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
