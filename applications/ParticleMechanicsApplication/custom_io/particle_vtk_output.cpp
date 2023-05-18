//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//

// System includes

// External includes

// Project includes
#include "particle_mechanics_application_variables.h"
#include "particle_vtk_output.h"
#include "includes/kratos_filesystem.h"

namespace Kratos
{

ParticleVtkOutput::ParticleVtkOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : VtkOutput(rModelPart,ThisParameters)
{
    mOutputSettings["output_path"].SetString("Particle_VTK_Output");
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteNodesToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    ) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    if (GetEntityType(rModelPart) == EntityType::ELEMENT) {
        rFileStream << "POINTS " << r_local_mesh.NumberOfElements() << " float\n";
        for (auto itr_element = r_local_mesh.ElementsBegin(); itr_element != r_local_mesh.ElementsEnd(); itr_element++) {
            std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };
            itr_element->CalculateOnIntegrationPoints(MP_COORD, mp_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mp_coord[0], rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    } else if (GetEntityType(rModelPart) == EntityType::CONDITION) {
        rFileStream << "POINTS " << r_local_mesh.NumberOfConditions() << " float\n";
        for (auto itr_condition = r_local_mesh.ConditionsBegin(); itr_condition != r_local_mesh.ConditionsEnd(); itr_condition++) {
            std::vector<array_1d<double, 3>> mpc_coord = { ZeroVector(3) };
            itr_condition->CalculateOnIntegrationPoints(MPC_COORD, mpc_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mpc_coord[0], rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void ParticleVtkOutput::WritePropertiesIdsToFile(
    const TContainerType& rContainer,
    std::ofstream& rFileStream
    ) const
{
    rFileStream << "PROPERTIES_ID" << " 1 " << rContainer.size() << "  int\n";
    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.GetProperties().Id(), rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void ParticleVtkOutput::WriteIdsToFile(
    const TContainerType& rContainer,
    const std::string& DataName,
    std::ofstream& rFileStream
    ) const
{
    rFileStream << DataName << " 1 " << rContainer.size() << "  int\n";
    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.Id(), rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

} // namespace Kratos
