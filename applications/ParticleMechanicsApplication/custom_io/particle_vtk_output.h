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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "input_output/vtk_output.h"

namespace Kratos
{
/** \brief ParticleVtkOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) ParticleVtkOutput : public VtkOutput
{
public:

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Pointer definition of ParticleVtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(ParticleVtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    explicit ParticleVtkOutput(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})" )
        );

    /// Destructor.
    virtual ~ParticleVtkOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " ParticleVtkOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " ParticleVtkOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Write the nodes in the rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodesToFile(
        const ModelPart& rModelPart,
        std::ofstream& rFileStream
        ) const override;

    ///@}

private:
    ///@name Operations
    ///@{

    /**
     * @brief Prints the Properties Id as an integer variable in each element/condition
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rContainer the container which is being output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WritePropertiesIdsToFile(
        const TContainerType& rContainer,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Prints the Ids of the container entities as an integer variable in entity (e.g. node, element, condition)
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rContainer the container which is being output
     * @param DataName name of the data in the vtk file
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteIdsToFile(
        const TContainerType& rContainer,
        const std::string& DataName,
        std::ofstream& rFileStream
        ) const;

    ///@}
};

} // namespace Kratos
