//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/io.h"
#include "includes/model_part.h"
#include "containers/variable.h"
#include "containers/flags.h"
#include "containers/container_expression/container_expression.h"
#include "utilities/xml_utilities/xml_nodes.h"
#include "utilities/xml_utilities/xml_entities.h"


namespace Kratos {
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    KRATOS_DEFINE_LOCAL_FLAG( NODES );
    KRATOS_DEFINE_LOCAL_FLAG( CONDITIONS );
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS );

    ///@}
    ///@name Life cycle
    ///@{

    VtuOutput() {}

    ///@}
    ///@name Public operations
    ///@{

    void SetEntities(
        const ModelPart& rModelPart,
        const Flags& EntityFlags,
        const bool IsInitialConfiguration = true,
        const IndexType Precision = 9);

    template<class TDataType>
    void AddHistoricalVariable(
        const Variable<TDataType>& rVariable,
        const IndexType Precision = 9);

    template<class TDataType>
    void AddNonHistoricalVariable(
        const Variable<TDataType>& rVariable,
        const Flags& EntityFlags,
        const IndexType Precision = 9);

    void AddFlagVariable(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        const Flags& EntityFlags);

    template <class TContainerType>
    void AddContainerExpression(
        const std::string& rExpressionName,
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression,
        const IndexType Precision = 9);

    void Write(const std::string& rOutputFileName) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart const* mpModelPart = nullptr;

    XmlElement::Pointer mpXmlVtkFile = nullptr;

    XmlElement::Pointer mpXmlPiece = nullptr;

    XmlNodes::Pointer mpXmlPoints = nullptr;

    typename XmlEntities<ModelPart::ConditionType>::Pointer mpXmlConditions = nullptr;

    typename XmlEntities<ModelPart::ElementType>::Pointer mpXmlElements = nullptr;

    XmlElement::Pointer mpXmlPointData = nullptr;

    XmlElement::Pointer mpXmlCellData = nullptr;

    ///@}
    ///@name Private operations
    ///@{

    template<class TEntityType>
    XmlElement::Pointer GetEntityData();

    template<class TEntityType, class TDataWriterType>
    void VtuOutput::AddEntityData(
        const std::string& rDataName,
        const IndexType NumberOfComponents,
        const TDataWriterType& rDataWriter,
        const IndexType Precision);

    ///@}
};
} // namespace Kratos