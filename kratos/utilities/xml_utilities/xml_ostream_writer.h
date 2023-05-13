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
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/xml_utilities/xml_writer.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) XmlOStreamWriter : public XmlWriter
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    XmlOStreamWriter(
        std::ostream& rOStream,
        const IndexType Precision);

    ///@}
    ///@name Public operations
    ///@{

    void WriteElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const bool IsEmptyElement) override;

    void CloseElement(
        const std::string& rTagName,
        const IndexType Level) override;

    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<Expression::Pointer>& rExpressions,
        const std::vector<IndexType> rNumberOfEntities,
        const IndexType Level) override;

    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<LiteralFlatExpression<int>::Pointer>& rExpressions,
        const std::vector<IndexType> rNumberOfEntities,
        const IndexType Level) override;

    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<LiteralFlatExpression<double>::Pointer>& rExpressions,
        const std::vector<IndexType> rNumberOfEntities,
        const IndexType Level) override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream;

    ///@}
    ///@name Private operations
    ///@{

    void WriteAttributes(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level);

    static const std::string GetTabbing(const IndexType Level);

    ///@}
};

} // namespace Kratos