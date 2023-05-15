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
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) XmlWriter {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    virtual ~XmlWriter() = default;

    ///@}
    ///@name Public operations
    ///@{

    virtual void WriteElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const bool IsEmptyElement) = 0;

    virtual void CloseElement(
        const std::string& rTagName,
        const IndexType Level) = 0;

    virtual void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<Expression::Pointer>& rExpressions,
        const IndexType Level) = 0;

    ///@}
};

} // namespace Kratos