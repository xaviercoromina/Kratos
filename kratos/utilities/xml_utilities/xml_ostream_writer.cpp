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

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/xml_utilities/xml_writer.h"

// Include base h
#include "xml_ostream_writer.h"

namespace Kratos {

const std::string XmlOStreamWriter::GetTabbing(const IndexType Level)
{
    std::stringstream ss_tabbing;
    for (IndexType i = 0; i < Level; ++i) {
        ss_tabbing << "   ";
    }
    return ss_tabbing.str();
}

XmlOStreamWriter::XmlOStreamWriter(std::ostream& rOStream)
    : mrOStream(rOStream)
{
}

void XmlOStreamWriter::WriteAttributes(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "<" << rTagName;
    if (rAttributes.size() > 0) {
        for (const auto& r_pair : rAttributes) {
            mrOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
        }
    }
}

void XmlOStreamWriter::WriteElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level,
    const bool IsEmptyElement)
{
    WriteAttributes(rTagName, rAttributes, Level);

    if (IsEmptyElement) {
        mrOStream << "/>\n";
    } else {
        mrOStream << ">\n";
    }

}

void XmlOStreamWriter::CloseElement(
    const std::string& rTagName,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "</" << rTagName << ">\n";
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<Expression::Pointer>& rExpressions,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    mrOStream << " Format=\"ascii\">\n";

    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);

    for (IndexType i = 0; i < rExpressions.size(); ++i) {
        auto& r_expression = *rExpressions[i];
        auto number_of_entities = rNumberOfEntities[i];
        const IndexType number_of_components = r_expression.GetFlattenedSize();

        mrOStream << tabbing;
        for (IndexType entity_index = 0; entity_index < number_of_entities; ++entity_index) {
            const IndexType entity_start_index = entity_index * number_of_components;
            for (IndexType component_index = 0; component_index < number_of_components; ++component_index) {
                mrOStream << "  " << r_expression.Evaluate(entity_index, entity_start_index, component_index);
            }
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<LiteralFlatExpression<int>::Pointer>& rExpressions,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    mrOStream << " Format=\"ascii\">\n";

    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);

    for (IndexType i = 0; i < rExpressions.size(); ++i) {
        auto& r_expression = *rExpressions[i];
        auto number_of_entities = rNumberOfEntities[i];
        const IndexType number_of_components = r_expression.GetFlattenedSize();
        const IndexType data_size = number_of_entities * number_of_components;

        mrOStream << tabbing;
        for (IndexType i = 0; i < data_size; ++i) {
            mrOStream << "  " << r_expression[i];
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<LiteralFlatExpression<double>::Pointer>& rExpressions,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    mrOStream << " Format=\"ascii\">\n";

    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);

    for (IndexType i = 0; i < rExpressions.size(); ++i) {
        auto& r_expression = *rExpressions[i];
        auto number_of_entities = rNumberOfEntities[i];
        const IndexType number_of_components = r_expression.GetFlattenedSize();
        const IndexType data_size = number_of_entities * number_of_components;

        mrOStream << tabbing;
        for (IndexType i = 0; i < data_size; ++i) {
            mrOStream << "  " << r_expression[i];
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

} // namespace Kratos