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
#include <utility>
#include <vector>
#include <algorithm>

// External includes
#include "includes/define.h"

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/xml_utilities/xml_ostream_writer.h"

// Include base h
#include "xml_element.h"

namespace Kratos {

XmlElement::XmlElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

XmlElement::XmlElement(
    const std::string& rDataName,
    const std::vector<Expression::Pointer>& rExpressions,
    const std::vector<IndexType>& rNumberOfEntities)
    : mTagName("DataArray"),
      mExpressions(rExpressions),
      mNumberOfEntities(rNumberOfEntities)
{
    KRATOS_ERROR_IF(rExpressions.size() == 0)
        << "Empty expression lists are not allowed.";

    KRATOS_ERROR_IF(rExpressions.size() != rNumberOfEntities.size())
        << "Expressions list and size of number of entities list mismatch "
        << "[ number of expressions = " << mExpressions.size()
        << ", size of number entities list = "
        << mNumberOfEntities.size() << " ].";

    std::string data_type;
    IndexType number_of_components = 0;

    for (IndexType i = 0; i < rExpressions.size(); ++i) {
        Expression* p_expression = &*(rExpressions[i]);
        std::string current_data_type;
        if (dynamic_cast<LiteralFlatExpression<int>*>(p_expression)) {
            current_data_type = "Int32";
        } else {
            current_data_type = "Float32";
        }

        if (data_type == "") {
            data_type = current_data_type;
        }

        if (data_type != current_data_type) {
            data_type = "Float32";
        }

        if (number_of_components == 0) {
            number_of_components = p_expression->GetFlattenedSize();
        }

        KRATOS_ERROR_IF(number_of_components != p_expression->GetFlattenedSize())
            << "Found expressions with mismatching shapes.";
    }

    AddAttribute("type", data_type);
    AddAttribute("Name", rDataName);
    AddAttribute("NumberOfComponents", std::to_string(number_of_components));
}

const std::string XmlElement::GetTagName() const
{
    return mTagName;
}

void XmlElement::AddAttribute(
    const std::string& rName,
    const std::string& rValue)
{
    if (mExpressions.size() == 0) {
        for (const auto& r_attribute : mAttributes) {
            KRATOS_ERROR_IF(r_attribute.first == rName)
                << "There exists an attribute named \"" << rName
                << "\" in the xml element with value = \""
                << r_attribute.second << "\" [ given new value = \""
                << rValue << "\" ].\n";
        }
        mAttributes.push_back(std::make_pair(rName, rValue));
    } else {
        KRATOS_ERROR << "Cannot add attributes to an Xml element which has data.\n";
    }
}

const std::vector<std::pair<const std::string,const std::string>>& XmlElement::GetAttributes() const
{
    return mAttributes;
}

void XmlElement::AddElement(const XmlElement::Pointer pXmlElement)
{
    if (mExpressions.size() == 0) {
        for (const auto& p_element : mXmlElements) {
            KRATOS_ERROR_IF(&*(p_element) == &*(pXmlElement))
                << "The xml element is already aded.";
        }
        mXmlElements.push_back(pXmlElement);
    } else {
        KRATOS_ERROR << "Cannot add Xml elements to an Xml element which has data.\n";
    }
}

const std::vector<XmlElement::Pointer>& XmlElement::GetElements() const
{
    return mXmlElements;
}

void XmlElement::ClearElements()
{
    mXmlElements.clear();
}

void XmlElement::Write(
    XmlWriter& rWriter,
    const IndexType Level) const
{
    if (mExpressions.size() == 0 && mXmlElements.size() == 0) {
        rWriter.WriteElement(GetTagName(), GetAttributes(), Level, true);
    } else if (mXmlElements.size() > 0) {
        rWriter.WriteElement(GetTagName(), GetAttributes(), Level, false);
        for (const auto p_element : mXmlElements) {
            p_element->Write(rWriter, Level + 1);
        }
        rWriter.CloseElement(GetTagName(), Level);
    } else if (mExpressions.size() > 0) {
        if (std::all_of(
                mExpressions.begin(),
                mExpressions.end(),
                [](const auto& pExpression) {
                    return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression);})) {
            std::vector<LiteralFlatExpression<int>::Pointer> int_flat_expressions(mExpressions.size());
            std::transform(
                mExpressions.begin(),
                mExpressions.end(),
                int_flat_expressions.begin(),
                [](auto pExpression) {
                    return LiteralFlatExpression<int>::Pointer(dynamic_cast<LiteralFlatExpression<int>*>(&*(pExpression)));
                }
            );
            rWriter.WriteDataElement(GetTagName(), GetAttributes(),  int_flat_expressions, mNumberOfEntities, Level);
        } else if (std::all_of(
                        mExpressions.begin(),
                        mExpressions.end(),
                        [](const auto& pExpression) {
                            return dynamic_cast<LiteralFlatExpression<double>*>(&*pExpression);})) {
            std::vector<LiteralFlatExpression<double>::Pointer> double_flat_expressions(mExpressions.size());
            std::transform(
                mExpressions.begin(),
                mExpressions.end(),
                double_flat_expressions.begin(),
                [](auto pExpression) {
                    return LiteralFlatExpression<double>::Pointer(dynamic_cast<LiteralFlatExpression<double>*>(&*(pExpression)));
                }
            );
            rWriter.WriteDataElement(GetTagName(), GetAttributes(),  double_flat_expressions, mNumberOfEntities, Level);
        } else {
            rWriter.WriteDataElement(GetTagName(), GetAttributes(), mExpressions, mNumberOfEntities, Level);
        }
    }
}

std::string XmlElement::Print(const IndexType Level) const
{
    std::stringstream msg;
    XmlOStreamWriter writer(msg);
    this->Write(writer, Level);
    return msg.str();

}

} // namespace Kratos