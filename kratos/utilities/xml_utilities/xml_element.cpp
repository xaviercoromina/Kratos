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

// External includes
#include "includes/define.h"

// Project includes

// Include base h
#include "xml_element.h"

namespace Kratos {

XmlElement::XmlElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

void XmlElement::AddAttribute(
    const std::string& rName,
    const std::string& rValue)
{
    mAttributes.push_back(std::make_pair(rName, rValue));
}

const std::vector<std::pair<const std::string,const std::string>>& XmlElement::GetAttributes() const
{
    return mAttributes;
}

const std::vector<XmlElement::Pointer>& XmlElement::GetElements() const
{
    return mXmlElements;
}

const std::string XmlElement::GetTagName() const
{
    return mTagName;
}

void XmlElement::AddElement(const XmlElement::Pointer pXmlElement)
{
    if (!HasElementData()) {
        mXmlElements.push_back(pXmlElement);
    } else {
        KRATOS_ERROR << "Cannot add Xml elements to an Xml element which has data.\n";
    }
}

void XmlElement::Write(
    std::ostream& rOutputStream,
    const IndexType Level) const
{
    std::stringstream ss_tabbing;
    for (IndexType i = 0; i < Level; ++i) {
        ss_tabbing << "   ";
    }
    const std::string& tabbing = ss_tabbing.str();

    rOutputStream << tabbing << "<" << mTagName;

    if (mAttributes.size() > 0) {
        for (const auto& r_pair : mAttributes) {
            rOutputStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
        }
    }

    if (HasElementData()) {
        rOutputStream << ">\n";
        WrtieElementData(rOutputStream, tabbing);
        rOutputStream << tabbing << "</" << mTagName << ">";
    } else {
        if (mXmlElements.size() == 0) {
            rOutputStream << "/>";
        } else {
            rOutputStream << ">";
            for (const auto& p_xml_element : mXmlElements) {
                rOutputStream << "\n";
                p_xml_element->Write(rOutputStream, Level + 1);
            }
            rOutputStream << "\n" << tabbing << "</" << mTagName << ">";
        }
    }
}

std::string XmlElement::Print(const IndexType Level) const
{
    std::stringstream msg;
    this->Write(msg, Level);
    return msg.str();

}

bool XmlElement::HasElementData() const
{
    return false;
}

void XmlElement::WrtieElementData(
    std::ostream& rOutputStream,
    const std::string& rTabbing) const
{
}

} // namespace Kratos