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
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) XmlElement {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlElement(const std::string& rTagName);

    virtual ~XmlElement() = default;

    ///@}
    ///@name Public operations
    ///@{

    void AddAttribute(
        const std::string& rName,
        const std::string& rValue);

    const std::vector<std::pair<const std::string,const std::string>>& GetAttributes() const;

    void AddElement(const XmlElement::Pointer pXmlElement);

    const std::vector<XmlElement::Pointer>& GetElements() const;

    const std::string GetTagName() const;

    void Write(
        std::ostream& rOutputStream,
        const IndexType Level = 0) const;

    std::string Print(const IndexType Level = 0) const;

    ///@}

protected:
    ///@name Protected operations
    ///@{

    virtual bool HasElementData() const;

    virtual void WrtieElementData(
        std::ostream& rOutputStream,
        const std::string& rTabbing) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mTagName;

    std::vector<XmlElement::Pointer> mXmlElements;

    std::vector<std::pair<const std::string, const std::string>> mAttributes;

    ///@}
};
} // namespace Kratos