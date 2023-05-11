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
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/xml_utilities/xml_element.h"

namespace Kratos {
///@name Kratos Classes
///@{

class XmlNodes : public XmlElement {
public:
    ///@name Type definitions
    ///@{

    using BaseType = XmlElement;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlNodes);

    ///@}
    ///@name Life cycle
    ///@{

    XmlNodes(
        const std::vector<ModelPart::NodeType const*>& rNodes,
        const bool IsInitialConfiguration = true,
        const IndexType Precision = 8);

    ///@}
    ///@name Public operations
    ///@{

    const std::vector<ModelPart::NodeType const*>& GetEntities() const;

    const std::unordered_map<ModelPart::NodeType const*, IndexType>& GetNodeIdMap() const;

    const std::string GetDataName() const;

    ///@}

private:
    ///@name Private classes
    ///@{

    class XmlPositions : public XmlElement
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(XmlPositions);

        ///@}
        ///@name Life cycle
        ///@{

        XmlPositions(
            const std::vector<ModelPart::NodeType const*>& rNodes,
            const bool IsInitialConfiguration = true,
            const IndexType Precision = 8);

        ///@}
    private:
        ///@name Private member variables
        ///@{

        const bool mIsInitialConfiguration;

        const IndexType mPrecision;

        std::unordered_map<ModelPart::NodeType const*, IndexType> mNodeIdMap;

        std::vector<ModelPart::NodeType const*> mNodes;

        ///@}
        ///@name Private opeations
        ///@{

        bool HasElementData() const override;

        void WrtieElementData(
            std::ostream& rOutputStream,
            const std::string& rTabbing) const override;

        ///@}
        ///@name Frient classes
        ///@{

        friend class XmlNodes;

        ///@}
    };

    ///@}
    ///@name Private member variables
    ///@{

    XmlPositions::Pointer mpXmlPositions;

    ///@}
};

///@}
} // namespace Kratos
