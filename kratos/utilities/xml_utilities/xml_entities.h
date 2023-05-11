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
#include "utilities/xml_utilities/xml_nodes.h"

namespace Kratos {

template<class TEntityType>
class XmlEntities : public XmlElement {
public:
    ///@name Type definitions
    ///@{

    using BaseType = XmlElement;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlEntities);

    ///@}
    ///@name Life cycle
    ///@{

    XmlEntities(
        const XmlNodes::Pointer pXmlNodes,
        const std::vector<TEntityType const*>& rEntities);

    ///@}
    ///@name Public operations
    ///@{

    const std::string GetDataName() const;

    const XmlNodes::Pointer GetXmlNodes() const;

    const std::vector<TEntityType const*>& GetEntities() const;

    ///@}

private:
    ///@name Private classes
    ///@{

    // forward declaration
    class XmlConnectivities;

    class XmlOffsets : public XmlElement
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(XmlOffsets);

        ///@}
        ///@name Life cycle
        ///@{

        XmlOffsets(const typename XmlConnectivities::Pointer pXmlConnectivities);

        ///@}
    private:
        ///@name Private member variables
        ///@{

        const typename XmlConnectivities::Pointer mpXmlConnectivities;

        ///@}
        ///@name Private opeations
        ///@{

        bool HasElementData() const override;

        void WrtieElementData(
            std::ostream& rOutputStream,
            const std::string& rTabbing) const override;

        ///@}
    };

    class XmlGeometryTypes : public XmlElement
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(XmlGeometryTypes);

        ///@}
        ///@name Life cycle
        ///@{

        XmlGeometryTypes(const typename XmlConnectivities::Pointer pXmlConnectivities);

        ///@}
    private:
        ///@name Private member variables
        ///@{

        const typename XmlConnectivities::Pointer mpXmlConnectivities;

        ///@}
        ///@name Private opeations
        ///@{

        bool HasElementData() const override;

        void WrtieElementData(
            std::ostream& rOutputStream,
            const std::string& rTabbing) const override;

        ///@}
    };

    class XmlConnectivities : public XmlElement
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(XmlConnectivities);

        ///@}
        ///@name Life cycle
        ///@{

        XmlConnectivities(
            const XmlNodes::Pointer pXmlNodes,
            const std::vector<TEntityType const*>& rEntities);

        ///@}
    private:
        ///@name Private member variables
        ///@{

        const XmlNodes::Pointer mpXmlNodes;

        const std::vector<TEntityType const*> mEntities;

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

        template<class T>
        friend class XmlEntities;

        friend class XmlOffsets;

        friend class XmlGeometryTypes;

        ///@}
    };

    ///@}
    ///@name Private member variables
    ///@{

    typename XmlConnectivities::Pointer mpXmlConnectivities;

    typename XmlOffsets::Pointer mpXmlOffsets;

    typename XmlGeometryTypes::Pointer mpXmlGeometryTypes;

    ///@}
};
} // namespace Kratos
