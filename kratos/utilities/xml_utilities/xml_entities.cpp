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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/xml_utilities/xml_nodes.h"

// Include base h
#include "xml_entities.h"

namespace Kratos {

static const std::map<GeometryData::KratosGeometryType, int> geo_type_vtu_cell_type_map = {
    { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
    { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
    { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
    { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
    { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 },
    { GeometryData::KratosGeometryType::Kratos_Prism3D6,         13 },
    { GeometryData::KratosGeometryType::Kratos_Pyramid3D5,       14 },
    { GeometryData::KratosGeometryType::Kratos_Line2D3,          21 },
    { GeometryData::KratosGeometryType::Kratos_Line3D3,          21 },
    { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      22 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      22 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23 },
    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   24 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    25 },
    { GeometryData::KratosGeometryType::Kratos_Prism3D15,        26 },
    { GeometryData::KratosGeometryType::Kratos_Pyramid3D13,      27 }
};

template<class TEntityType>
XmlEntities<TEntityType>::XmlConnectivities::XmlConnectivities(
    const XmlNodes::Pointer pXmlNodes,
    const std::vector<TEntityType const*>& rEntities)
    : BaseType("DataArray"),
      mpXmlNodes(pXmlNodes),
      mEntities(rEntities)
{
    AddAttribute("type", "Int32");
    AddAttribute("Name", "connectivity");
    AddAttribute("Format", "ascii");
}

template<class TEntityType>
bool XmlEntities<TEntityType>::XmlConnectivities::HasElementData() const
{
    return true;
}

template<class TEntityType>
void XmlEntities<TEntityType>::XmlConnectivities::WrtieElementData(
    std::ostream& rOutputStream,
    const std::string& rTabbing) const
{
    rOutputStream << std::fixed;

    for (const auto p_entity : mEntities) {
        const auto& r_geometry = p_entity->GetGeometry();
        const auto& node_id_map = mpXmlNodes->GetNodeIdMap();

        rOutputStream << rTabbing;

        for (const auto& r_node : r_geometry) {
            auto p_itr = node_id_map.find(&r_node);
            if (p_itr != node_id_map.end()) {
                rOutputStream << "  " << p_itr->second;
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found in the XmlNodes.";
            }
        }

        rOutputStream << "\n";
    }
}

template<class TEntityType>
XmlEntities<TEntityType>::XmlOffsets::XmlOffsets(
    const typename XmlConnectivities::Pointer pXmlConnectivities)
    : BaseType("DataArray"),
      mpXmlConnectivities(pXmlConnectivities)
{
    AddAttribute("type", "Int32");
    AddAttribute("Name", "offsets");
    AddAttribute("Format", "ascii");
}

template<class TEntityType>
bool XmlEntities<TEntityType>::XmlOffsets::HasElementData() const
{
    return true;
}

template<class TEntityType>
void XmlEntities<TEntityType>::XmlOffsets::WrtieElementData(
    std::ostream& rOutputStream,
    const std::string& rTabbing) const
{
    rOutputStream << std::fixed;

    IndexType offset = 0;
    rOutputStream << rTabbing;

    for (const auto p_entity : mpXmlConnectivities->mEntities) {
        offset += p_entity->GetGeometry().size();
        rOutputStream << "  " << offset;
    }

    rOutputStream << "\n";
}

template<class TEntityType>
XmlEntities<TEntityType>::XmlGeometryTypes::XmlGeometryTypes(
    const typename XmlConnectivities::Pointer pXmlConnectivities)
    : BaseType("DataArray"),
      mpXmlConnectivities(pXmlConnectivities)
{
    AddAttribute("type", "Int32");
    AddAttribute("Name", "types");
    AddAttribute("Format", "ascii");
}

template<class TEntityType>
bool XmlEntities<TEntityType>::XmlGeometryTypes::HasElementData() const
{
    return true;
}

template<class TEntityType>
void XmlEntities<TEntityType>::XmlGeometryTypes::WrtieElementData(
    std::ostream& rOutputStream,
    const std::string& rTabbing) const
{
    rOutputStream << std::fixed;

    rOutputStream << rTabbing;

    for (const auto p_entity : mpXmlConnectivities->mEntities) {
        const auto p_itr = geo_type_vtu_cell_type_map.find(p_entity->GetGeometry().GetGeometryType());
        if (p_itr != geo_type_vtu_cell_type_map.end()) {
            rOutputStream << "  " << p_itr->second;
        } else {
            KRATOS_ERROR << "Element with id " << p_entity->Id() << " has unsupported geometry.";
        }
    }

    rOutputStream << "\n";
}

template<class TEntityType>
XmlEntities<TEntityType>::XmlEntities(
    const XmlNodes::Pointer pXmlNodes,
    const std::vector<TEntityType const*>& rEntities)
    : BaseType("Cells")
{
    mpXmlConnectivities = Kratos::make_shared<XmlConnectivities>(pXmlNodes, rEntities);
    mpXmlOffsets = Kratos::make_shared<XmlOffsets>(mpXmlConnectivities);
    mpXmlGeometryTypes = Kratos::make_shared<XmlGeometryTypes>(mpXmlConnectivities);

    AddElement(mpXmlConnectivities);
    AddElement(mpXmlOffsets);
    AddElement(mpXmlGeometryTypes);
}

template<class TEntityType>
const std::string XmlEntities<TEntityType>::GetDataName() const
{
    return "PointData";
}

template<class TEntityType>
const XmlNodes::Pointer XmlEntities<TEntityType>::GetXmlNodes() const
{
    return mpXmlConnectivities->mpXmlNodes;
}

template<class TEntityType>
const std::vector<TEntityType const*>& XmlEntities<TEntityType>::GetEntities() const
{
    return mpXmlConnectivities->mEntities;
}

// template instantiations
template class XmlEntities<ModelPart::ConditionType>;
template class XmlEntities<ModelPart::ElementType>;

} // namespace Kratos
