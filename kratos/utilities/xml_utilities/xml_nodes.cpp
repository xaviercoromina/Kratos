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

// Include base h
#include "xml_nodes.h"

namespace Kratos {

XmlNodes::XmlPositions::XmlPositions(
    const std::vector<ModelPart::NodeType const*>& rNodes,
    const bool IsInitialConfiguration,
    const IndexType Precision)
    : BaseType("DataArray"),
      mIsInitialConfiguration(IsInitialConfiguration),
      mPrecision(Precision)
{
    mNodes.resize(rNodes.size());
    for (IndexType i = 0; i < rNodes.size(); ++i) {
        auto p_node = *(rNodes.begin() + i);
        mNodeIdMap[p_node] = i;
        mNodes[i]  = p_node;
    }

    AddAttribute("type", "Float32");
    AddAttribute("name", "Position");
    AddAttribute("NumberOfComponents", "3");
    AddAttribute("Format", "ascii");
}

bool XmlNodes::XmlPositions::HasElementData() const
{
    return true;
}

void XmlNodes::XmlPositions::WrtieElementData(
    std::ostream& rOutputStream,
    const std::string& rTabbing) const
{
    const auto current_precision = rOutputStream.precision();
    rOutputStream << std::setprecision(mPrecision) << std::scientific;

    if (mIsInitialConfiguration) {
        for (const auto p_node : mNodes) {
            const auto& coordinates = p_node->GetInitialPosition();
            rOutputStream << rTabbing << "  " << coordinates[0] << "  " << coordinates[1] << "  " << coordinates[0] << "\n";
        }
    } else {
        for (const auto p_node : mNodes) {
            const auto& coordinates = p_node->Coordinates();
            rOutputStream << rTabbing << "  " << coordinates[0] << "  " << coordinates[1] << "  " << coordinates[0] << "\n";
        }
    }

    rOutputStream << std::setprecision(current_precision) << std::fixed;
}

XmlNodes::XmlNodes(
    const std::vector<ModelPart::NodeType const*>& rNodes,
    const bool IsInitialConfiguration,
    const IndexType Precision)
    : BaseType("Points")
{
    mpXmlPositions = Kratos::make_shared<XmlPositions>(rNodes, IsInitialConfiguration, Precision);
    AddElement(mpXmlPositions);
}

const std::string XmlNodes::GetDataName() const
{
    return "PointData";
}

const std::vector<ModelPart::NodeType const*>& XmlNodes::GetEntities() const
{
    return mpXmlPositions->mNodes;
}

const std::unordered_map<ModelPart::NodeType const*, IndexType>& XmlNodes::GetNodeIdMap() const
{
    return mpXmlPositions->mNodeIdMap;
}

} // namespace Kratos
