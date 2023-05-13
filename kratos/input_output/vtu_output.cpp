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
#include <vector>

// External includes

// Project includes
#include "includes/io.h"
#include "includes/data_communicator.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/string_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/xml_utilities/xml_element.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "vtu_output.h"

namespace Kratos {

    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, NODES,   1);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, CONDITIONS,  2);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, ELEMENTS, 3);

namespace VtuOutputHelperUtilities {

const std::map<GeometryData::KratosGeometryType, int> KratosVtuGeometryTypes = {
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

Expression::Pointer CreatePositionsExpression(
    const ModelPart::NodesContainerType& rNodes,
    const bool IsInitialConfiguration)
{
    auto p_position_expression = LiteralFlatExpression<double>::Create(rNodes.size(), {3});
    auto& r_position_expression = *p_position_expression;

    if (IsInitialConfiguration) {
        IndexPartition<IndexType>(rNodes.size()).for_each([&r_position_expression, &rNodes](const IndexType Index) {
            const auto& r_coordinates = (rNodes.begin() + Index)->GetInitialPosition().Coordinates();
            const IndexType start_index = Index * 3;
            r_position_expression.SetData(start_index, 0, r_coordinates[0]);
            r_position_expression.SetData(start_index, 1, r_coordinates[1]);
            r_position_expression.SetData(start_index, 2, r_coordinates[2]);
        });
    } else {
        IndexPartition<IndexType>(rNodes.size()).for_each([&r_position_expression, &rNodes](const IndexType Index) {
            const auto& r_coordinates = (rNodes.begin() + Index)->Coordinates();
            const IndexType start_index = Index * 3;
            r_position_expression.SetData(start_index, 0, r_coordinates[0]);
            r_position_expression.SetData(start_index, 1, r_coordinates[1]);
            r_position_expression.SetData(start_index, 2, r_coordinates[2]);
        });
    }

    return p_position_expression;
}

XmlElement::Pointer CreatePointsXmlElement(
    const ModelPart& rModelPart,
    const bool IsInitialConfiguration)
{
    auto p_points_xml_element = Kratos::make_shared<XmlElement>("Points");

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
    const auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

    std::vector<Expression::Pointer> p_position(2);

    p_position[0] = CreatePositionsExpression(r_local_nodes, IsInitialConfiguration);
    p_position[1] = CreatePositionsExpression(r_ghost_nodes, IsInitialConfiguration);

    auto p_positions_xml_element = Kratos::make_shared<XmlElement>(
        "Position", p_position,
        std::vector<IndexType>({r_local_nodes.size(), r_ghost_nodes.size()}));

    p_points_xml_element->AddElement(p_positions_xml_element);
    return p_points_xml_element;
}

template<class TContainerType>
Expression::Pointer CreateOffsetsExpression(
    int& rTotalOffset,
    const TContainerType& rContainer)
{
    auto p_offsets_expression = LiteralFlatExpression<int>::Create(rContainer.size(), {});
    rTotalOffset = 0;
    for (IndexType i = 0; i < rContainer.size(); ++i) {
        rTotalOffset += (rContainer.begin() + i)->GetGeometry().size();
        p_offsets_expression[i] = rTotalOffset;
    }
    return p_offsets_expression;
}

template<class TContainerType>
Expression::Pointer CreateGeometryTypesExpression(const TContainerType& rContainer)
{
    auto p_geometry_expression = LiteralFlatExpression<int>::Create(rContainer.size(), {});
    for (IndexType i = 0; i < rContainer.size(); ++i) {
        const auto p_itr = KratosVtuGeometryTypes.find((rContainer.begin() + i)->GetGeometry().GetGeometryType());
        if (p_itr != KratosVtuGeometryTypes.end()) {
            p_geometry_expression[i] = p_itr->second;
        } else {
            KRATOS_ERROR << "Element with id " << (rContainer.begin() + i)->Id() << " has unsupported geometry.";
        }
    }
    return p_geometry_expression;
}

template<class TContainerType>
Expression::Pointer CreateConnectivityExpression(
    const IndexType TotalOffset,
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto p_connectivity_expression = LiteralFlatExpression<int>::Create(TotalOffset, {});
    auto& r_connectivity_expression = *p_connectivity_expression;
    IndexType local_index = 0;
    for (const auto& r_entity : rContainer) {
        for (const auto& r_node : r_entity.GetGeometry()) {
            const auto p_itr = rKratosVtuIndicesMap.find(r_node.Id());
            if (p_itr != rKratosVtuIndicesMap.end()) {
                r_connectivity_expression[local_index++] = p_itr.second;
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found in nodes list.";
            }
        }
    }
    return p_connectivity_expression;
}

template<class TContainerType>
XmlElement::Pointer CreateCellsXmlElement(
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlElement>("Cells");

    int total_offset;
    auto p_offsets_expression = CreateOffsetsExpression(total_offset, rContainer);
    auto p_connectivity_expression = CreateConnectivityExpression(total_offset, rContainer, rKratosVtuIndicesMap);
    auto p_geometry_type_expression = CreateGeometryTypesExpression(rContainer);

    auto p_connectivity_xml_element = Kratos::make_shared<XmlElement>(
        "connectivity", std::vector<Expression::Pointer>({p_connectivity_expression}),
        std::vector<IndexType>({total_offset}));
    auto p_offsets_xml_element = Kratos::make_shared<XmlElement>(
        "offsets", std::vector<Expression::Pointer>({p_offsets_expression}),
        std::vector<IndexType>({rContainer.size()}));
    auto p_types_xml_element = Kratos::make_shared<XmlElement>(
        "types", std::vector<Expression::Pointer>({p_geometry_type_expression}),
        std::vector<IndexType>({rContainer.size()}));

    p_cells_xml_element->AddElement(p_connectivity_xml_element);
    p_cells_xml_element->AddElement(p_offsets_xml_element);
    p_cells_xml_element->AddElement(p_types_xml_element);
    return p_cells_xml_element;
}

template<class TDataType, class TContainerType, class TContainerDataIOTag>
XmlElement::Pointer CreateVariableDataXmlElement(
    const Variable<TDataType>& rVariable,
    ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        SpecializedContainerExpression<TContainerType, TContainerDataIOTag, MeshType::Local> local_container(rModelPart);
        local_container.Read(rVariable);
        SpecializedContainerExpression<TContainerType, TContainerDataIOTag, MeshType::Ghost> ghost_container(rModelPart);
        ghost_container.Read(rVariable);

        std::vector<Expression::Pointer> expressions({local_container.pGetExpression(), ghost_container.pGetExpression()});
        std::vector<IndexType> number_of_entities({local_container.GetContainer().size(), ghost_container.GetContainer().size()});

        return Kratos::make_shared<XmlElement>(rVariable.Name(), expressions, number_of_entities);
    } else {
        SpecializedContainerExpression<TContainerType, TContainerDataIOTag> local_container(rModelPart);
        local_container.Read(rVariable);

        std::vector<Expression::Pointer> expressions({local_container.pGetExpression()});
        std::vector<IndexType> number_of_entities({local_container.GetContainer().size()});

        return Kratos::make_shared<XmlElement>(rVariable.Name(), expressions, number_of_entities);
    }
}

template<class TContainerType>
Expression::Pointer CreateContainerFlagExpression(
    const TContainerType& rContainer,
    const Flags& rFlag)
{
    auto p_flag_expression = LiteralFlatExpression<int>::Create(rContainer.size(), {});
    auto& r_flag_expression = *p_flag_expression;

    IndexPartition<IndexType>(rContainer.size()).for_each([&r_flag_expression, &rContainer, &rFlag](const IndexType Index) {
        r_flag_expression[Index] = (rContainer.begin() + Index)->Is(rFlag);
    });

    return p_flag_expression;
}

template<class TContainerType>
XmlElement::Pointer CreateFlagDataXmlElement(
    const std::string& rFlagName,
    const Flags& rFlags,
    const ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_communicator = rModelPart.GetCommunicator();
        const auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
        const auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

        std::vector<Expression::Pointer> expressions({CreateContainerFlagExpression(r_local_nodes, rFlags), CreateContainerFlagExpression(r_ghost_nodes, rFlags)});
        std::vector<IndexType> number_of_entities({r_local_nodes.size(), r_ghost_nodes.size()});

        return Kratos::make_shared<XmlElement>(rFlagName, expressions, number_of_entities);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {

        std::vector<Expression::Pointer> expressions({CreateContainerFlagExpression(rModelPart.Conditions(), rFlags)});
        std::vector<IndexType> number_of_entities({rModelPart.NumberOfConditions()});

        return Kratos::make_shared<XmlElement>(rFlagName, expressions, number_of_entities);
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {

        std::vector<Expression::Pointer> expressions({CreateContainerFlagExpression(rModelPart.Elements(), rFlags)});
        std::vector<IndexType> number_of_entities({rModelPart.NumberOfElements()});

        return Kratos::make_shared<XmlElement>(rFlagName, expressions, number_of_entities);
    }
}

void CopyAttributes(
    XmlElement& rOutputElement,
    const XmlElement& rInputElement)
{
    for (const auto& r_attribute_data : rInputElement.GetAttributes()) {
        rOutputElement.AddAttribute(r_attribute_data.first, r_attribute_data.second);
    }
}

void CreatePDataArrays(
    XmlElement& rOutputElement,
    const XmlElement& rInputElement)
{
    for (const auto& p_element : rInputElement.GetElements()) {
        auto new_pdata_array = Kratos::make_shared<XmlElement>("PDataArray");
        CopyAttributes(*new_pdata_array, *p_element);
        rOutputElement.AddElement(new_pdata_array);
    }
}

}; // namespace VtuOutputHelperUtilities

VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    const bool IsInitialConfiguration,
    const IndexType Precision)
    : mrModelPart(rModelPart),
      mIsInitialConfiguration(IsInitialConfiguration),
      mPrecision(Precision)
{
    const auto& r_communicator = rModelPart.GetCommunicator();

    mIsConditionsConsidered = r_communicator.GlobalNumberOfConditions() > 0;
    mIsElementsConsidered = r_communicator.GlobalNumberOfElements() > 0;

    KRATOS_WARNING_IF("VtuOutput", mIsElementsConsidered && mIsConditionsConsidered)
        << "Conditions and Elements vtu output chosen for " << mrModelPart.FullName()
        << " which is not supported. Giving priority to elements.";

    mIsConditionsConsidered = mIsElementsConsidered ? false : mIsConditionsConsidered;

    if (mIsConditionsConsidered || mIsElementsConsidered) {
        // we first always use the local mesh
        IndexType vtu_index = 0;
        for (const auto& r_node : mrModelPart.GetCommunicator().LocalMesh().Nodes()) {
            mKratosVtuIndicesMap[r_node.Id()] = vtu_index++;
        }

        // then we add the ghost mesh
        for (const auto& r_node : mrModelPart.GetCommunicator().GhostMesh().Nodes()) {
            mKratosVtuIndicesMap[r_node.Id()] = vtu_index++;
        }
    }
}

template<class TDataType>
void VtuOutput::AddHistoricalVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_ERROR_IF(mNonHistoricalNodalVariablesList.find(&rVariable) !=
                    mNonHistoricalNodalVariablesList.end())
        << "The same \"" << rVariable.Name()
        << "\" exists in the non-historical nodal variables list.\n";

    mHistoricalVariablesList.insert(&rVariable);
}

template<class TDataType>
void VtuOutput::AddNonHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        KRATOS_ERROR_IF(mHistoricalVariablesList.find(&rVariable) !=
                        mHistoricalVariablesList.end())
            << "The same \"" << rVariable.Name()
            << "\" exists in the historical nodal variables list.\n";

        mNonHistoricalNodalVariablesList.insert(&rVariable);
    } else if (rEntityFlags.Is(CONDITIONS)) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Condition variable \"" << rVariable.Name()
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mNonHistoricalCellVariablesList.insert(&rVariable);
    } else if (rEntityFlags.Is(ELEMENTS)) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Element variable \"" << rVariable.Name()
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mNonHistoricalCellVariablesList.insert(&rVariable);
    }
}

void VtuOutput::AddFlagVariable(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        mNodalFlagsList[rFlagName] = &rFlagVariable;
    } else if (rEntityFlags.Is(CONDITIONS)) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Condition flag \"" << rFlagName
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mCellFlagsList[rFlagName] = &rFlagVariable;
    } else if (rEntityFlags.Is(ELEMENTS)) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Element flag \"" << rFlagName
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mCellFlagsList[rFlagName] = &rFlagVariable;
    }
}

template <class TContainerType>
void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    const typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Conditions container expression \"" << rExpressionName
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Elements container expression \"" << rExpressionName
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
    }

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        KRATOS_ERROR_IF(mPointContainerExpressionsList.find(rExpressionName) !=
                        mPointContainerExpressionsList.end())
            << "A container expression named \"" << rExpressionName
            << "\" already exists in point data.\n";

        mPointContainerExpressionsList[rExpressionName] = pContainerExpression;
    } else {
        KRATOS_ERROR_IF(mCellContainerExpressionsList.find(rExpressionName) !=
                        mCellContainerExpressionsList.end())
            << "A container expression named \"" << rExpressionName
            << "\" already exists in cell data.\n";
        mCellContainerExpressionsList[rExpressionName] = pContainerExpression;
    }
}

// template instantiations
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<int>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<double>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 3>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 4>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 6>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 9>>&);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<int>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<double>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 3>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 4>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 6>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 9>>&, const Flags&);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>(const std::string&, const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer);

} // namespace Kratos