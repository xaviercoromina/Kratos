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
#include "includes/io.h"
#include "includes/data_communicator.h"
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

    auto p_local_position_expression = CreatePositionsExpression(r_local_nodes, IsInitialConfiguration);
    auto p_ghost_position_expression = CreatePositionsExpression(r_ghost_nodes, IsInitialConfiguration);

    auto p_positions_xml_element = Kratos::make_shared<XmlElement>(
                                            "Position",
                                            {
                                                p_local_position_expression,
                                                p_ghost_position_expression
                                            },
                                            {
                                                r_local_nodes.size(),
                                                r_ghost_nodes.size()
                                            });

    p_points_xml_element->AddElement(p_positions_xml_element);
    return p_points_xml_element;
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
    mHistoricalVariablesList.push_back(&rVariable);
}

template<class TDataType>
void VtuOutput::AddNonHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        mNonHistoricalNodalVariablesList.push_back(&rVariable);
    } else if (rEntityFlags.Is(CONDITIONS)) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Condition variable \"" << rVariable.Name()
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mNonHistoricalCellVariablesList.push_back(&rVariable);
    } else if (rEntityFlags.Is(ELEMENTS)) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Element variable \"" << rVariable.Name()
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mNonHistoricalCellVariablesList.push_back(&rVariable);
    }
}

void VtuOutput::AddFlagVariable(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        mNodalFlagsList.push_back(std::make_pair(rFlagName, &rFlagVariable));
    } else if (rEntityFlags.Is(CONDITIONS)) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Condition flag \"" << rFlagName
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mCellFlagsList.push_back(std::make_pair(rFlagName, &rFlagVariable));
    } else if (rEntityFlags.Is(ELEMENTS)) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Element flag \"" << rFlagName
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
        mCellFlagsList.push_back(std::make_pair(rFlagName, &rFlagVariable));
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

    mContainerExpressionsList.push_back(pContainerExpression);
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