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
#include <utility>
#include <set>
#include <unordered_map>

// External includes

// Project includes
#include "includes/io.h"
#include "includes/data_communicator.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/string_utilities.h"
#include "utilities/xml_utilities/xml_entity_data.h"
#include "utilities/pointer_communicator.h"
#include "utilities/global_pointer_utilities.h"

// Include base h
#include "vtu_output.h"

namespace Kratos {

    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, NODES,   1);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, CONDITIONS,  2);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, ELEMENTS, 3);

namespace VtuOutputHelperUtilities {

template<class TDataType, std::size_t... TIndex>
void WriteArray(
    std::ostream& rOutputStream,
    const TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    ((rOutputStream << "  " << rValue[TIndex]), ...);
}

void WriteValue(
    std::ostream& rOutputStream,
    const int Value)
{
    rOutputStream << "  " << Value;
}

void WriteValue(
    std::ostream& rOutputStream,
    const double Value)
{
    rOutputStream << "  " << Value;
}

template<long unsigned int TSize>
void WriteValue(
    std::ostream& rOutputStream,
    const array_1d<double, TSize>& rValue)
{
    WriteArray(rOutputStream, rValue, std::make_index_sequence<TSize>{});
}

void WriteExpression(
    std::ostream& rOutputStream,
    const Expression& rExpression,
    const IndexType EntityIndex)
{
    const IndexType number_of_components = rExpression.GetFlattenedSize();
    const IndexType entity_data_begin_index = EntityIndex * number_of_components;
    for (IndexType i = 0; i < number_of_components; ++i) {
        rOutputStream << "  " << rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
    }
}

void FillGhostNodalValues(
    std::vector<double>& rOutput,
    const DataCommunicator& rDataCommunicator,
    const Expression& rExpression,
    const ModelPart::NodesContainerType& rLocalNodes,
    const ModelPart::NodesContainerType& rGhostNodes)
{
    const IndexType number_of_ghost_nodes = rGhostNodes.size();

    std::vector<int> ghost_indices(number_of_ghost_nodes);
    std::transform(rGhostNodes.begin(), rGhostNodes.end(), ghost_indices.begin(), [](const auto& rNode) { return rNode.Id(); });
    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(rGhostNodes, ghost_indices, rDataCommunicator);

    GlobalPointerCommunicator<ModelPart::NodeType> pointer_comm(rDataCommunicator, gp_list.ptr_begin(), gp_list.ptr_end());

    const IndexType number_of_components = rExpression.GetFlattenedSize();

    auto values_proxy = pointer_comm.Apply(
        [&rLocalNodes, &rExpression, number_of_components](GlobalPointer<ModelPart::NodeType>& rGP) -> std::vector<double> {
            IndexType entity_index;
            for (entity_index = 0; entity_index < rLocalNodes.size(); ++entity_index) {
                if ((rLocalNodes.begin() + entity_index)->Id() == rGP->Id()) {
                    break;
                }
            }

            const IndexType enitity_data_begin_index = entity_index * number_of_components;
            std::vector<double> values(number_of_components);
            for (IndexType i = 0; i < number_of_components; ++i) {
                values[i] = rExpression.Evaluate(entity_index, enitity_data_begin_index, i);
            }

            return values;
        }
    );

    const IndexType number_of_values = number_of_ghost_nodes * number_of_components;

    if (rOutput.size() != number_of_values) {
        rOutput.resize(number_of_values);
    }

    IndexType local_index = 0;
    for(IndexType i = 0; i < number_of_ghost_nodes; ++i) {
        const auto& r_gp_value = values_proxy.Get(gp_list(i));
        for (IndexType j = 0; j < number_of_components; ++j) {
            rOutput[local_index++] = r_gp_value[j];
        }
    }
}

template<class TDataType>
constexpr std::size_t NumberOfComponents() {
    if constexpr(std::is_same_v<TDataType, int> || std::is_same_v<TDataType, bool> || std::is_same_v<TDataType, double>) {
        return 0;
    } else {
        return std::tuple_size_v<typename TDataType::array_type>;
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

template<class TDataType, class TEntityType, class TContainerDataIOTag>
class VariableDataRetrieval
{
public:
    ///@name Type definitions
    ///@{

    using DataType = std::conditional_t<std::is_same_v<TDataType, int>, int, double>;

    ///@}
    ///@name Life cycle
    ///@{

    VariableDataRetrieval(const Variable<TDataType>& rVariable)
        : mpVariable(&rVariable)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    void WriteData(
        std::ostream& rOutput,
        const std::vector<TEntityType const*>& rEntities,
        const std::string& rTabbing) const
    {
        for (const auto& p_entity : rEntities) {
            rOutput << rTabbing;
            WriteValue(rOutput, ContainerDataIO<TContainerDataIOTag>::GetValue(
                                    *p_entity, *mpVariable));
            rOutput << "\n";
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Variable<TDataType>* mpVariable;

    ///@}
};

template<class TEntityType>
class VariableFlagRetrieval
{
public:
    ///@name Type definitions
    ///@{

    using DataType = bool;

    ///@}
    ///@name Life cycle
    ///@{

    VariableFlagRetrieval(const Flags& rFlag)
        : mpFlag(&rFlag)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    void WriteData(
        std::ostream& rOutput,
        const std::vector<TEntityType const*>& rEntities,
        const std::string& rTabbing) const
    {
        for (const auto& p_entity : rEntities) {
            rOutput << rTabbing << p_entity->Is(*mpFlag) << "\n";
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Flags* mpFlag;

    ///@}
};

template<class TContainerType>
class ContainerExpressionRetrievalBase
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    ContainerExpressionRetrievalBase(
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression,
        const std::vector<typename TContainerType::value_type const*>& rEntities)
        : mpContainerExpression(pContainerExpression)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    void WriteData(
        std::ostream& rOutput,
        const std::vector<typename TContainerType::value_type const*>& rEntities,
        const std::string& rTabbing) const
    {
        for (IndexType i = 0; i < rEntities.size(); ++i) {
            rOutput << rTabbing;
            WriteExpression(rOutput, mpContainerExpression->GetExpression(), i);
            rOutput << "\n";
        }
    }

    ///@}

protected:
    ///@name Private member variables
    ///@{

    const typename ContainerExpression<TContainerType>::Pointer mpContainerExpression;

    ///@}
};

template<class TContainerType>
class ContainerExpressionRetrieval: public ContainerExpressionRetrievalBase<TContainerType>
{
public:
    ///@name Type definitions
    ///@{

    using DataType = double;

    using BaseType = ContainerExpressionRetrievalBase<TContainerType>;

    ///@}
    ///@name Life cycle
    ///@{

    ContainerExpressionRetrieval(
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression,
        const std::vector<typename TContainerType::value_type const*>& rEntities,
        const Communicator& rCommunicator)
        : BaseType(pContainerExpression, rEntities)
    {
    }

    ///@}
};

template<>
class ContainerExpressionRetrieval<ModelPart::NodesContainerType>: public ContainerExpressionRetrievalBase<ModelPart::NodesContainerType>
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DataType = double;

    using BaseType = ContainerExpressionRetrievalBase<ModelPart::NodesContainerType>;

    ///@}
    ///@name Life cycle
    ///@{

    ContainerExpressionRetrieval(
        const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer pContainerExpression,
        const std::vector<ModelPart::NodeType const*>& rEntities,
        const Communicator& rCommunicator)
        : BaseType(pContainerExpression, rEntities),
          mrCommunicator(rCommunicator)
    {
        const auto& r_local_nodes = mrCommunicator.LocalMesh().Nodes();
        for (IndexType i = 0; i < r_local_nodes.size(); ++i) {
            mLocalNodeIndexMap[&*(r_local_nodes.begin() + i)] = i;
        }

        const auto& r_ghost_nodes = mrCommunicator.GhostMesh().Nodes();
        for (IndexType i = 0; i < r_ghost_nodes.size(); ++i) {
            mGhostNodeIndexMap[&*(r_ghost_nodes.begin() + i)] = i;
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    void WriteData(
        std::ostream& rOutput,
        const std::vector<ModelPart::NodeType const*>& rNodes,
        const std::string& rTabbing) const
    {
        const auto& r_expression = mpContainerExpression->GetExpression();
        const IndexType number_of_components = r_expression.GetFlattenedSize();

        std::vector<double> mGhostValues;

        FillGhostNodalValues(mGhostValues, mrCommunicator.GetDataCommunicator(),
                             r_expression, mrCommunicator.LocalMesh().Nodes(),
                             mrCommunicator.GhostMesh().Nodes());

        for (const auto& p_node : rNodes) {
            rOutput << rTabbing;
            const auto p_local_itr = mLocalNodeIndexMap.find(p_node);
            if (p_local_itr != mLocalNodeIndexMap.end()) {
                WriteExpression(rOutput, r_expression, p_local_itr->second);
            } else {
                const auto p_ghost_itr = mGhostNodeIndexMap.find(p_node);
                if (p_ghost_itr != mGhostNodeIndexMap.end()) {
                    const IndexType entity_data_begin_index = p_ghost_itr->second * number_of_components;
                    for (IndexType i = 0; i < number_of_components; ++i) {
                        rOutput << "  " << mGhostValues[entity_data_begin_index + i];
                    }
                } else {
                    KRATOS_ERROR << "Node with id " << p_node->Id()
                                 << " not found in the list of output nodes.\n";
                }
            }
            rOutput << "\n";
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Communicator& mrCommunicator;

    std::unordered_map<ModelPart::NodeType const*, IndexType> mLocalNodeIndexMap;

    std::unordered_map<ModelPart::NodeType const*, IndexType> mGhostNodeIndexMap;

    ///@}

};

}; // namespace VtuOutputHelperUtilities

void VtuOutput::SetEntities(
    const ModelPart& rModelPart,
    const Flags& EntityFlags,
    const bool IsInitialConfiguration,
    const IndexType Precision)
{
    KRATOS_ERROR_IF(mpModelPart)
        << "The VtuOutput I/O is already initialized with " << mpModelPart->FullName()
        << " [ new model part name = " << rModelPart.FullName() << " ].\n";

    mpModelPart = &rModelPart;

    const auto& r_communicator = rModelPart.GetCommunicator();

    const bool consider_nodes = EntityFlags.Is(NODES) && r_communicator.GlobalNumberOfNodes() > 0;
    const bool consider_conditions = EntityFlags.Is(CONDITIONS) && r_communicator.GlobalNumberOfConditions() > 0;
    const bool consider_elements = EntityFlags.Is(ELEMENTS) && r_communicator.GlobalNumberOfElements() > 0;

    KRATOS_WARNING_IF("VtuOutput", consider_elements && consider_conditions)
        << "Conditions and Elements vtu output chosen for " << mpModelPart->FullName()
        << " which is not supported. Giving priority to elements.";

    using p_node_type = ModelPart::NodeType const*;

    mpXmlVtkFile = Kratos::make_shared<XmlElement>("VTKFile");
    mpXmlVtkFile->AddAttribute("type", "UnstructuredGrid");
    mpXmlVtkFile->AddAttribute("version", "0.1");
    mpXmlVtkFile->AddAttribute("byte_order", "BigEndian");

    if (consider_elements) {
        const IndexType number_of_entities = rModelPart.Elements().size();
        std::set<p_node_type> nodes_set;
        std::vector<ModelPart::ElementType const*> elements(number_of_entities);
        for (IndexType i = 0; i < number_of_entities; ++i) {
            const auto& r_element = *(rModelPart.ElementsBegin() + i);
            for (const auto& r_node : r_element.GetGeometry()) {
                nodes_set.insert(&r_node);
            }
            elements[i] = &r_element;
        }

        std::vector<p_node_type> nodes(nodes_set.size());
        std::copy(nodes_set.begin(), nodes_set.end(), nodes.begin());

        mpXmlPoints = Kratos::make_shared<XmlNodes>(nodes, IsInitialConfiguration, Precision);
        mpXmlElements = Kratos::make_shared<XmlEntities<ModelPart::ElementType>>(mpXmlPoints, elements);

    } else if (consider_conditions) {
        const IndexType number_of_entities = rModelPart.Conditions().size();
        std::set<p_node_type> nodes_set;
        std::vector<ModelPart::ConditionType const*> conditions(number_of_entities);
        for (IndexType i = 0; i < number_of_entities; ++i) {
            const auto& r_condition = *(rModelPart.ConditionsBegin() + i);
            for (const auto& r_node : r_condition.GetGeometry()) {
                nodes_set.insert(&r_node);
            }
            conditions[i] = &r_condition;
        }

        std::vector<p_node_type> nodes(nodes.size());
        std::copy(nodes_set.begin(), nodes_set.end(), nodes.begin());

        mpXmlPoints = Kratos::make_shared<XmlNodes>(nodes, IsInitialConfiguration, Precision);
        mpXmlConditions = Kratos::make_shared<XmlEntities<ModelPart::ConditionType>>(mpXmlPoints, conditions);

    } else if (consider_nodes) {
        std::vector<p_node_type> nodes;
        for (const auto& r_node : rModelPart.Nodes()) {
            nodes.push_back(&r_node);
        }
        mpXmlPoints = Kratos::make_shared<XmlNodes>(nodes, IsInitialConfiguration, Precision);

    } else {
        KRATOS_ERROR << "No nodes, conditions or elements found or selected for output.\n";
    }

    // create piece
    mpXmlPiece = Kratos::make_shared<XmlElement>("Piece");
    mpXmlPiece->AddAttribute("NumberOfPoints", std::to_string(mpXmlPoints->GetEntities().size()));
    mpXmlPiece->AddElement(mpXmlPoints);

    if (mpXmlElements) {
        mpXmlPiece->AddAttribute("NumberOfCells", std::to_string(mpXmlElements->GetEntities().size()));
        mpXmlPiece->AddElement(mpXmlElements);
    } else if (mpXmlConditions) {
        mpXmlPiece->AddAttribute("NumberOfCells", std::to_string(mpXmlConditions->GetEntities().size()));
        mpXmlPiece->AddElement(mpXmlConditions);
    }

    // create unstructured grid
    auto unstructured_grid = Kratos::make_shared<XmlElement>("UnstructuredGrid");
    unstructured_grid->AddElement(mpXmlPiece);

    mpXmlVtkFile->AddElement(unstructured_grid);
}

template<>
XmlElement::Pointer VtuOutput::GetEntityData<ModelPart::NodeType>()
{
    KRATOS_ERROR_IF_NOT(mpXmlPoints)
        << "No points xml element found to create the entity data xml element.";

    if (mpXmlPointData) {
        return mpXmlPointData;
    }

    mpXmlPointData = Kratos::make_shared<XmlElement>("PointData");
    mpXmlPiece->AddElement(mpXmlPointData);
    return mpXmlPointData;
}

template<>
XmlElement::Pointer VtuOutput::GetEntityData<ModelPart::ConditionType>()
{
    KRATOS_ERROR_IF_NOT(mpXmlConditions)
        << "No conditions xml element found to create the entity data xml element.";

    if (mpXmlCellData) {
        return mpXmlCellData;
    }

    mpXmlCellData = Kratos::make_shared<XmlElement>("CellData");
    mpXmlPiece->AddElement(mpXmlCellData);
    return mpXmlCellData;
}

template<>
XmlElement::Pointer VtuOutput::GetEntityData<ModelPart::ElementType>()
{
    KRATOS_ERROR_IF_NOT(mpXmlElements)
        << "No elements xml element found to create the entity data xml element.";

    if (mpXmlCellData) {
        return mpXmlCellData;
    }

    mpXmlCellData = Kratos::make_shared<XmlElement>("CellData");
    mpXmlPiece->AddElement(mpXmlCellData);
    return mpXmlCellData;
}

template<class TEntityType, class TDataWriterType>
void VtuOutput::AddEntityData(
    const std::string& rDataName,
    const IndexType NumberOfComponents,
    const TDataWriterType& rDataWriter,
    const IndexType Precision)
{
    auto xml_entity_data = GetEntityData<TEntityType>();
    if constexpr(std::is_same_v<TEntityType, ModelPart::NodeType>) {
        xml_entity_data->AddElement(Kratos::make_shared<XmlEntityData<XmlNodes, TDataWriterType>>(
            mpXmlPoints, rDataName, NumberOfComponents, rDataWriter, Precision));
    } else if constexpr(std::is_same_v<TEntityType, ModelPart::ConditionType>) {
        xml_entity_data->AddElement(Kratos::make_shared<XmlEntityData<XmlEntities<TEntityType>, TDataWriterType>>(
            mpXmlConditions, rDataName, NumberOfComponents, rDataWriter, Precision));
    } else if constexpr(std::is_same_v<TEntityType, ModelPart::ElementType>) {
        xml_entity_data->AddElement(Kratos::make_shared<XmlEntityData<XmlEntities<TEntityType>, TDataWriterType>>(
            mpXmlElements, rDataName, NumberOfComponents, rDataWriter, Precision));
    }
}

template<class TDataType>
void VtuOutput::AddHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const IndexType Precision)
{
    using data_retriever_type =
        VtuOutputHelperUtilities::VariableDataRetrieval<TDataType, ModelPart::NodeType, ContainerDataIOTags::Historical>;
    AddEntityData<TDataType, data_retriever_type>(
        rVariable.Name(), VtuOutputHelperUtilities::NumberOfComponents<TDataType>(),
        data_retriever_type(rVariable), Precision);
}

template<class TDataType>
void VtuOutput::AddNonHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const Flags& EntityFlags,
    const IndexType Precision)
{
    if (EntityFlags.Is(NODES)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableDataRetrieval<TDataType, ModelPart::NodeType, ContainerDataIOTags::NonHistorical>;
        AddEntityData<TDataType, data_retriever_type>(
            rVariable.Name(), VtuOutputHelperUtilities::NumberOfComponents<TDataType>(),
            data_retriever_type(rVariable), Precision);
    }

    if (EntityFlags.Is(CONDITIONS)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableDataRetrieval<TDataType, ModelPart::ConditionType, ContainerDataIOTags::NonHistorical>;
        AddEntityData<TDataType, data_retriever_type>(
            rVariable.Name(), VtuOutputHelperUtilities::NumberOfComponents<TDataType>(),
            data_retriever_type(rVariable), Precision);
    }

    if (EntityFlags.Is(ELEMENTS)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableDataRetrieval<TDataType, ModelPart::ElementType, ContainerDataIOTags::NonHistorical>;
        AddEntityData<TDataType, data_retriever_type>(
            rVariable.Name(), VtuOutputHelperUtilities::NumberOfComponents<TDataType>(),
            data_retriever_type(rVariable), Precision);
    }
}

void VtuOutput::AddFlagVariable(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    const Flags& EntityFlags)
{
    if (EntityFlags.Is(NODES)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableFlagRetrieval<ModelPart::NodeType>;
        AddEntityData<bool, data_retriever_type>(
            rFlagName, 1, data_retriever_type(rFlagVariable), 1);
    }

    if (EntityFlags.Is(CONDITIONS)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableFlagRetrieval<ModelPart::ConditionType>;
        AddEntityData<bool, data_retriever_type>(
            rFlagName, 1, data_retriever_type(rFlagVariable), 1);
    }

    if (EntityFlags.Is(ELEMENTS)) {
        using data_retriever_type =
            VtuOutputHelperUtilities::VariableFlagRetrieval<ModelPart::ElementType>;
        AddEntityData<bool, data_retriever_type>(
            rFlagName, 1, data_retriever_type(rFlagVariable), 1);
    }
}

template <class TContainerType>
void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    const typename ContainerExpression<TContainerType>::Pointer pContainerExpression,
    const IndexType Precision)
{
    static_assert(std::is_same_v<TContainerType, ModelPart::NodesContainerType> ||
                  std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> ||
                  std::is_same_v<TContainerType, ModelPart::ElementsContainerType>,
                  "Only supports nodal, condition or elemental containers.");

    using data_retriever_type = VtuOutputHelperUtilities::ContainerExpressionRetrieval<TContainerType>;

    using data_retriever_type =
        VtuOutputHelperUtilities::VariableFlagRetrieval<ModelPart::ElementType>;
    AddEntityData<bool, data_retriever_type>(
        rFlagName, 1, data_retriever_type(rFlagVariable), 1);

    if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {


        KRATOS_ERROR_IF_NOT(mpXmlElements) << "No elemental xml element found.";
        const data_retriever_type data_retriever(rExpressionName, pContainerExpression, mpXmlElements->GetEntities(), mpModelPart->GetCommunicator());
        GetCellData()->AddElement(Kratos::make_shared<XmlEntityData<XmlEntities<ModelPart::ElementType>, data_retriever_type>>(mpXmlElements, data_retriever, Precision));
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        KRATOS_ERROR_IF_NOT(mpXmlConditions) << "No condition xml element found.";
        const data_retriever_type data_retriever(rExpressionName, pContainerExpression, mpXmlConditions->GetEntities(), mpModelPart->GetCommunicator());
        GetCellData()->AddElement(Kratos::make_shared<XmlEntityData<XmlEntities<ModelPart::ConditionType>, data_retriever_type>>(mpXmlConditions, data_retriever, Precision));
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>){
        KRATOS_ERROR_IF_NOT(mpXmlPoints) << "No points xml element found.";
        const data_retriever_type data_retriever(rExpressionName, pContainerExpression, mpXmlPoints->GetEntities(), mpModelPart->GetCommunicator());
        GetPointData()->AddElement(Kratos::make_shared<XmlEntityData<XmlNodes, data_retriever_type>>(mpXmlPoints, data_retriever, Precision));
    }
}

void VtuOutput::Write(const std::string& rOutputFileName) const
{
    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileName;
    if (mpModelPart->GetCommunicator().IsDistributed()) {
        output_vtu_file_name << "_" << mpModelPart->GetCommunicator().MyPID();
    }
    output_vtu_file_name << ".vtu";

    // write mesh and data
    std::ofstream output_file;
    output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
    mpXmlVtkFile->Write(output_file);
    output_file.close();

    // create the pvtu file
    // this should be written by the rank zero only.
    const auto& r_data_communicator = mpModelPart->GetCommunicator().GetDataCommunicator();

    std::stringstream list_of_file_names;
    list_of_file_names << output_vtu_file_name.str() << "\n";
    if (r_data_communicator.Rank() == 0) {
        for (int rank = 1; rank < r_data_communicator.Size(); ++rank) {
            std::string msg;
            r_data_communicator.Recv(msg, rank);
            list_of_file_names << msg;
        }
    } else {
        r_data_communicator.Send(list_of_file_names.str(), 0);
    }
    r_data_communicator.Barrier();


    if (r_data_communicator.Rank() == 0) {
        const auto& r_file_names = StringUtilities::SplitStringByDelimiter(list_of_file_names.str(), '\n');

        // create PUnstructuredGrid
        XmlElement unstructured_grid("PUnstructuredGrid");
        unstructured_grid.AddAttribute("GhostLevel", "0");

        // create PPoints
        XmlElement ppoints("PPoints");
        VtuOutputHelperUtilities::CreatePDataArrays(ppoints, *mpXmlPoints);
        unstructured_grid.AddElement(Kratos::make_shared<XmlElement>(ppoints));

        // create PCells
        XmlElement pcells("PCells");
        if (mpXmlElements) {
            VtuOutputHelperUtilities::CreatePDataArrays(pcells, *mpXmlElements);
            unstructured_grid.AddElement(Kratos::make_shared<XmlElement>(pcells));
        } else if (mpXmlConditions) {
            VtuOutputHelperUtilities::CreatePDataArrays(pcells, *mpXmlConditions);
            unstructured_grid.AddElement(Kratos::make_shared<XmlElement>(pcells));
        }

        // create point data
        XmlElement ppoint_data("PPointData");
        if (mpXmlPointData) {
            VtuOutputHelperUtilities::CreatePDataArrays(ppoint_data, *mpXmlPointData);
            unstructured_grid.AddElement(Kratos::make_shared<XmlElement>(ppoint_data));
        }

        // create cell data
        XmlElement pcell_data("PCellData");
        if (mpXmlCellData) {
            VtuOutputHelperUtilities::CreatePDataArrays(pcell_data, *mpXmlCellData);
            unstructured_grid.AddElement(Kratos::make_shared<XmlElement>(pcell_data));
        }

        // now add pieces
        for (const auto& r_file_name : r_file_names) {
            auto piece = Kratos::make_shared<XmlElement>("Piece");
            piece->AddAttribute("Source", r_file_name);
            unstructured_grid.AddElement(piece);
        }

        // create vtk_file
        XmlElement vtk_file("VTKFile");
        vtk_file.AddAttribute("type", "PUnstructuredGrid");
        vtk_file.AddAttribute("version", "0.1");
        vtk_file.AddAttribute("byte_order", "BigEndian");
        vtk_file.AddElement(Kratos::make_shared<XmlElement>(unstructured_grid));

        // writing to file
        std::stringstream output_pvtu_file_name;
        output_pvtu_file_name << rOutputFileName << ".pvtu";

        std::ofstream output_file;
        output_file.open(output_pvtu_file_name.str(), std::ios::out | std::ios::trunc);

        vtk_file.Write(output_file);

        output_file.close();
    }

}

// template instantiations
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<int>&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<double>&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 3>>&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 4>>&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 6>>&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 9>>&, const IndexType);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<int>&, const Flags&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<double>&, const Flags&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 3>>&, const Flags&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 4>>&, const Flags&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 6>>&, const Flags&, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 9>>&, const Flags&, const IndexType);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>(const std::string&, const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer, const IndexType);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer, const IndexType);

} // namespace Kratos