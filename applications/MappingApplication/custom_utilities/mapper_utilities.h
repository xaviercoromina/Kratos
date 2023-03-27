//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED)
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED

// System includes
#include <array>
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "mapping_application_variables.h"
#include "mappers/mapper_flags.h"
#include "searching/search_local_system/search_local_system.h"

namespace Kratos {
namespace MapperUtilities {

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Node<3> NodeType;

static void FillFunction(const NodeType& rNode,
                         const Variable<double>& rVariable,
                         double& rValue)
{
    rValue = rNode.FastGetSolutionStepValue(rVariable);
}

static void FillFunctionNonHist(const NodeType& rNode,
                                const Variable<double>& rVariable,
                                double& rValue)
{
    rValue = rNode.GetValue(rVariable);
}

static inline std::function<void(const NodeType&, const Variable<double>&, double&)>
GetFillFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL))
        return &FillFunctionNonHist;
    return &FillFunction;
}

static void UpdateFunction(NodeType& rNode,
                           const Variable<double>& rVariable,
                           const double Value,
                           const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) = Value * Factor;
}

static void UpdateFunctionWithAdd(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) += Value * Factor;
}

static void UpdateFunctionNonHist(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) = Value * Factor;
}

static void UpdateFunctionNonHistWithAdd(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) += Value * Factor;
}

static inline std::function<void(NodeType&, const Variable<double>&, const double, const double)>
GetUpdateFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES) && rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHistWithAdd;
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
        return &UpdateFunctionWithAdd;
    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHist;
    return &UpdateFunction;
}

template<class TVectorType, bool TParallel=true>
void UpdateSystemVectorFromModelPart(
    TVectorType& rVector,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction(rMappingOptions);

    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    // necessary bcs the Trilinos Vector is not threadsafe in the default configuration
    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    KRATOS_ERROR_IF(!rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL) && !rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Solution step variable \"" << rVariable.Name() << "\" missing in ModelPart \"" << rModelPart.FullName() << "\"!" << std::endl;

    IndexPartition<std::size_t>(num_local_nodes, num_threads).for_each([&](const std::size_t i){
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
    });

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateModelPartFromSystemVector(
    const TVectorType& rVector,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;

    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    // necessary bcs the Trilinos Vector is not threadsafe in the default configuration
    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    KRATOS_ERROR_IF(!rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL) && !rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Solution step variable \"" << rVariable.Name() << "\" missing in ModelPart \"" << rModelPart.FullName() << "\"!" << std::endl;

    IndexPartition<std::size_t>(num_local_nodes, num_threads).for_each([&](const std::size_t i){
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
    });

    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL)) {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    KRATOS_CATCH("");
}

/**
* @brief Assigning INTERFACE_EQUATION_IDs to the nodes, with and without MPI
* This function assigns the INTERFACE_EQUATION_IDs to the nodes, which
* act as EquationIds for the MappingMatrix. This work with and without MPI,
* in MPI a ScanSum is performed with the local number of nodes
* @param rModelPartCommunicator The Modelpart-Communicator to be used
* @author Philipp Bucher
*/
void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

void CreateSearchLocalSystemsFromNodes(const SearchLocalSystem& rSearchLocalSystemPrototype,
                                       const Communicator& rModelPartCommunicator,
                                       std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems);

void CreateSearchLocalSystemsFromGeometries(const SearchLocalSystem& rSearchLocalSystemPrototype,
                                            const Communicator& rModelPartCommunicator,
                                            std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems);

template <class T1, class T2, class T3>
bool PointsAreCollinear(
    const T1& rP1,
    const T2& rP2,
    const T3& rP3)
{
    // TODO this can probably be optimized
    const double a = MathUtils<double>::Norm3(rP1-rP2);
    const double b = MathUtils<double>::Norm3(rP2-rP3);
    const double c = MathUtils<double>::Norm3(rP3-rP1);

    const double s = (a+b+c) / 2.0;

    return (std::sqrt(s*(s-a)*(s-b)*(s-c))) < 1e-12;
}

void CheckInterfaceModelParts(const int CommRank);

void KRATOS_API(MAPPING_APPLICATION) SaveCurrentConfiguration(ModelPart& rModelPart);
void KRATOS_API(MAPPING_APPLICATION) RestoreCurrentConfiguration(ModelPart& rModelPart);

template<class TDataType>
void EraseNodalVariable(ModelPart& rModelPart, const Variable<TDataType>& rVariable)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.GetData().Erase(rVariable);
    });

    KRATOS_CATCH("");
}

}  // namespace MapperUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined
