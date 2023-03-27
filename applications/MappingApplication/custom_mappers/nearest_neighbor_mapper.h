//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"
#include "searching/search_interface_info/nearest_neighbor_interface_info.h"
#include "searching/search_local_system/nearest_neighbor_local_system.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Nearest Neighbor Mapper
/** This class implements the Nearest Neighbor Mapping technique.
* Each node on the destination side gets assigned is's closest neighbor on the other side of the interface.
* In the mapping phase every node gets assigned the value of it's neighbor
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) NearestNeighborMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::SearchInterfaceInfoUniquePointerType SearchInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination)
                          : BaseType(rModelPartOrigin, rModelPartDestination) {}

    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
                          : BaseType(rModelPartOrigin,
                                     rModelPartDestination,
                                     JsonParameters)
    {
        KRATOS_TRY;

        auto check_has_nodes = [](const ModelPart& rModelPart){
            if (rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) {
                KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfNodes() == 0) << "No nodes exist in ModelPart \"" << rModelPart.FullName() << "\"" << std::endl;
            }
        };
        check_has_nodes(rModelPartOrigin);
        check_has_nodes(rModelPartDestination);

        this->ValidateInput();
        this->Initialize();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~NearestNeighborMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<NearestNeighborMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{

    int AreMeshesConforming() const override
    {
        KRATOS_WARNING_ONCE("Mapper") << "Developer-warning: \"AreMeshesConforming\" is deprecated and will be removed in the future" << std::endl;
        return BaseType::mMeshesAreConforming;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NearestNeighborMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    ///@name Private Operations
    ///@{

    void CreateSearchLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateSearchLocalSystemsFromNodes(
            NearestNeighborLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    SearchInterfaceInfoUniquePointerType GetSearchInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestNeighborInterfaceInfo>();
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class NearestNeighborMapper

///@} addtogroup block
}  // namespace Kratos.