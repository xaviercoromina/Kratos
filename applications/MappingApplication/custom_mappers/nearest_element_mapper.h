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
#include "searching/search_local_system/nearest_element_local_system.h"
#include "searching/search_interface_info/nearest_element_interface_info.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Interpolative Mapper
/** This class implements the Nearest Element Mapping technique.
* Each node on the destination side gets assigned is's closest condition or element (distance to center)
* on the other side of the interface.
* In the mapping phase every node gets assigned the interpolated value of the condition/element.
* The interpolation is done with the shape funcitons
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) NearestElementMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::SearchInterfaceInfoUniquePointerType SearchInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : BaseType(rModelPartOrigin,rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : BaseType(rModelPartOrigin,
                                    rModelPartDestination,
                                    JsonParameters)
    {
        KRATOS_TRY;

        this->ValidateInput();

        mLocalCoordTol = JsonParameters["local_coord_tolerance"].GetDouble();
        KRATOS_ERROR_IF(mLocalCoordTol < 0.0) << "The local-coord-tolerance cannot be negative" << std::endl;

        this->Initialize();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~NearestElementMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NearestElementMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:
    ///@name Member Variables
    ///@{

    double mLocalCoordTol;

    ///@}

    ///@name Private Operations
    ///@{

    void CreateSearchLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<SearchLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateSearchLocalSystemsFromNodes(
            NearestElementLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    SearchInterfaceInfoUniquePointerType GetSearchInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestElementInterfaceInfo>(mLocalCoordTol);
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "local_coord_tolerance"        : 0.25,
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class NearestElementMapper

}  // namespace Kratos.
