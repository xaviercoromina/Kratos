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
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "mappers/mapper_flags.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "utilities/builtin_timer.h"
#include "searching/configures/interface_object_configure.h"
#include "searching/search_local_system/search_local_system.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Object for exchanging data on the Interface
/** Mapping requires knowledge about the "other" side of the Interface. This class communicates the
 * data required by the mappers, hence it also includes the (local) searching
*/
class KRATOS_API(MAPPING_APPLICATION) InterfaceCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicator);

    /// The type of the InterfaceObject
    using SearchInterfaceInfoUniquePointerType = Kratos::unique_ptr<SearchInterfaceInfo>;

    /// The type of the InterfaceObject pointer
    using SearchInterfaceInfoPointerType = Kratos::shared_ptr<SearchInterfaceInfo>;

    /// The type of the InterfaceObject pointer vector
    using SearchInterfaceInfoPointerVectorType = std::vector<std::vector<SearchInterfaceInfoPointerType>>;

    /// The type of the SearchLocalSystem pointer
    using SearchLocalSystemPointer = Kratos::unique_ptr<SearchLocalSystem>;

    /// The type of the SearchLocalSystem pointer vector
    using SearchLocalSystemPointerVector = std::vector<SearchLocalSystemPointer>;

    /// The type of the BinsObjectDynamic
    using BinsUniquePointerType = Kratos::unique_ptr<BinsObjectDynamic<InterfaceObjectConfigure>>;

    /// The type of the InterfaceObjectContainer
    using InterfaceObjectContainerType = InterfaceObjectConfigure::ContainerType;

    /// The type of the InterfaceObjectContainer pointer
    using InterfaceObjectContainerUniquePointerType = Kratos::unique_ptr<InterfaceObjectContainerType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Interface Communicator object
     * @param rModelPartOrigin The ModelPart from which the Interface is taken
     * @param rSearchLocalSystems The SearchLocalSystems for the Interface
     * @param SearchSettings The settings for the search
     */
    InterfaceCommunicator(ModelPart& rModelPartOrigin,
                          SearchLocalSystemPointerVector& rSearchLocalSystems,
                          Parameters SearchSettings);

    /// Destructor.
    virtual ~InterfaceCommunicator() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function excahnges the interface data
     * @param rComm The Communicator
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    void ExchangeInterfaceData(const Communicator& rComm,
                               const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This function checks if the meshes are conforming
     * @return int 1 if conforming, 0 if not
     */
    int AreMeshesConforming() {
        return mMeshesAreConforming;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "InterfaceCommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPartOrigin;                                        /// The ModelPart from which the Interface is taken
    const SearchLocalSystemPointerVector& mrSearchLocalSystems;          /// The SearchLocalSystems for the Interface
    SearchInterfaceInfoPointerVectorType mSearchInterfaceInfosContainer; /// This contains the InterfaceInfos for all ranks! => needed to do the async communication

    BinsUniquePointerType mpLocalBinStructure;                           /// The local bin structure

    InterfaceObjectContainerUniquePointerType mpInterfaceObjectsOrigin;  /// The InterfaceObjects of the origin

    Parameters mSearchSettings;                                          /// The settings for the search
    double mSearchRadius = -1.0;                                         /// The search radius

    int mEchoLevel = 0;                                                  /// Echo level for the search
    int mMeshesAreConforming = 0;                                        /// 1 if conforming, 0 if not

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This function initializes the search
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    virtual void InitializeSearch(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    /**
     * @brief This function finalizes the search
     */
    virtual void FinalizeSearch();

    /**
     * @brief This function initializes the search iteration
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    virtual void InitializeSearchIteration(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    /**
     * @brief This function finalizes the search iteration
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    virtual void FinalizeSearchIteration(const SearchInterfaceInfoUniquePointerType& rpInterfaceInfo);

    /**
     * @brief This function filters the InterfaceInfos for the successful search
     */
    void FilterInterfaceInfosSuccessfulSearch();

    /**
     * @brief This function assigns the InterfaceInfos
     */
    void AssignInterfaceInfos();

    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
     * @brief This function computes the max edge length of the entities in the container
     * @tparam TContainer The type of the container
     * @param rEntityContainer The container
     * @return double The max edge length
     */
    template <typename TContainer>
    static double ComputeMaxEdgeLengthLocal(const TContainer& rEntityContainer);

    /**
     * @brief This function computes the max edge length of the entities in the container
     * @tparam TContainer The type of the container
     * @param rEntityContainer The container
     * @return double The max edge length
     */
    static double ComputeSearchRadius(const ModelPart& rModelPart, int EchoLevel);

    /**
     * @brief This function computes the max edge length of the entities in the container
     * @tparam TContainer The type of the container
     * @param rEntityContainer The container
     * @return double The max edge length
     */
    static double ComputeSearchRadius(const ModelPart& rModelPart1, const ModelPart& rModelPart2, const int EchoLevel);

    /**
     * @brief This function conducts the local search
     */
    void ConductLocalSearch();

    /**
     * @brief This function creates the InterfaceObjects for the origin
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    void CreateInterfaceObjectsOrigin(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    /**
     * @brief This function updates the InterfaceObjects for the origin
     */
    void UpdateInterfaceObjectsOrigin();

    /**
     * @brief This function initializes the bins search structure
     */
    void InitializeBinsSearchStructure();

    /**
     * @brief This function performs the search and the exchange of the data on the interface
     * @param rpRefInterfaceInfo The reference InterfaceInfo
     */
    void ConductSearchIteration(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    /**
     * @brief This function returns if all neighbors have been found
     * @param rComm The Communicator
     * @return true All neighbors have been found
     */
    bool AllNeighborsFound(const Communicator& rComm) const;

    /**
     * @brief This function prints the info about the current search
     * @param rComm The Communicator
     * @param rTimer The timer
     */
    void PrintInfoAboutCurrentSearchSuccess(
        const Communicator& rComm,
        const BuiltinTimer& rTimer) const;

    ///@}

}; // Class InterfaceCommunicator

///@}

///@}

///@} addtogroup block

}  // namespace Kratos.
