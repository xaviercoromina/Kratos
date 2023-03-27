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
#include <vector>
#include <array>

// External includes

// Project includes
#include "includes/define.h"
#include "searching/search_local_system/search_local_system.h"

namespace Kratos {
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
// forward declaring ModelPart to be avoid including heavy header here
class ModelPart;

/**
 * @namespace SearchUtilities
 * @ingroup KratosCore
 * @brief This class implement some utilities for searching
 * @todo Use OBB to improve it
 * @details Basically bounding boxes
 * @author Philipp Bucher, Jordi Cotela
 */
class KRATOS_API(KRATOS_CORE) SearchUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the bounding box array
    using BoundingBoxType = std::array<double, 6>;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    using SearchInterfaceInfoUniquePointerType = Kratos::unique_ptr<SearchInterfaceInfo>;

    using SearchInterfaceInfoPointerType = Kratos::shared_ptr<SearchInterfaceInfo>;
    using SearchInterfaceInfoPointerVectorType = std::vector<std::vector<SearchInterfaceInfoPointerType>>;

    using SearchLocalSystemPointer = Kratos::unique_ptr<SearchLocalSystem>;
    using SearchLocalSystemPointerVector = std::vector<SearchLocalSystemPointer>;
    using SearchLocalSystemPointerVectorPointer = Kratos::shared_ptr<SearchLocalSystemPointerVector>;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the bounding box of a given model part (global)
     * @param rModelPart The model part to compute the bounding box
     * @return The bounding box
     */
    static BoundingBoxType ComputeLocalBoundingBox(const ModelPart& rModelPart);

    /**
     * @brief This method computes the bounding box of a given model part (local)
     * @param rModelPart The model part to compute the bounding box
     * @return The bounding box
     */
    static BoundingBoxType ComputeGlobalBoundingBox(const ModelPart& rModelPart);

    /**
     * @brief This method computes several bounding boxes with a given tolerance from a set of bounding boxes
     * @param rBoundingBoxes The bounding boxes to compute the bounding boxes with tolerance
     * @param Tolerance The tolerance to compute the bounding boxes with tolerance
     * @param rBoundingBoxesWithTolerance The bounding boxes with tolerance
     */
    static void ComputeBoundingBoxesWithTolerance(
        const std::vector<double>& rBoundingBoxes,
        const double Tolerance,
        std::vector<double>& rBoundingBoxesWithTolerance
        );

    /**
     * @brief This method returns the stream of a bounding box
     * @param rBoundingBox The bounding box to compute the stream
     * @return The stream of the bounding box
     */
    static std::string BoundingBoxStringStream(const BoundingBoxType& rBoundingBox);

    /**
     * @brief This method checks if a point is inside a bounding box
     * @param rBoundingBox The bounding box
     * @param rCoords The coordinates of the point
     * @return True if the point is inside the bounding box, false otherwise
     */
    static bool PointIsInsideBoundingBox(
        const BoundingBoxType& rBoundingBox,
        const array_1d<double, 3>& rCoords
        );

    /**
     * @brief This fills the buffer before the local search
     * @param rSearchLocalSystems The local search systems
     * @param rBoundingBoxes The bounding boxes
     * @param BufferSizeEstimate The estimated size of the buffer
     * @param rSendBuffer The buffer to be filled
     * @param rSendSizes The sizes of the serialized SearchLocalSystems
     */
    static void FillBufferBeforeLocalSearch(
        const SearchLocalSystemPointerVector& rSearchLocalSystems,
        const std::vector<double>& rBoundingBoxes,
        const SizeType BufferSizeEstimate,
        std::vector<std::vector<double>>& rSendBuffer,
        std::vector<int>& rSendSizes
        );

    /**
     * @brief This method creates the SearchInterfaceInfos from the received buffer
     * @param rRecvBuffer The buffer containing the serialized SearchInterfaceInfos
     * @param rpRefInterfaceInfo The reference SearchInterfaceInfo
     * @param CommRank The rank of the process
     * @param rSearchInterfaceInfosContainer The container to be filled with the SearchInterfaceInfos
     */
    static void CreateSearchInterfaceInfosFromBuffer(
        const std::vector<std::vector<double>>& rRecvBuffer,
        const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
        const int CommRank,
        SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer
        );

    /**
     * @brief This method fills the buffer with the serialized SearchInterfaceInfos
     * @param rSearchInterfaceInfosContainer The container containing the SearchInterfaceInfos
     * @param rpRefInterfaceInfo The reference SearchInterfaceInfo
     * @param CommRank The rank of the process
     * @param rSendBuffer The buffer to be filled with the serialized SearchInterfaceInfos
     * @param rSendSizes The sizes of the serialized SearchInterfaceInfos
     * @note The SearchInterfaceInfos are serialized in the same order as they are stored in the container
     */
    static void FillBufferAfterLocalSearch(
        SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer,
        const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
        const int CommRank,
        std::vector<std::vector<char>>& rSendBuffer,
        std::vector<int>& rSendSizes
        );

    /**
     * @brief This method assigns the SearchInterfaceInfos to the SearchLocalSystems
     * @param rSearchInterfaceInfosContainer The container containing the SearchInterfaceInfos
     * @param rpSearchLocalSystems The container containing the SearchLocalSystems
     * @note The SearchInterfaceInfos are assigned to the SearchLocalSystems based on the
     *      SearchInterfaceInfo::mInterfaceId
     * @note The SearchInterfaceInfos are assigned to the SearchLocalSystems in the same order
     *     as the SearchLocalSystems are stored in the container
     */
    static void AssignInterfaceInfosAfterRemoteSearch(
        const SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer,
        SearchLocalSystemPointerVectorPointer& rpSearchLocalSystems
        );

    /**
     * @brief This method deserializes the SearchInterfaceInfos
     * @param rSendBuffer The buffer containing the serialized SearchInterfaceInfos
     * @param rpRefInterfaceInfo The reference SearchInterfaceInfo
     * @param CommRank The rank of the process
     * @param rSearchInterfaceInfosContainer The container where the SearchInterfaceInfos will be stored
     */
    static void DeserializeSearchInterfaceInfosFromBuffer(
        const std::vector<std::vector<char>>& rSendBuffer,
        const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
        const int CommRank,
        SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer
        );

    ///@}
};

/**
 * @class SearchInterfaceInfoSerializer
 * @ingroup KratosCore
 * @brief Helper class to serialize/deserialize a vector containing SearchInterfaceInfos
 * @details This class serializes the vector containing the SearchInterfaceInfos (Shared Ptrs)
 * The goal of this class is to have a more efficient/faster implementation than the
 * one of the Serializer by avoiding the casting that is done in the serializer when pointers
 * are serialized
 * @TODO test the performance against the Serializer
 * @author Philipp Bucher
 */
class KRATOS_API(KRATOS_CORE) SearchInterfaceInfoSerializer
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the unique pointer type for SearchInterfaceInfo
    using SearchInterfaceInfoUniquePointerType = Kratos::unique_ptr<SearchInterfaceInfo>;

    /// Definition of the shared pointer type for SearchInterfaceInfo
    using SearchInterfaceInfoPointerType = Kratos::shared_ptr<SearchInterfaceInfo>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rSearchInterfaceInfosContainer The vector containing the SearchInterfaceInfos
     * @param rpRefInterfaceInfo A pointer to the SearchInterfaceInfo to be used for deserialization
     */
    SearchInterfaceInfoSerializer(
        std::vector<SearchInterfaceInfoPointerType>& rSearchInterfaceInfosContainer,
        const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo
        ) : mrInterfaceInfos(rSearchInterfaceInfosContainer)
          , mrpRefInterfaceInfo(rpRefInterfaceInfo->Create())
        {
            // Empty constructor
        }
    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<SearchInterfaceInfoPointerType>& mrInterfaceInfos; /// The vector containing the SearchInterfaceInfos
    SearchInterfaceInfoPointerType mrpRefInterfaceInfo;            /// A pointer to the SearchInterfaceInfo to be used for deserialization

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;
    virtual void load(Serializer& rSerializer);

    ///@}
};

}
