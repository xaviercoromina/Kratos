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
#include "searching/search_interface_info/search_interface_info.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class NearestNeighborInterfaceInfo
 * @ingroup KratosCore
 * @brief This class is used to store the information of the nearest neighbor
 * @details This class is used to store the information of the nearest neighbor
 * @author Philipp Bucher
 */
class KRATOS_API(KRATOS_CORE) NearestNeighborInterfaceInfo
    : public SearchInterfaceInfo
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestNeighborInterfaceInfo
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborInterfaceInfo);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NearestNeighborInterfaceInfo() {}

    /**
     * @brief Constructor with coordinates
     * @param rCoordinates The coordinates of the interface object
     * @param SourceLocalSystemIndex The local system index of the source
     * @param SourceRank The rank of the source
    */
    explicit NearestNeighborInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                 const IndexType SourceLocalSystemIndex,
                                 const IndexType SourceRank)
        : SearchInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates a new interface object
     * @return The new interface object pointer
     */
    SearchInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>();
    }

    /**
     * @brief This method creates a new interface object
     * @param rCoordinates The coordinates of the interface object
     * @param SourceLocalSystemIndex The local system index of the source
     * @param SourceRank The rank of the source
     * @return The new interface object pointer
     */
    SearchInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    /**
     * @brief This method returns the type of the interface object
     * @return The type of the interface object
     */
    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Node_Coords;
    }

    /**
     * @brief This method processes the search result
     * @param rInterfaceObject The interface object that was found
     */
    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    /**
     * @brief This method returns the nearest neighbor Id
     * @param rValue Then nearest neighbor Id
     * @param ValueType The type of the requested variable
     * @return The nearest neighbor Id
     */
    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborId;
    }

    /**
     * @brief This method returns the value  nearest neighbor distance
     * @param rValue The nearest neighbor distance
     * @param ValueType The type of the requested variable
     * @return The nearest neighbor distance
     */
    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborDistance;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<int> mNearestNeighborId = {};                             /// The Id of the nearest neighbor
    double mNearestNeighborDistance = std::numeric_limits<double>::max(); /// The distance to the nearest neighbor

    ///@}
    ///@name Un accessible methods
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SearchInterfaceInfo );
        rSerializer.save("NearestNeighborId", mNearestNeighborId);
        rSerializer.save("NearestNeighborDistance", mNearestNeighborDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SearchInterfaceInfo );
        rSerializer.load("NearestNeighborId", mNearestNeighborId);
        rSerializer.load("NearestNeighborDistance", mNearestNeighborDistance);
    }

    ///@}
};

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block
}  // namespace Kratos.