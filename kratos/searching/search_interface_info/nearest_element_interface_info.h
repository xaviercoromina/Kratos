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
#include "utilities/projection_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) NearestElementInterfaceInfo
    : public SearchInterfaceInfo
{
public:

    /// Default constructor.
    explicit NearestElementInterfaceInfo(const double LocalCoordTol=0.0) : mLocalCoordTol(LocalCoordTol) {}

    /**
     * @brief Constructor with coordinates
     * @param rCoordinates The coordinates of the interface object
     * @param SourceLocalSystemIndex The local system index of the source
     * @param SourceRank The rank of the source
     * @param LocalCoordTol The local coordinate tolerance
     */
    explicit NearestElementInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                const IndexType SourceLocalSystemIndex,
                                const IndexType SourceRank,
                                         const double LocalCoordTol=0.0)
        : SearchInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank), mLocalCoordTol(LocalCoordTol) {}

    /**
     * @brief This method is used to create a new interface object
     * @return The new interface object pointer
     */
    SearchInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(mLocalCoordTol);
    }

    /**
     * @brief This method is used to create a new interface object
     * @param rCoordinates The coordinates of the interface object
     * @param SourceLocalSystemIndex The local system index of the source
     * @param SourceRank The rank of the source 
     * @return The new interface object pointer
     */
    SearchInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mLocalCoordTol);
    }

    /**
     * @brief This method is used to get the type of the interface object
     * @return InterfaceObject::ConstructionType The type of the interface object
     */
    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    /**
     * @brief This method is used to process the search result
     * @param rInterfaceObject The interface object
     */
    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    /**
     * @brief This method is used to process the search result for approximation
     * @param rInterfaceObject The interface object
     */
    void ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject) override;

    /**
     * @brief Returns the node ids
     * @param rValue The node ids
     * @param ValueType The type of the value
     */
    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNodeIds;
    }

    /**
     * @brief Returns the shape function values
     * @param rValue The shape function values
     * @param ValueType The type of the value
     */
    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mShapeFunctionValues;
    }

    /**
     * @brief Returns the distance of the closest projection
     * @param rValue The distance of the closest projection
     * @param ValueType The type of the value
     */
    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mClosestProjectionDistance;
    }

    /**
     * @brief Returns the pairing index
     * @param rValue The pairing index
     * @param ValueType The type of the value
     */
    void GetValue(int& rValue,
                  const InfoType ValueType) const override
    {
        rValue = (int)mPairingIndex;
    }

    /**
     * @brief Returns the number of search results  
     * @return std::size_t The number of search results
     */
    std::size_t GetNumSearchResults() const { return mNumSearchResults; }

private:
    ///@name Member Variables
    ///@{

    std::vector<int> mNodeIds;                                                                        /// The node ids of the element
    std::vector<double> mShapeFunctionValues;                                                         /// The shape function values of the element
    double mClosestProjectionDistance = std::numeric_limits<double>::max();                           /// The distance of the closest projection
    ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified; /// The pairing index of the element
    double mLocalCoordTol;                                                                            /// Local coordinate tolerance for the search. NOTE: This is not needed after searching, hence no need to serialize it
    std::size_t mNumSearchResults = 0;                                                                /// The number of search results

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Saves the search result in the member variables
     * @param rInterfaceObject The interface object containing the search result
     * @param ComputeApproximation If true, the search result is saved for the approximation
     */
    void SaveSearchResult(const InterfaceObject& rInterfaceObject,
                          const bool ComputeApproximation);

    ///@}
    ///@name Un accessible methods
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SearchInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("SFValues", mShapeFunctionValues);
        rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
        rSerializer.save("PairingIndex", (int)mPairingIndex);
        rSerializer.save("NumSearchResults", mNumSearchResults);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SearchInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("SFValues", mShapeFunctionValues);
        rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
        int temp;
        rSerializer.load("PairingIndex", temp);
        mPairingIndex = (ProjectionUtilities::PairingIndex)temp;
        rSerializer.load("NumSearchResults", mNumSearchResults);
    }

    ///@}
};

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block
}  // namespace Kratos.