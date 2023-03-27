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
#include "includes/define.h"
#include "searching/interface_object.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Object for storing data that is needed to construct the local-info-system, for example the local-mapper-system
 * @details When constructing th local-mapper-system, some data from the "other" side of the
 * Interface is needed which is different for every sharing info/mapper.
 * This class stores the data needed and provides it later when the local-system is
 * constructed.
 * For the data-exchange in MPI the data is saved using serialization
 * @author Philipp Bucher, Jordi Cotela
 */
class SearchInterfaceInfo
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SearchInterfaceInfo
    KRATOS_CLASS_POINTER_DEFINITION(SearchInterfaceInfo);

    /// The type of the index
    using IndexType = std::size_t;

    /// The type of the coordinates
    using CoordinatesArrayType = typename InterfaceObject::CoordinatesArrayType;

    /// The type of the node
    using NodeType = InterfaceObject::NodeType;

    /// The type of the geometry
    using GeometryType = InterfaceObject::GeometryType;

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief Enum for the different types of SearchInterfaceInfo
     * @details This enum is used to identify the type of the SearchInterfaceInfo
     * @note This is needed for the serialization
     */
    enum class InfoType
    {
        Dummy
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SearchInterfaceInfo() = default;

    explicit SearchInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                        const IndexType SourceLocalSystemIndex,
                        const IndexType SourceRank)
        : mSourceLocalSystemIndex(SourceLocalSystemIndex),
          mCoordinates(rCoordinates),
          mSourceRank(SourceRank)
    {}

    /// Destructor.
    virtual ~SearchInterfaceInfo() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Processing the result of the search
    * This function processes the results of the search, i.e. it extracts
    * the information needed from the neighbor and saves it such that it can be
    * later accesse by the SearchLocalSystem for assembling it's local system.
    * This happens in the remote partition.
    * @param rpInterfaceObject The InterfaceObject found by the search
    * @author Philipp Bucher
    */
    virtual void ProcessSearchResult(const InterfaceObject& rInterfaceObject) = 0;

    /**
    * @brief Processing the result of the search for computing an approximation
    * This function processes the results of the search for computing an approximation.
    * This can be necessary if e.g. a projection fails.
    * In case an approximation is found the function "SetIsApproximation" has to be
    * called to set internal flags for data-exchange and the computation of the
    * local system by the SearchLocalSystem
    * This happens in the remote partition.
    * It's implementation is optional
    * @param rpInterfaceObject The InterfaceObject found by the search
    * @see ProcessSearchResult
    * @see SetIsApproximation
    * @author Philipp Bucher
    */
    virtual void ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject) {}

    virtual SearchInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                                const IndexType SourceLocalSystemIndex,
                                                const IndexType SourceRank) const = 0;

    // needed for serialization
    virtual SearchInterfaceInfo::Pointer Create() const = 0;

    /**
    * @brief returning the type of construction for the InterfaceObject
    * The returned type is used to create the object for the bin-search
    * It is tightly coupled to this class since it decides which information
    * can be extracted by this class
    * @author Philipp Bucher
    */
    virtual InterfaceObject::ConstructionType GetInterfaceObjectType() const = 0;

    /**
     * @brief This function computes the distance from another SearchInterfaceInfo
     * @param rOther The other SearchInterfaceInfo
     * @return The distance
     */
    double ComputeDistance(const SearchInterfaceInfo& rOther)
    {
        return std::sqrt( std::pow(mCoordinates[0] - rOther.mCoordinates[0] , 2) +
                          std::pow(mCoordinates[1] - rOther.mCoordinates[1] , 2) +
                          std::pow(mCoordinates[2] - rOther.mCoordinates[2] , 2) );
    }

    /**
     * @brief This function computes the distance from a set of coordinates
     * @param rOtherCoordinates The coordinates
     * @return The distance
     */
    double ComputeDistance(const CoordinatesArrayType& rOtherCoordinates)
    {
        return std::sqrt( std::pow(mCoordinates[0] - rOtherCoordinates[0] , 2) +
                          std::pow(mCoordinates[1] - rOtherCoordinates[1] , 2) +
                          std::pow(mCoordinates[2] - rOtherCoordinates[2] , 2) );
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This function is used to get the local system index
     * @return The local system index of the source
     */
    IndexType GetLocalSystemIndex() const { return mSourceLocalSystemIndex; }

    /**
     * @brief This function is used to get the rank of the source
     * @return The rank of the source
     */
    IndexType GetSourceRank() const { return mSourceRank; }

    /**
     * @brief This function is used to get the flag if the search was successful
     * @return The flag if the search was successful
     */
    bool GetLocalSearchWasSuccessful() const { return mLocalSearchWasSuccessful; }

    /**
     * @brief This function is used to get if an approximation is considered
     * @return The flag if an approximation is considered
     */
    bool GetIsApproximation() const { return mIsApproximation; }

    /**
     * @brief This function is to get the coordinates of the source
     * @return The coordinates of the source
     */
    CoordinatesArrayType& Coordinates()
    {
        return mCoordinates;
    }

    virtual void GetValue(int& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::size_t& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(double& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(bool& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(GeometryType& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }

    virtual void GetValue(std::vector<int>& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<std::size_t>& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<double>& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<bool>& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<GeometryType>& rValue, const InfoType ValueType) const { KRATOS_ERROR << "Base class function called!" << std::endl; }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "SearchInterfaceInfo";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    /* These variables need serialization */
    IndexType mSourceLocalSystemIndex;  /// The local system index of the source

    /* These variables are NOT being serialized bcs they are not needed after searching! */
    CoordinatesArrayType mCoordinates;  /// The coordinates of the source
    IndexType mSourceRank = 0;          /// The rank of the source

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This function is used to set the flag if the local search was successful
     */
    void SetLocalSearchWasSuccessful()
    {
        mLocalSearchWasSuccessful = true;
        mIsApproximation = false;
    }

    /**
     * @brief This function is used to set the flag if an approximation is considered
     */
    void SetIsApproximation()
    {
        // If an approximation is found also means that the local search has been successful!
        // this is needed otherwise it won't be properly processes by the search
        // the SearchLocalSystem has to take care of this!
        mLocalSearchWasSuccessful = true;

        mIsApproximation = true;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    bool mIsApproximation = false;           /// Flag if an approximation is considered

    bool mLocalSearchWasSuccessful = false;  /// Flag if the local search was successful. NOTE: This is not being serialized since it is not needed after mpi-data-exchange!

    ///@}
    ///@name Private  Access
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("LocalSysIdx", mSourceLocalSystemIndex);
        rSerializer.save("IsApproximation", mIsApproximation);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("LocalSysIdx", mSourceLocalSystemIndex);
        rSerializer.load("IsApproximation", mIsApproximation);
    }

    ///@}

}; // Class SearchInterfaceInfo

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.
