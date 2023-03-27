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
#include "searching/search_interface_info/search_interface_info.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// This is the "Condition" of the search
/** This class assembles the local system for the searches, using the Information that
 * was provided by the "SearchInterfaceInfo"
*/
class SearchLocalSystem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SearchLocalSystem
    KRATOS_CLASS_POINTER_DEFINITION(SearchLocalSystem);

    using SearchInterfaceInfoPointerType = Kratos::shared_ptr<SearchInterfaceInfo>;
    using SearchLocalSystemUniquePointer = Kratos::unique_ptr<SearchLocalSystem>;

    using CoordinatesArrayType = typename SearchInterfaceInfo::CoordinatesArrayType;

    using MatrixType = Matrix;
    using EquationIdVectorType = std::vector<int>; // int bcs of mpi

    using NodePointerType = InterfaceObject::NodePointerType;
    using GeometryPointerType = InterfaceObject::GeometryPointerType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class PairingStatus
    {
        NoInterfaceInfo,
        Approximation,
        InterfaceInfoFound
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~SearchLocalSystem() = default;

    ///@}
    ///@name Operations
    ///@{

    void EquationIdVectors(EquationIdVectorType& rOriginIds,
                           EquationIdVectorType& rDestinationIds)
    {
        if (!mIsComputed) {
            CalculateAll(mLocalMappingMatrix, mOriginIds, mDestinationIds, mPairingStatus);
            mIsComputed = true;
        }

        rOriginIds      = mOriginIds;
        rDestinationIds = mDestinationIds;
    }

    void CalculateLocalSystem(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds) const
    {
        if (mIsComputed) {
            // This will be called if the EquationIdVectors have been querried before
            // i.e. matrix-based mapping
            rLocalMappingMatrix = mLocalMappingMatrix;
            rOriginIds      = mOriginIds;
            rDestinationIds = mDestinationIds;
        }
        else {
            // This will be called if the EquationIdVectors have NOT been querried before
            // i.e. matrix-free mapping
            CalculateAll(rLocalMappingMatrix, rOriginIds, rDestinationIds, mPairingStatus);
        }
    }

    /**
    * @brief Resizing the output if no InterfaceInfo is available
    * This function resizes the system vectors to zero and also sets that no valid
    * Information from the other side could be found to compute the local system
    * @param rLocalMappingMatrix The vector conatining the mapping weights
    * @param rOriginIds The vector containing the ids on the origin
    * @param rDestinationIds The vector containing the ids on the destination
    * @param rPairingStatus The pairingstatus of the SearchLocalSystem
    * @see CalculateAll
    * @author Philipp Bucher
    */
    void ResizeToZero(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      SearchLocalSystem::PairingStatus& rPairingStatus) const
    {
        rPairingStatus = SearchLocalSystem::PairingStatus::NoInterfaceInfo;

        rLocalMappingMatrix.resize(0, 0, false);
        rOriginIds.resize(0);
        rDestinationIds.resize(0);
    }

    virtual CoordinatesArrayType& Coordinates() const = 0;


    void AddInterfaceInfo(SearchInterfaceInfoPointerType pInterfaceInfo) // TODO pass by const ref?
    {
        mInterfaceInfos.push_back(pInterfaceInfo);
    }

    bool HasInterfaceInfo() const
    {
        return mInterfaceInfos.size() > 0;
    }

    bool HasInterfaceInfoThatIsNotAnApproximation() const
    {
        for (const auto& r_info : mInterfaceInfos) {
            if (!r_info->GetIsApproximation()) {
                return true;
            }
        }
        return false;
    }

    virtual bool IsDoneSearching() const
    {
        return HasInterfaceInfoThatIsNotAnApproximation();
    }

    virtual SearchLocalSystemUniquePointer Create(NodePointerType pNode) const
    {
        KRATOS_ERROR << "Create is not implemented for NodePointerType!" << std::endl;
    }

    virtual SearchLocalSystemUniquePointer Create(GeometryPointerType pGeometry) const
    {
        KRATOS_ERROR << "Create is not implemented for GeometryPointerType!" << std::endl;
    }

    virtual void Clear()
    {
        mInterfaceInfos.clear();
        mLocalMappingMatrix.clear();
        mOriginIds.clear();
        mDestinationIds.clear();
        mIsComputed = false;
    }

    PairingStatus GetPairingStatus() const
    {
        return mPairingStatus;
    }

    virtual void SetPairingStatusForPrinting()
    {
        KRATOS_ERROR << "SetPairingStatusForPrinting is not implemented!" << std::endl;
    }

    ///@}
    ///@name Input and output
    ///@{

    virtual void PairingInfo(std::ostream& rOStream, const int EchoLevel) const = 0;

    /// Turn back information as a string.
    virtual std::string Info() const {return "SearchLocalSystem";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected Life Cycle
    ///@{

    SearchLocalSystem() = default; // only accessbíble by derived classes

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<SearchInterfaceInfoPointerType> mInterfaceInfos;

    bool mIsComputed = false;

    MatrixType mLocalMappingMatrix;
    EquationIdVectorType mOriginIds;
    EquationIdVectorType mDestinationIds;

    mutable PairingStatus mPairingStatus = PairingStatus::NoInterfaceInfo;

    ///@}
    ///@name Protected Operations
    ///@{

    // This function calculates the components necessary for the mapping
    // Note that it is "const", therefore it can NOT modify its members
    // Whether members are to be saved is decided in other functions of this class
    virtual void CalculateAll(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds,
                              SearchLocalSystem::PairingStatus& rPairingStatus) const = 0;

    ///@}

}; // Class SearchLocalSystem

///@}

///@} addtogroup block

}  // namespace Kratos.
