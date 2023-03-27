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
#include "searching/interface_communicator.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object for exchanging data on the Interface in MPI
/** It extends it's baseclass by remote-searching capabilities. I.e. before the local search,
 * data is exchanged among the partitions to be used by the local-search
*/
class KRATOS_API(MAPPING_APPLICATION) InterfaceCommunicatorMPI : public InterfaceCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceCommunicatorMPI
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicatorMPI);

    /// The type of the buffer type for doubles
    using BufferTypeDouble = std::vector<std::vector<double>>;

    /// The type of the buffer type for chars
    using BufferTypeChar = std::vector<std::vector<char>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Interface Communicator MPI object
     * @param rModelPartOrigin The ModelPart from which the InterfaceObjects are constructed
     * @param rSearchLocalSystems The SearchLocalSystems used for the local search
     * @param SearchSettings The settings for the search
     */
    InterfaceCommunicatorMPI(ModelPart& rModelPartOrigin,
                             SearchLocalSystemPointerVector& rSearchLocalSystems,
                             Parameters SearchSettings);

    /// Destructor.
    virtual ~InterfaceCommunicatorMPI()
    {
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceCommunicatorMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceCommunicatorMPI";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief This function initializes the search
     * @param rpRefInterfaceInfo The InterfaceInfo of the reference interface
     */
    void InitializeSearch(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    /**
     * @brief This function constructs the InterfaceObjects on the Destination
     * @details In serial it only does it once, whereas in MPI this involves Data-Exchange!
     *          Imagine a sliding interface, there the partitions might change!
     * @param rpRefInterfaceInfo The InterfaceInfo of the reference interface
     */
    void InitializeSearchIteration(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    /**
     * @brief This function finalizes the search
     * @param rpRefInterfaceInfo The InterfaceInfo of the reference interface
     */
    void FinalizeSearchIteration(const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<double> mGlobalBoundingBoxes;  /// The global bounding boxes of the local search systems xmax, xmin,  ymax, ymin,  zmax, zmin

    int mCommRank;                             /// The rank of the current MPI-process
    int mCommSize;                             /// The size of the MPI-communicator

    std::vector<int> mSendSizes;               /// The sizes of the send buffers
    std::vector<int> mRecvSizes;               /// The sizes of the receive buffers

    BufferTypeDouble mSendBufferDouble;        /// The send buffer for doubles
    BufferTypeDouble mRecvBufferDouble;        /// The receive buffer for doubles

    BufferTypeChar mSendBufferChar;            /// The send buffer for chars
    BufferTypeChar mRecvBufferChar;            /// The receive buffer for chars

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function estimates the buffer size for the data exchange
    */
    std::size_t GetBufferSizeEstimate() const
    {
        return mrSearchLocalSystems.size() / mCommSize;
    }

    /**
     * @brief This function computes the global bounding boxes
     * @details The global bounding boxes are used for the remote-search
    */
    void ComputeGlobalBoundingBoxes();

    /**
     * @brief This function excahnges the data between the partitions
     * @tparam TDataType The type of the data to be exchanged
     * @param rSendBuffer The send buffer
     * @param rRecvBuffer The receive buffer
     */
    template< typename TDataType >
    int ExchangeDataAsync(
        const std::vector<std::vector<TDataType>>& rSendBuffer,
        std::vector<std::vector<TDataType>>& rRecvBuffer);

    ///@}

}; // Class InterfaceCommunicatorMPI

///@}

///@} addtogroup block

}  // namespace Kratos.
