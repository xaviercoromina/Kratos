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

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "searching/search_local_system/nearest_neighbor_local_system.h"
#include "searching/search_interface_info/nearest_neighbor_interface_info.h"

namespace Kratos
{
///@name Kratos Classes
///@{

void NearestNeighborLocalSystem::CalculateAll(
    MatrixType& rLocalMappingMatrix,
    EquationIdVectorType& rOriginIds,
    EquationIdVectorType& rDestinationIds,
    SearchLocalSystem::PairingStatus& rPairingStatus
    ) const
{
    if (mInterfaceInfos.size() > 0) {
        rPairingStatus = SearchLocalSystem::PairingStatus::InterfaceInfoFound;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

        std::vector<int> nearest_neighbor_id;
        double nearest_neighbor_distance;
        mInterfaceInfos[0]->GetValue(nearest_neighbor_id, SearchInterfaceInfo::InfoType::Dummy);
        mInterfaceInfos[0]->GetValue(nearest_neighbor_distance, SearchInterfaceInfo::InfoType::Dummy);
        rOriginIds = nearest_neighbor_id;

        for (std::size_t i=1; i<mInterfaceInfos.size(); ++i) {
            // no check if this InterfaceInfo is an approximation is necessary
            // bcs this does not exist for NearestNeighbor
            double distance;
            mInterfaceInfos[i]->GetValue(distance, SearchInterfaceInfo::InfoType::Dummy);

            if (distance < nearest_neighbor_distance) {
                nearest_neighbor_distance = distance;
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id, SearchInterfaceInfo::InfoType::Dummy);
                rOriginIds = nearest_neighbor_id;
            } else if (distance == nearest_neighbor_distance) {
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id, SearchInterfaceInfo::InfoType::Dummy);
                rOriginIds.insert(rOriginIds.end(), nearest_neighbor_id.begin(), nearest_neighbor_id.end()); // appending to end
            }
        }

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != rOriginIds.size()) {
            rLocalMappingMatrix.resize(1, rOriginIds.size(), false);
        }
        const double mapping_weight = 1.0/rOriginIds.size();
        for (IndexType i=0; i<rOriginIds.size(); ++i) {
            rLocalMappingMatrix(0,i) = mapping_weight;
        }

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    } else {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void NearestNeighborLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    rOStream << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 3) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}
/***********************************************************************************/
/***********************************************************************************/

void NearestNeighborLocalSystem::SetPairingStatusForPrinting()
{
    if (mPairingStatus == SearchLocalSystem::PairingStatus::Approximation) {
        mpNode->SetValue(PAIRING_STATUS, 0);
    } else {
        mpNode->SetValue(PAIRING_STATUS, -1);
    }
}

///@} addtogroup block
}  // namespace Kratos.