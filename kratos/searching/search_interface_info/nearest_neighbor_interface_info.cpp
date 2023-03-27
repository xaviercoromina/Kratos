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
#include "searching/search_interface_info/nearest_neighbor_interface_info.h"

namespace Kratos
{
///@name Kratos Classes
///@{

void NearestNeighborInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SetLocalSearchWasSuccessful();

    const double neighbor_distance = this->ComputeDistance(rInterfaceObject);

    if (neighbor_distance < mNearestNeighborDistance) {
        mNearestNeighborDistance = neighbor_distance;
        mNearestNeighborId.resize(1);
        mNearestNeighborId[0] = rInterfaceObject.pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    } else if (neighbor_distance == mNearestNeighborDistance) {
        mNearestNeighborId.push_back(rInterfaceObject.pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID));
    }
}

///@} addtogroup block
}  // namespace Kratos.