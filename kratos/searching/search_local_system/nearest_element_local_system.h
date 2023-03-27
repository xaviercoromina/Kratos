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
#include "utilities/projection_utilities.h"
#include "searching/search_local_system/search_local_system.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) NearestElementLocalSystem
    : public SearchLocalSystem
{
public:

    using NodePointerType = SearchLocalSystem::NodePointerType;

    explicit NearestElementLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      SearchLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        return mpNode->Coordinates();
    }

    SearchLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestElementLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

    bool IsDoneSearching() const override;

private:
    NodePointerType mpNode;
    mutable ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;

};

///@} addtogroup block
}  // namespace Kratos.