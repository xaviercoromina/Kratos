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
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
///@name Kratos namespaces
///@{

/**
 * @brief This namespace defines the utilities for the projection
 */
namespace ProjectionUtilities
{
///@name  Enum's
///@{

/**
 * @brief This enum is used to identify the type of pairing
 */
enum class PairingIndex
{
    Volume_Inside   = -1,
    Volume_Outside  = -2,
    Surface_Inside  = -3,
    Surface_Outside = -4,
    Line_Inside     = -5,
    Line_Outside    = -6,
    Closest_Point   = -7,
    Unspecified     = -8
};

///@}
///@name Type Definitions
///@{

/// The index type definition
using SizeType = std::size_t;

/// The index type definition
using IndexType = std::size_t;

/// The geometry type definition
using GeometryType = Geometry<Node<3>>;

///@}
///@name Operations
///@{

/**
 * @brief This method projects a point onto a line
 * @param rGeometry The geometry to be projected onto
 * @param rPointToProject The point to be projected
 * @param LocalCoordTol The local coordinate tolerance
 * @param rShapeFunctionValues The shape function values
 * @param rEquationIds The equation ids
 * @param rProjectionDistance The projection distance
 * @param ComputeApproximation If the approximation should be computed
 * @return The pairing index
 */
PairingIndex KRATOS_API(KRATOS_CORE) ProjectOnLine(const GeometryType& rGeometry,
                           const Point& rPointToProject,
                           const double LocalCoordTol,
                           Vector& rShapeFunctionValues,
                           std::vector<int>& rEquationIds,
                           double& rProjectionDistance,
                           const bool ComputeApproximation=true);

/**
 * @brief This method projects a point onto a surface
 * @param rGeometry The geometry to be projected onto
 * @param rPointToProject The point to be projected
 * @param LocalCoordTol The local coordinate tolerance
 * @param rShapeFunctionValues The shape function values
 * @param rEquationIds The equation ids
 * @param rProjectionDistance The projection distance
 * @param ComputeApproximation If the approximation should be computed
 * @return The pairing index
 */
PairingIndex KRATOS_API(KRATOS_CORE) ProjectOnSurface(const GeometryType& rGeometry,
                     const Point& rPointToProject,
                     const double LocalCoordTol,
                     Vector& rShapeFunctionValues,
                     std::vector<int>& rEquationIds,
                     double& rProjectionDistance,
                     const bool ComputeApproximation=true);

/**
 * @brief This method projects a point into a volume
 * @param rGeometry The geometry to be projected onto
 * @param rPointToProject The point to be projected
 * @param LocalCoordTol The local coordinate tolerance
 * @param rShapeFunctionValues The shape function values
 * @param rEquationIds The equation ids
 * @param rProjectionDistance The projection distance
 * @param ComputeApproximation If the approximation should be computed
 * @return The pairing index
 */
PairingIndex KRATOS_API(KRATOS_CORE) ProjectIntoVolume(const GeometryType& rGeometry,
                               const Point& rPointToProject,
                               const double LocalCoordTol,
                               Vector& rShapeFunctionValues,
                               std::vector<int>& rEquationIds,
                               double& rProjectionDistance,
                               const bool ComputeApproximation=true);

/**
 * @brief This method projects a point onto a geometry
 * @param rGeometry The geometry to be projected onto
 * @param rPointToProject The point to be projected
 * @param LocalCoordTol The local coordinate tolerance
 * @param rShapeFunctionValues The shape function values
 * @param rEquationIds The equation ids
 * @param rProjectionDistance The projection distance
 * @param ComputeApproximation If the approximation should be computed
 * @return The pairing index
 */
bool KRATOS_API(KRATOS_CORE) ComputeProjection(const GeometryType& rGeometry,
                       const Point& rPointToProject,
                       const double LocalCoordTol,
                       Vector& rShapeFunctionValues,
                       std::vector<int>& rEquationIds,
                       double& rProjectionDistance,
                       PairingIndex& rPairingIndex,
                       const bool ComputeApproximation=true);
///@}

}  // namespace ProjectionUtilities.

}  // namespace Kratos.
