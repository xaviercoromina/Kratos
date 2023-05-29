//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/element.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos::Testing 
{

/** @brief Test case for default construction of SpatialSearchResult.
 *  @tparam GeometricalObject The type of geometrical object.
 */
KRATOS_TEST_CASE_IN_SUITE(DefaultConstructionSpatialSearchResult, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();

    // Check initial conditions
    KRATOS_CHECK_EQUAL(result.IsObjectFound(), false);
    KRATOS_CHECK_EQUAL(result.IsDistanceCalculated(), false);
    KRATOS_CHECK_EQUAL(result.Get(), nullptr);
    KRATOS_CHECK_EQUAL(result.GetDistance(), 0.0);
}

/** @brief Test case for setting the object found flag of SpatialSearchResult.
 *  @tparam GeometricalObject The type of geometrical object.
 */
KRATOS_TEST_CASE_IN_SUITE(ObjectFoundSpatialSearchResult, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetIsObjectFound(true);

    // Check if the object found flag is set
    KRATOS_CHECK_EQUAL(result.IsObjectFound(), true);
}

/** @brief Test case for setting the distance calculated flag of SpatialSearchResult.
 *  @tparam GeometricalObject The type of geometrical object.
 */
KRATOS_TEST_CASE_IN_SUITE(DistanceCalculatedSpatialSearchResult, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetIsDistanceCalculated(true);

    // Check if the distance calculated flag is set
    KRATOS_CHECK_EQUAL(result.IsDistanceCalculated(), true);
}

/** @brief Test case for assigning a pointer to SpatialSearchResult.
 *  @tparam GeometricalObject The type of geometrical object.
 */
KRATOS_TEST_CASE_IN_SUITE(PointerSpatialSearchResult, KratosCoreFastSuite)
{
    Element element = Element();
    auto result = SpatialSearchResult<GeometricalObject>(&element);
    GeometricalObject* ptr = &element;

    // Check if the assigned pointer is correct
    KRATOS_CHECK_EQUAL(result.Get().get(), ptr);
}

/** @brief Test case for setting the distance of SpatialSearchResult.
 *  @tparam GeometricalObject The type of geometrical object.
 */
KRATOS_TEST_CASE_IN_SUITE(DistanceSpatialSearchResult, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetDistance(3.14);

    // Check if the distance is set correctly
    KRATOS_CHECK_EQUAL(result.GetDistance(), 3.14);
}

}  // namespace Kratos::Testing