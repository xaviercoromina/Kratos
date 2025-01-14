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
#include "spatial_containers/spatial_search_result.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultDefaultConstruction, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    KRATOS_CHECK_EQUAL(result.IsObjectFound(), false);
    KRATOS_CHECK_EQUAL(result.IsDistanceCalculated(), false);
    KRATOS_CHECK_EQUAL(result.Get(), nullptr);
    KRATOS_CHECK_EQUAL(result.GetDistance(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultObjectFound, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetIsObjectFound(true);
    KRATOS_CHECK_EQUAL(result.IsObjectFound(), true);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultDistanceCalculated, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetIsDistanceCalculated(true);
    KRATOS_CHECK_EQUAL(result.IsDistanceCalculated(), true);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultPointer, KratosCoreFastSuite)
{
    Element element = Element();
    auto result = SpatialSearchResult<GeometricalObject>(&element);
    GeometricalObject* ptr = &element;
    KRATOS_CHECK_EQUAL(result.Get().get(), ptr);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultDistance, KratosCoreFastSuite)
{
    auto result = SpatialSearchResult<GeometricalObject>();
    result.SetDistance(3.14);
    KRATOS_CHECK_EQUAL(result.GetDistance(), 3.14);
}

}  // namespace Kratos::Testing
