//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    /// Factory functions

    /** Generates a sample Quadrilateral2D4
     * Generates a quadrilateral defined by four random
     * points in the space.
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateQuadrilateral2D4(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointD = GeneratePoint<TPointType>()
        )
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        PointA,
        PointB,
        PointC,
        PointD
      ));
    }

    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral on the origin
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateRightQuadrilateral2D4()
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0)
      ));
    }

    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral centered on the origin
     * and 45deg rotated
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateDiagQuadrilateral2D4()
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0),
        GeneratePoint<TPointType>(-1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 0.0,-1.0, 0.0)
      ));
    }

    /** Generates a sample Quadrilateral2D4
     * Generates a quadrilateral with trapezoidal shape
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateTrapezoidalQuadrilateral2D4()
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 4.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 2.0, 3.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 3.0, 0.0)
      ));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node>();
        KRATOS_CHECK_NEAR(geom->Area(), 1.0, TOLERANCE);
    }

    /** Test a box and quadrilateral HasIntersection which should give true
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4NodeBoxIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node>();
        Point point_1 (-0.3, 0.8, 0.0);
        Point point_2 ( 0.2, 1.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3 ( 0.9, 1.5, 0.0);
        Point point_4 ( 1.1, 0.8, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5 ( 0.9,-0.8, 0.0);
        Point point_6 ( 1.1, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));

        Point point_7 (-0.3, 0.1, 0.0);
        Point point_8 ( 0.2,-0.6, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_7, point_8));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4EdgeBoxIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateDiagQuadrilateral2D4<Node>();
        Point point_1 ( 0.2, 0.2, 0.0 );
        Point point_2 ( 1.0, 1.0, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3 (-0.2, 0.2, 0.0 );
        Point point_4 (-0.9, 0.9, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5 (-0.2,-0.2, 0.0 );
        Point point_6 (-0.9,-0.9, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));

        Point point_7 ( 0.2,-0.2, 0.0 );
        Point point_8 ( 1.0,-1.0, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_7, point_8));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4BoxNoIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateDiagQuadrilateral2D4<Node>();
        Point point_1 ( 0.7, 0.4, 0.0 );
        Point point_2 ( 1.0, 1.2, 0.0 );
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

    /** Tests the PointLocalCoordinates for Quadrilateral2D4.
     * Tests the PointLocalCoordinates for Quadrilateral2D4.
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node>();

        Point TestPointA(1.0, 1.0, 0.0);
        Point TestPointB(0.5, 0.5, 0.0);
        Point TestResultA(0.0, 0.0, 0.0);
        Point TestResultB(0.0, 0.0, 0.0);

        geom->PointLocalCoordinates(TestResultA, TestPointA);
        geom->PointLocalCoordinates(TestResultB, TestPointB);

        // Test transformation in the edge
        KRATOS_CHECK_NEAR(TestResultA[0], 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultA[1], 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultA[2], 0.0, TOLERANCE);

        // Test transformation in the center
        KRATOS_CHECK_NEAR(TestResultB[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultB[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultB[2], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateRightQuadrilateral2D4<Node>();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 4.0;
      coord[2] = 0.0;
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.09375, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.28125, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.46875, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.15625, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node>();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node>();

        Vector lumping_factors(4);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.25, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.25, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.25, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.25, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4LumpingFactorsTrapezoidalShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateTrapezoidalQuadrilateral2D4<Node>();

        Vector lumping_factors(4);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.3, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.3, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.2, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.2, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.325, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.325, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.175, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.175, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.4, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.4, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.1, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.1, TOLERANCE);
    }

} // namespace Testing
} // namespace Kratos
