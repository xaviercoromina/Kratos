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
#include "utilities/projection_utilities.h"
#include "searching/search_interface_info/nearest_element_interface_info.h"

namespace Kratos
{
///@name Kratos Classes
///@{

void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject, false);
    mNumSearchResults++;
}

/***********************************************************************************/
/***********************************************************************************/

void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NearestElementInterfaceInfo::SaveSearchResult(const InterfaceObject& rInterfaceObject,
                                                   const bool ComputeApproximation)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    double proj_dist;

    const Point point_to_proj(this->Coordinates());

    Vector shape_function_values;
    std::vector<int> eq_ids;

    ProjectionUtilities::PairingIndex pairing_index;

    const bool is_full_projection = ProjectionUtilities::ComputeProjection(*p_geom, point_to_proj, mLocalCoordTol, shape_function_values, eq_ids, proj_dist, pairing_index, ComputeApproximation);

    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        if (!ComputeApproximation) {
            return;
        } else {
            SetIsApproximation();
        }
    }

    const std::size_t num_values = shape_function_values.size();
    KRATOS_ERROR_IF_NOT(num_values == eq_ids.size()) << "Number of equation-ids is not the same as the number of ShapeFunction values, something went wrong!" << std::endl;

    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        if (mShapeFunctionValues.size() != num_values) mShapeFunctionValues.resize(num_values);
        for (std::size_t i=0; i<num_values; ++i) {
            mShapeFunctionValues[i] = shape_function_values[i];
        }
    }
}

///@} addtogroup block
}  // namespace Kratos.