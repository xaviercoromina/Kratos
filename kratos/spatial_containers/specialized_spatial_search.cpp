//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |  (   | |  (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "spatial_containers/specialized_spatial_search.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

template<>
void PointObject<Node<3>>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void PointObject<Condition>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void PointObject<Element>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
}

template class PointObject<Node<3>>;
template class PointObject<Condition>;
template class PointObject<Element>;

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Element>;
    using PointTypePointer = typename PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureElements, rInputElements, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputElements, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputElements, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusExclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Node<3>>;
    using PointTypePointer = typename PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureNodes, rInputNodes, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputNodes, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputNodes, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusExclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Condition>;
    using PointTypePointer = typename PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureConditions, rInputConditions, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputConditions, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        ParallelSearch(rInputConditions, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
Parameters SpecializedSpatialSearch<TSearchBackend>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "allocation_size" : 1000,
        "bucket_size"     : 4
    })" );

    return default_parameters;
}

template class SpecializedSpatialSearch<SpatialContainer::KDTree>;
template class SpecializedSpatialSearch<SpatialContainer::Octree>;
template class SpecializedSpatialSearch<SpatialContainer::BinsStatic>;
template class SpecializedSpatialSearch<SpatialContainer::BinsDynamic>;
// template class SpecializedSpatialSearch<SpatialContainer::BinsStaticObjects>;
// template class SpecializedSpatialSearch<SpatialContainer::BinsDynamicObjects>;
} // namespace Kratos.