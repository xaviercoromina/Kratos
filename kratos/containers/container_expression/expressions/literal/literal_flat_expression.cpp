//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <sstream>

// Project includes

// Include base h
#include "literal_flat_expression.h"

namespace Kratos {

template<class TDataType>
LiteralFlatExpression<TDataType>::LiteralFlatExpression(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : mShape(rShape),
      mData(NumberOfEntities * this->GetFlattenedSize())
{
}

template<class TDataType>
LiteralFlatExpression<TDataType>::LiteralFlatExpression(
    TDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : mShape(rShape),
      mData(pDataBegin)
{
}

template<class TDataType>
typename LiteralFlatExpression<TDataType>::Pointer LiteralFlatExpression<TDataType>::Create(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<TDataType>>(NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<TDataType>>(NumberOfEntities, rShape);
    }
}

template<class TDataType>
typename LiteralFlatExpression<TDataType>::Pointer LiteralFlatExpression<TDataType>::Create(
    TDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<TDataType>>(pDataBegin, NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<TDataType>>(pDataBegin, NumberOfEntities, rShape);
    }
}

template<class TDataType>
void LiteralFlatExpression<TDataType>::SetData(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex,
    const TDataType Value)
{
    mData[EntityDataBeginIndex + ComponentIndex] = Value;
}

template<class TDataType>
const std::vector<std::size_t> LiteralFlatExpression<TDataType>::GetShape() const
{
    return mShape;
}

template<>
std::string LiteralFlatExpression<int>::Info() const
{
    std::stringstream msg;
    msg << "IntVec" << mShape;
    return msg.str();
}

template<>
std::string LiteralFlatExpression<double>::Info() const
{
    std::stringstream msg;
    msg << "DoubleVec" << mShape;
    return msg.str();
}

template<class TDataType>
double LiteralScalarFlatExpression<TDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mData[EntityIndex];
}

template<class TDataType>
double LiteralNonScalarFlatExpression<TDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mData[EntityDataBeginIndex + ComponentIndex];
}

// template instantiations
template class LiteralFlatExpression<int>;
template class LiteralScalarFlatExpression<int>;
template class LiteralNonScalarFlatExpression<int>;

template class LiteralFlatExpression<double>;
template class LiteralScalarFlatExpression<double>;
template class LiteralNonScalarFlatExpression<double>;
} // namespace Kratos