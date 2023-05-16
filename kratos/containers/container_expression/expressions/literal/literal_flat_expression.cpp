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

template<class RawTDataType>
LiteralFlatExpression<RawTDataType>::LiteralFlatExpression(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : Expression(NumberOfEntities),
      mShape(rShape),
      mData(NumberOfEntities * this->GetFlattenedShapeSize())
{
}

template<class RawTDataType>
LiteralFlatExpression<RawTDataType>::LiteralFlatExpression(
    RawTDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : Expression(NumberOfEntities),
      mShape(rShape),
      mData(pDataBegin, NumberOfEntities * this->GetFlattenedShapeSize())
{
}

template<class RawTDataType>
typename LiteralFlatExpression<RawTDataType>::Pointer LiteralFlatExpression<RawTDataType>::Create(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<RawTDataType>>(NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<RawTDataType>>(NumberOfEntities, rShape);
    }
}

template<class RawTDataType>
typename LiteralFlatExpression<RawTDataType>::Pointer LiteralFlatExpression<RawTDataType>::Create(
    RawTDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<RawTDataType>>(pDataBegin, NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<RawTDataType>>(pDataBegin, NumberOfEntities, rShape);
    }
}

template<class RawTDataType>
void LiteralFlatExpression<RawTDataType>::SetData(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex,
    const RawTDataType Value)
{
    *(mData.data_begin() + EntityDataBeginIndex + ComponentIndex) = Value;
}

template<class RawTDataType>
const std::vector<std::size_t> LiteralFlatExpression<RawTDataType>::GetShape() const
{
    return mShape;
}

template<>
std::string LiteralFlatExpression<char>::Info() const
{
    std::stringstream msg;
    msg << "CharVec" << mShape;
    return msg.str();
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

template<class RawTDataType>
double LiteralScalarFlatExpression<RawTDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return *(this->mData.data_begin() + EntityIndex);
}

template<class RawTDataType>
double LiteralNonScalarFlatExpression<RawTDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return *(this->mData.data_begin() + EntityDataBeginIndex + ComponentIndex);
}

// template instantiations
template class LiteralFlatExpression<char>;
template class LiteralScalarFlatExpression<char>;
template class LiteralNonScalarFlatExpression<char>;

template class LiteralFlatExpression<int>;
template class LiteralScalarFlatExpression<int>;
template class LiteralNonScalarFlatExpression<int>;

template class LiteralFlatExpression<double>;
template class LiteralScalarFlatExpression<double>;
template class LiteralNonScalarFlatExpression<double>;
} // namespace Kratos