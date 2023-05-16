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

#pragma once

// System includes
#include <vector>
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Expression to hold a literal with a flattened data structure
 *
 * This expression holds data of a literal with rShape by
 * flattening the data structure to a vector. This allocates
 * a Vector with large sizes.
 *
 * @tparam RawTDataType        Raw data type of the literal. Can be char, int or double
 */
template<class RawTDataType = double>
class KRATOS_API(KRATOS_CORE) LiteralFlatExpression : public Expression {
public:
    ///@name Type definition
    ///@{

    using DataType = RawTDataType;

    using Pointer = Kratos::intrusive_ptr<LiteralFlatExpression<RawTDataType>>;

    ///@}
    ///@name Life cycle
    ///@{

    LiteralFlatExpression(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    LiteralFlatExpression(
        RawTDataType* pDataBegin,
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief A specific create method is used in here to return a pointer to the LiteralFlatExpression
     *
     * This method is used to return the derrived class LiteralFlatExpression pointer
     * so that additional SetData method can be accessed to initialize
     * the vector data.
     *
     * @param NumberOfEntities                  Number of entities (eg. nodes, elements, etc.) in the container.
     * @param rShape                            Shape of the data in each entity.
     * @return LiteralFlatExpression::Pointer   Returns an intrusive pointer to LiteralFlatExpression.
     */
    static LiteralFlatExpression<RawTDataType>::Pointer Create(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    static LiteralFlatExpression<RawTDataType>::Pointer Create(
        RawTDataType* pDataBegin,
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    void SetData(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex,
        const RawTDataType Value);

    const std::vector<IndexType> GetShape() const override;

    inline IndexType Size() const noexcept { return mData.Size(); }

    inline RawTDataType* data_begin() noexcept { return mData.data_begin(); }

    inline RawTDataType const* data_begin() const noexcept { return mData.data_begin(); }

    inline RawTDataType const* data_end() const noexcept { return mData.data_end(); }

    std::string Info() const override;

    ///@}
protected:
    ///@name Protected classes
    ///@{

    /**
     * @brief This class is used represent the flattened data array
     *
     * Instances of this class used to represent the data in the flat array.
     * This representation can be done in two ways:
     *
     *      1) If the data is moved, then this instance will hold a pointer
     *         to the double flattened data array. In this case it is important to
     *         make sure that the double flattened data array lives until this expression
     *         is destructed. In this case mIsManaged will be false.
     *      2) If the data is copied, then a flattened array for the required size is
     *         created at the construction in the heap, and will be destructed when the
     *         @ref Data instance is destructed. Hence, the copy of the data is stored
     *         and life time is managed by this instance.
     */
    class Data
    {
    public:
        ///@}
        ///@name Life cycle
        ///@{

        /**
         * @brief Construct a new Data object where the array memory is managed by the object.
         *
         * @param Size      Size of the allocated array.
         */
        Data(const IndexType Size): mpBegin(new RawTDataType[Size]), mIsManaged(true), mSize(Size) {}

        /**
         * @brief Construct a new Data object, where the underlying array memory is not managed by the object.
         *
         * @param pBegin    Pointer to the memory array.
         */
        Data(RawTDataType* pBegin, const IndexType Size): mpBegin(pBegin), mIsManaged(false), mSize(Size) {}

        ~Data() { if (mIsManaged) { delete[] mpBegin; } }

        ///@}
        ///@name Operators
        ///@{

        inline RawTDataType* data_begin() noexcept { return mpBegin; }

        inline RawTDataType const* data_begin() const noexcept { return mpBegin; }

        inline RawTDataType const* data_end() const noexcept { return mpBegin + mSize; }

        inline IndexType Size() const noexcept { return mSize; }

        ///@}
    private:
        ///@name Private member variables
        ///@{

        RawTDataType* mpBegin;

        const bool mIsManaged;

        const IndexType mSize;

        ///@}
    };

    ///@}
    ///@name Protected member variables
    ///@{

    const std::vector<IndexType> mShape;

    Data mData;

    ///@}
};

template<class RawTDataType = double>
class LiteralScalarFlatExpression : public LiteralFlatExpression<RawTDataType>
{
public:

    using IndexType = std::size_t;

    using LiteralFlatExpression<RawTDataType>::LiteralFlatExpression;

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

};

template<class RawTDataType = double>
class LiteralNonScalarFlatExpression : public LiteralFlatExpression<RawTDataType>
{
public:

    using IndexType = std::size_t;

    using LiteralFlatExpression<RawTDataType>::LiteralFlatExpression;

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

};

} // namespace Kratos