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
#include <string>
#include <variant>
#include <vector>
#include <utility>

// External includes

// Project includes
#include "includes/io.h"
#include "includes/model_part.h"
#include "containers/variable.h"
#include "containers/flags.h"
#include "containers/container_expression/container_expression.h"


namespace Kratos {
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using SupportedVariables = std::variant<
                                    const Variable<int>*,
                                    const Variable<double>*,
                                    const Variable<array_1d<double, 3>>*,
                                    const Variable<array_1d<double, 4>>*,
                                    const Variable<array_1d<double, 6>>*,
                                    const Variable<array_1d<double, 9>>*>;

    using SupportedContainerExpressions = std::variant<
                                            ContainerExpression<ModelPart::NodesContainerType>::Pointer,
                                            ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
                                            ContainerExpression<ModelPart::ElementsContainerType>::Pointer>;

    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    KRATOS_DEFINE_LOCAL_FLAG( NODES );
    KRATOS_DEFINE_LOCAL_FLAG( CONDITIONS );
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS );

    ///@}
    ///@name Life cycle
    ///@{

    VtuOutput(
        ModelPart& rModelPart,
        const bool IsInitialConfiguration,
        const IndexType Precision);

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    void AddHistoricalVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void AddNonHistoricalVariable(
        const Variable<TDataType>& rVariable,
        const Flags& rEntityFlags);

    void AddFlagVariable(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        const Flags& rEntityFlags);

    template <class TContainerType>
    void AddContainerExpression(
        const std::string& rExpressionName,
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression);

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    const bool mIsInitialConfiguration;

    const IndexType mPrecision;

    bool mIsConditionsConsidered;

    bool mIsElementsConsidered;

    std::unordered_map<IndexType, IndexType> mKratosVtuIndicesMap;

    std::vector<SupportedVariables> mHistoricalVariablesList;

    std::vector<SupportedVariables> mNonHistoricalNodalVariablesList;

    std::vector<SupportedVariables> mNonHistoricalCellVariablesList;

    std::vector<std::pair<std::string, const Flags*>> mNodalFlagsList;

    std::vector<std::pair<std::string, const Flags*>> mCellFlagsList;

    std::vector<SupportedContainerExpressions> mContainerExpressionsList;

    ///@}
};
} // namespace Kratos