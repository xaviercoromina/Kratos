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
#include <type_traits>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/xml_utilities/xml_element.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class XmlEntitiesType, class TDataWriterType>
class XmlEntityData : public XmlElement
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = XmlElement;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlEntityData);

    ///@}
    ///@name Life cycle
    ///@{

    XmlEntityData(
        const typename XmlEntitiesType::Pointer mXmlEntities,
        const std::string& rDataName,
        const IndexType NumberOfComponents,
        const TDataWriterType& rDataWriter,
        const IndexType Precision = 9)
        : BaseType("DataArray"),
          mpXmlEntities(mXmlEntities),
          mDataWriter(rDataWriter),
          mPrecision(Precision)
    {
        static_assert(std::is_same_v<typename TDataWriterType::DataType, bool> ||
                      std::is_same_v<typename TDataWriterType::DataType, int> ||
                      std::is_same_v<typename TDataWriterType::DataType, double>,
                      "TDataRetrievalFunctor::DataType should be bool, int or double.");

        if constexpr(std::is_same_v<typename TDataWriterType::DataType, bool>) {
            AddAttribute("type", "UInt8");
        } else if constexpr(std::is_same_v<typename TDataWriterType::DataType, int>) {
            AddAttribute("type", "Int32");
        } else if constexpr(std::is_same_v<typename TDataWriterType::DataType, double>) {
            AddAttribute("type", "Float32");
        }

        AddAttribute("Name", rDataName);

        if (NumberOfComponents > 1) {
            AddAttribute("NumberOfComponents", std::to_string(NumberOfComponents));
        }

        AddAttribute("Format", "ascii");
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const typename XmlEntitiesType::Pointer mpXmlEntities;

    const TDataWriterType mDataWriter;

    const IndexType mPrecision;

    ///@}
    ///@name Private opeations
    ///@{

    bool HasElementData() const override
    {
        return true;
    }

    void WrtieElementData(
        std::ostream& rOutputStream,
        const std::string& rTabbing) const override
    {
        if constexpr(std::is_same_v<typename TDataRetrievalFunctor::DataType, double>) {
            rOutputStream << std::setprecision(mPrecision) << std::scientific;
        } else {
            rOutputStream << std::fixed;
        }

        mDataWriter.WriteData(rOutputStream, mXmlEntities->GetEntities(), rTabbing);
    }

    ///@}
};

} // namespace Kratos