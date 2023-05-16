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
#include <string>
#include <vector>
#include <iomanip>
#include <numeric>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/expression_iterator.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/xml_utilities/xml_writer.h"

// Include base h
#include "xml_ostream_writer.h"

namespace Kratos {


const std::string XmlOStreamWriter::GetTabbing(const IndexType Level)
{
    std::stringstream ss_tabbing;
    for (IndexType i = 0; i < Level; ++i) {
        ss_tabbing << "   ";
    }
    return ss_tabbing.str();
}

template<class TExpressionType>
void XmlOStreamWriter::WriteDataElementAscii(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level,
    const std::vector<Expression::Pointer>& rExpressions)
{
    using exp_itr_type = ExpressionIterator<TExpressionType>;

    using data_itr_type = typename exp_itr_type::ConstIteratorType;

    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << " format=\"ascii\">\n" << tabbing;

    std::vector<TExpressionType*> transformed_expressions(rExpressions.size());
    std::transform(rExpressions.begin(), rExpressions.end(),
                   transformed_expressions.begin(), [](auto& pExpression) {
                       return dynamic_cast<TExpressionType*>(&*(pExpression));
                   });

    for (const auto& p_expression : transformed_expressions) {
        exp_itr_type expression_iterator(p_expression);
        auto data_begin = expression_iterator.cbegin();
        auto data_end   = expression_iterator.cend();
        for (data_itr_type itr = data_begin; itr != data_end; ++itr) {
            if constexpr(std::is_same_v<typename exp_itr_type::DataType, char>) {
                mrOStream << "  " << static_cast<int>(*(itr));
            } else {
                mrOStream << "  " << *itr;
            }
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

template<class TExpressionType>
void XmlOStreamWriter::WriteDataElementBinary(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level,
    const std::vector<Expression::Pointer>& rExpressions)
{
    using exp_itr_type = ExpressionIterator<TExpressionType>;

    using data_type = typename exp_itr_type::DataType;

    using data_itr_type = typename exp_itr_type::ConstIteratorType;

    WriteAttributes(rTagName, rAttributes, Level);

    std::vector<TExpressionType*> transformed_expressions(rExpressions.size());
    std::transform(rExpressions.begin(), rExpressions.end(),
                   transformed_expressions.begin(), [](auto& pExpression) {
                       return dynamic_cast<TExpressionType*>(&*(pExpression));
                   });

    if (rExpressions.size() == 0) {
        mrOStream << " format=\"binary\"/>\n";
        return;
    } else {
        data_type min_value{std::numeric_limits<data_type>::max()}, max_value{std::numeric_limits<data_type>::lowest()};
        for (const auto& p_expression : transformed_expressions) {
            data_itr_type data_itr = exp_itr_type(p_expression).cbegin();
            const auto values = IndexPartition<IndexType>(p_expression->GetFlattenedShapeSize() * p_expression->NumberOfEntities()).for_each<CombinedReduction<MinReduction<data_type>, MaxReduction<data_type>>>([&data_itr](const IndexType Index) {
                const data_type value = *(data_itr + Index);
                return std::make_tuple(value, value);
            });
            min_value = std::min(std::get<0>(values), min_value);
            max_value = std::max(std::get<1>(values), max_value);
        }

        const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
        if constexpr(std::is_same_v<data_type, char>) {
            mrOStream << " format=\"binary\" RangeMin=\"" << static_cast<int>(min_value) << "\" RangeMax=\"" << static_cast<int>(max_value) << "\">\n" << tabbing << "  ";
        } else {
            mrOStream << " format=\"binary\" RangeMin=\"" << min_value << "\" RangeMax=\"" << max_value << "\">\n" << tabbing << "  ";
        }

    }

    constexpr char base64_map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    const IndexType total_entities = std::accumulate(rExpressions.begin(), rExpressions.end(), 0U, [](const IndexType LHS, const auto& pExpression) { return LHS + pExpression->NumberOfEntities();});

    using writing_data_type = data_type;
    constexpr IndexType data_type_size = sizeof(writing_data_type);
    constexpr IndexType size_type_size = sizeof(IndexType);
    const IndexType total_number_of_values = total_entities * rExpressions[0]->GetFlattenedShapeSize();
    const IndexType size_type_triplets = size_type_size / 3;
    const IndexType total_data_size = total_number_of_values * data_type_size;
    const IndexType total_number_of_triplets = (total_data_size + size_type_size)  / 3;

    IndexType byte_index = 0;
    auto p_expression = transformed_expressions.data();
    data_itr_type data_itr = exp_itr_type(*p_expression).cbegin();
    data_itr_type data_end = exp_itr_type(*p_expression).cend();
    writing_data_type current_value =  *data_itr;

    auto get_next_byte = [&]() -> char {
        if (byte_index == data_type_size) {
            byte_index = 0;
            ++data_itr;

            if (data_itr == data_end) {
                ++p_expression;
                data_itr = exp_itr_type(*p_expression).cbegin();
                data_end = exp_itr_type(*p_expression).cend();
            }

            current_value = *data_itr;
        }

        const char byte = *(reinterpret_cast<const char*>(&current_value) + byte_index++);

        return byte;
    };

    auto write_encoded_triplet = [&](const std::array<char, 3>& bytes, size_t padding) {
        char tmp[5] = {
            base64_map[(bytes[0] & 0xfc) >> 2],
            base64_map[((bytes[0] & 0x03) << 4) + ((bytes[1] & 0xf0) >> 4)],
            base64_map[((bytes[1] & 0x0f) << 2) + ((bytes[2] & 0xc0) >> 6)],
            base64_map[bytes[2] & 0x3f], '\0'};

        std::fill(tmp + 4 - padding, tmp + 4, '=');

        mrOStream << tmp;
    };

    // first write the total number of bytes in the array this will be of 8 bytes
    auto total_data_size_begin = reinterpret_cast<const char*>(&total_data_size);
    IndexType local_index = 0;
    IndexType initial_number_of_triplets = 0;
    for (IndexType i = 0; i < size_type_triplets; ++i) {
        write_encoded_triplet({*(total_data_size_begin + local_index++), *(total_data_size_begin + local_index++), *(total_data_size_begin + local_index++)}, 0);
        ++initial_number_of_triplets;
    }

    std::array<char, 3> data;
    switch (size_type_size % 3) {
        case 1:
            data[0] = *(total_data_size_begin + local_index++);
            data[1] = get_next_byte();
            data[2] = get_next_byte();
            ++initial_number_of_triplets;
            write_encoded_triplet(data, 0);
            break;
        case 2:
            data[0] = *(total_data_size_begin + local_index++);
            data[1] = *(total_data_size_begin + local_index++);
            data[2] = get_next_byte();
            ++initial_number_of_triplets;
            write_encoded_triplet(data, 0);
            break;
        default:
            break;
    }

    for (IndexType i = initial_number_of_triplets; i < total_number_of_triplets; ++i) {
        write_encoded_triplet({get_next_byte(), get_next_byte(), get_next_byte()}, 0);
    }

    const IndexType number_of_bytes_remaining = (total_data_size + size_type_size) % 3;
    if (number_of_bytes_remaining != 0) {
        std::array<char, 3> bytes{'\0', '\0', '\0'};

        for (IndexType i = 0; i < number_of_bytes_remaining; ++i) {
            bytes[i] = get_next_byte();
        }
        write_encoded_triplet(bytes, 3 - number_of_bytes_remaining);
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

XmlOStreamWriter::XmlOStreamWriter(
    std::ostream& rOStream,
    const WriterFormat OutputFormat,
    const IndexType Precision)
    : mrOStream(rOStream),
      mOutputFormat(OutputFormat)
{
    mrOStream << std::scientific << std::setprecision(Precision);
}

void XmlOStreamWriter::WriteAttributes(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "<" << rTagName;
    if (rAttributes.size() > 0) {
        for (const auto& r_pair : rAttributes) {
            mrOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
        }
    }
}

void XmlOStreamWriter::WriteElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const IndexType Level,
    const bool IsEmptyElement)
{
    WriteAttributes(rTagName, rAttributes, Level);

    if (IsEmptyElement) {
        mrOStream << "/>\n";
    } else {
        mrOStream << ">\n";
    }

}

void XmlOStreamWriter::CloseElement(
    const std::string& rTagName,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "</" << rTagName << ">\n";
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<Expression::Pointer>& rExpressions,
    const IndexType Level)
{

    switch (mOutputFormat){
        case ASCII:
            if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<char>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<double>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
            } else {
                WriteDataElementAscii<Expression>(rTagName, rAttributes, Level, rExpressions);
            }
            break;
        case BINARY:
            if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<char>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<double>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
            } else {
                WriteDataElementBinary<Expression>(rTagName, rAttributes, Level, rExpressions);
            }
            break;
    }
}

} // namespace Kratos