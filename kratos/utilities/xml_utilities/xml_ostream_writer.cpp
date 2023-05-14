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
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
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
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level,
    const std::vector<TExpressionType>& rExpressions)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    mrOStream << " Format=\"ascii\">\n";

    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);

    for (IndexType i = 0; i < rExpressions.size(); ++i) {
        auto& r_expression = *rExpressions[i];
        auto number_of_entities = rNumberOfEntities[i];
        const IndexType number_of_components = r_expression.GetFlattenedSize();

        mrOStream << tabbing;

        if constexpr(std::is_same_v<TExpressionType, Expression::Pointer>) {
            for (IndexType entity_index = 0; entity_index < number_of_entities; ++entity_index) {
                const IndexType entity_start_index = entity_index * number_of_components;
                for (IndexType component_index = 0; component_index < number_of_components; ++component_index) {
                    mrOStream << "  " << r_expression.Evaluate(entity_index, entity_start_index, component_index);
                }
            }
        } else {
            const IndexType data_size = number_of_entities * number_of_components;
            for (IndexType i = 0; i < data_size; ++i) {
                mrOStream << "  " << r_expression[i];
            }
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

// constexpr char base64Map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// template<typename Iterator>
// inline std::string base64Encode( Iterator begin, Iterator end )
// {
//     constexpr size_t size = sizeof( decltype( *begin ) );

//     size_t length = static_cast<size_t>( std::distance( begin, end ) );
//     size_t rawBytes = length * size;
//     size_t encodedBytes = ( rawBytes / 3 + 1 ) * 4;

//     std::string result;

//     result.reserve( encodedBytes );

//     auto it = begin;
//     size_t byteIndex = 0;

//     auto next = [&]( )
//     {
//         char byte = *( reinterpret_cast<const char*>( &( *it ) ) + byteIndex++ );

//         if( byteIndex == size )
//         {
//             it++;
//             byteIndex = 0;
//         }

//         return byte;
//     };

//     auto encodeTriplet = [&]( std::array<char, 3> bytes, size_t padding )
//     {
//         char tmp[5] = { base64Map[(   bytes[0] & 0xfc ) >> 2],
//                         base64Map[( ( bytes[0] & 0x03 ) << 4 ) + ( ( bytes[1] & 0xf0 ) >> 4 )],
//                         base64Map[( ( bytes[1] & 0x0f ) << 2 ) + ( ( bytes[2] & 0xc0 ) >> 6 )],
//                         base64Map[bytes[2] & 0x3f],
//                         '\0' };

//         std::fill( tmp + 4 - padding, tmp + 4, '=' );

//         result += tmp;
//     };

//     for( size_t i = 0; i < rawBytes / 3; ++i )
//     {
//         encodeTriplet( { next( ), next( ), next( ) }, 0 );
//     }

//     if( it != end )
//     {
//         std::array<char, 3> bytes { '\0', '\0', '\0' };

//         size_t remainder = static_cast<size_t>( std::distance( it, end ) ) * size - static_cast<size_t>( byteIndex );

//         for( size_t i = 0; i < remainder; ++i )
//         {
//             bytes[i] = next( );
//         }

//         encodeTriplet( bytes, 3 - remainder );
//     }

//     return result;
// }

template<class TExpressionType>
void XmlOStreamWriter::WriteDataElementBinary(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level,
    const std::vector<TExpressionType>& rExpressions)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    mrOStream << " Format=\"binary\">\n";
    if (rExpressions.size() == 0) {
        CloseElement(rTagName, Level);
        return;
    }

    constexpr char base64_map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    const IndexType total_entities =
        std::accumulate(rNumberOfEntities.begin(), rNumberOfEntities.end(), 0);

    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "  ";

    const IndexType number_of_components = rExpressions[0]->GetFlattenedSize();

    // following data type can be either int or double.
    using data_type = typename TExpressionType::element_type::DataType;
    constexpr IndexType data_type_size = sizeof(data_type);
    const IndexType number_of_values = total_entities * number_of_components;
    const IndexType total_data_size = number_of_values * data_type_size;
    const IndexType number_of_triplets = total_data_size / 3;

    auto itr_expression = rExpressions.begin();
    auto itr_number_of_entities = rNumberOfEntities.begin();
    IndexType byte_index = 0;
    IndexType current_value_index = 0;
    IndexType current_entity_index = 0;
    IndexType current_entity_begin_index = 0;
    IndexType current_component_index = 0;
    data_type current_value;

    if constexpr(std::is_same_v<TExpressionType, Expression::Pointer>) {
        current_value = (*itr_expression)->Evaluate(0, 0, 0);
        KRATOS_WATCH("Writing inefficient")
    } else {
        current_value = (*(*itr_expression))[0];
        KRATOS_WATCH("Writing efficient")
    }

    auto get_next_byte = [&]() -> char {
        if (byte_index == data_type_size) {
            byte_index = 0;
            ++current_component_index;
            ++current_value_index;

            if (current_component_index == number_of_components) {
                current_component_index = 0;
                ++current_entity_index;
                current_entity_begin_index += number_of_components;

                if (current_entity_index == *itr_number_of_entities) {
                    current_value_index = 0;
                    current_entity_index = 0;
                    current_entity_begin_index = 0;
                    ++itr_number_of_entities;
                    ++itr_expression;
                }
            }

            if constexpr(std::is_same_v<TExpressionType, Expression::Pointer>) {
                current_value = (*itr_expression)->Evaluate(
                    current_entity_index, current_entity_begin_index, current_component_index);
            } else {
                current_value = (*(*itr_expression))[current_value_index];
            }
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

    for (IndexType i = 0; i < number_of_triplets; ++i) {
        write_encoded_triplet({get_next_byte(), get_next_byte(), get_next_byte()}, 0);
    }

    const IndexType number_of_bytes_remaining = total_data_size % 3;
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
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    switch (mOutputFormat){
        case ASCII:
            WriteDataElementAscii(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
        case BINARY:
            WriteDataElementBinary(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
    }
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<LiteralFlatExpression<int>::Pointer>& rExpressions,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    switch (mOutputFormat){
        case ASCII:
            WriteDataElementAscii(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
        case BINARY:
            WriteDataElementBinary(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
    }
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
    const std::vector<LiteralFlatExpression<double>::Pointer>& rExpressions,
    const std::vector<IndexType> rNumberOfEntities,
    const IndexType Level)
{
    switch (mOutputFormat){
        case ASCII:
            WriteDataElementAscii(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
        case BINARY:
            WriteDataElementBinary(rTagName, rAttributes, rNumberOfEntities, Level, rExpressions);
            break;
    }
}

} // namespace Kratos