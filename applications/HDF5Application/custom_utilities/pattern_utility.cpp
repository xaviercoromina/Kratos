//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen, https://github.com/matekelemen
//

// Core includes
#include "ghc/filesystem.hpp"
#include "includes/define.h"

// Project includes
#include "pattern_utility.h"

// STL includes
#include <algorithm>


namespace Kratos
{


PlaceholderPattern::PlaceholderPattern(const std::string& rPattern,
                                       const PlaceholderMap& rPlaceholderMap)
    : mPattern(rPattern),
      mPlaceholderGroupMap(),
      mRegexString(FormatRegexLiteral(rPattern)),
      mRegex()
{
    KRATOS_TRY

    using PositionPair = std::pair<std::size_t,PlaceholderMap::const_iterator>;
    const auto number_of_placeholders = rPlaceholderMap.size();

    // Array for tracking position-placeholder pairs for assigning
    // regex groups to placeholders later.
    // {position_in_pattern, it_placeholder_regex_pair}
    std::vector<PositionPair> position_map;
    position_map.reserve(number_of_placeholders);

    // Replace placeholders with their corresponding regex strings
    // and note their positions in the pattern
    PlaceholderMap::const_iterator it_pair = rPlaceholderMap.begin();
    const auto it_end = rPlaceholderMap.end();

    for ( ; it_pair!=it_end; ++it_pair){
        // Wrap the regex in a group
        std::string regex;
        if (it_pair->second.front()!='(' || it_pair->second.back()!=')') {
            regex = "(" + it_pair->second + ")";
        }
        else {
            regex = it_pair->second;
        }

        const auto placeholder_size = it_pair->first.size();
        const auto regex_size = regex.size();
        const int size_difference = regex_size - placeholder_size;

        while (true) {
            // Find the next instance of the current placeholder
            const auto position_in_pattern = mRegexString.find(it_pair->first);

            // Replace it with its regex (if found)
            if (position_in_pattern != std::string::npos) {
                mRegexString.replace(position_in_pattern, placeholder_size, regex);
            }
            else {
                break;
            }

            // Update positions
            for (auto& r_pair : position_map) {
                if (position_in_pattern < r_pair.first) r_pair.first += size_difference;
            }

            position_map.emplace_back(position_in_pattern, it_pair);
        } // while placeholder in pattern
    } // for placeholder in rPlaceHolderMap

    // Replace positions with indices in ascending order (based on position)
    // in lieu of std::transform_if
    std::sort(
        position_map.begin(),
        position_map.end(),
        [](const PositionPair& rLeft, const PositionPair& rRight) {return rLeft.first < rRight.first;});

    std::size_t index = 0;
    for (auto& r_pair : position_map) r_pair.first = index++;

    // Populate the placeholder - group index map
    it_pair = rPlaceholderMap.begin();
    for ( ; it_pair!=it_end; ++it_pair) {
        // Move the placeholder string and construct an associated empty index array
        auto emplace_result = mPlaceholderGroupMap.emplace(it_pair->first, PlaceholderGroupMap::mapped_type());

        // Fill the index array with the group indices
        for (const auto& r_pair : position_map) {
            if (r_pair.second == it_pair) emplace_result.first->second.push_back(r_pair.first);
        }
    }

    // Construct the regex
    mRegexString = "^" + mRegexString + "$";
    mRegex = std::regex(mRegexString);

    KRATOS_CATCH("");
} // PlaceholderPattern::PlaceholderPattern


bool PlaceholderPattern::IsAMatch(const std::string& rString) const
{
    KRATOS_TRY

    return std::regex_match(rString, mRegex);

    KRATOS_CATCH("");
} // PlaceholderPattern::IsAMatch


PlaceholderPattern::MatchType PlaceholderPattern::Match(const std::string& rString) const
{
    KRATOS_TRY

    std::smatch results;
    MatchType output;

    // Perform regex search and extract matches
    if (std::regex_match(rString, results, mRegex)) {
        for (auto& r_pair : mPlaceholderGroupMap) {
            // Construct empty group matches
            auto emplace_result = output.emplace(r_pair.first, MatchType::value_type::second_type());

            // Collect matches for the current placeholder
            for (auto i_group : r_pair.second) {
                // First match (index 0) is irrelevant because it's the entire pattern,
                // the rest is offset by 1
                const auto i_group_match = i_group + 1;

                if (!results.str(i_group_match).empty()) {
                    emplace_result.first->second.push_back(results.str(i_group_match));
                }
            } // for i_group
        } // for placeholder, group_indices
    } // if regex_match
    else {
        KRATOS_ERROR << "'" << rString << "' is not a match for '" << this->GetRegexString() << "'";
    }

    return output;

    KRATOS_CATCH("");
}


std::string PlaceholderPattern::Apply(const PlaceholderMap& rPlaceholderValueMap) const
{
    KRATOS_TRY

    auto output = mPattern;

    for (const auto& r_pair : rPlaceholderValueMap) {
        while (true) {
            auto position = output.find(r_pair.first);
            if (position != output.npos) output.replace(position, r_pair.second.size(), r_pair.second);
            else break;
        } // while placeholder in output
    } // for placeholder, value in map

    return output;

    KRATOS_CATCH("");
}


std::vector<PlaceholderPattern::PathType> PlaceholderPattern::Glob() const
{
    KRATOS_TRY

    PathType pattern(mRegexString);
    std::vector<PathType> paths;

    // Decide where to begin globbing
    if (pattern.is_absolute()) paths.emplace_back(pattern.root_name());
    else paths.emplace_back(ghc::filesystem::current_path());

    // Compare the pattern parts to the globbed files'/directories' parts
    for (auto pattern_part : pattern) {
        std::regex pattern_part_regex(pattern_part.string());
        std::vector<PathType> tmp_paths;

        for (const auto& r_path : paths) {
            if (ghc::filesystem::is_directory(r_path)) {
                for (auto item : ghc::filesystem::directory_iterator(r_path)) {
                    if (std::regex_match(item.path().filename().string(), pattern_part_regex)) {
                        tmp_paths.emplace_back(item);
                    }
                } // for r_item in path.glob(*)
            } // if r_path.is_directory()
        } // for r_path in paths

        paths = tmp_paths;
    } // for pattern_part in pattern.parts

    return paths;

    KRATOS_CATCH("");
}


const std::regex& PlaceholderPattern::GetRegex() const
{
    return mRegex;
}


const std::string& PlaceholderPattern::GetRegexString() const
{
    return mRegexString;
}


std::string PlaceholderPattern::FormatRegexLiteral(const std::string& rLiteral)
{
    KRATOS_TRY

    auto output = rLiteral;

    for (char char_to_escape : R"(!$()*+-?[\]^)") {
        std::size_t position = 0;

        std::string escaped;
        escaped.reserve(2);
        escaped.push_back('\\');
        escaped.push_back(char_to_escape);

        while (true) {
            if (output.size() <= position) break;

            position = output.find(char_to_escape, position);
            if (position == output.npos) {
                break;
            }
            else {
                // Escape the sensitive character
                if (!position || output[position-1] != '\\') {
                    output.replace(position, 1, escaped);
                    position += 2;
                }
            } // if char_to_escape in output
        } // while True
    } // for char_to_escape

    return output;

    KRATOS_CATCH("");
}


const ModelPartPattern::PlaceholderMap ModelPartPattern::mModelPartPlaceholderMap = {
    {"<model_part_name>", ".+"},
    {"<step>", "[0-9]+"},
    {"<time>", "[+-]?[0-9]*[.]?[0-9]+"}
};


ModelPartPattern::ModelPartPattern(const std::string& rPattern)
    : PlaceholderPattern(rPattern, mModelPartPlaceholderMap)
{
}


} // namespace Kratos