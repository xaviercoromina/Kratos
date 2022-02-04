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

// External includes
#include "pybind11/stl.h"

// HDF5 includes
#include "custom_utilities/pattern_utility.h"

// Internal includes
#include "add_custom_utilities_to_python.h"

// STL includes
#include <iterator>
#include <algorithm>


namespace Kratos {
namespace Python {


namespace {
/// Convert an array of paths to an array of strings
std::vector<std::string> Glob (const PlaceholderPattern& rInstance) {
    std::vector<std::string> output;
    auto result = rInstance.Glob();
    std::transform(result.begin(),
                   result.end(),
                   std::back_inserter(output),
                   [](const PlaceholderPattern::PathType& rItem) {return rItem.string();});
    return output;
}
} // namespace


void AddCustomUtilitiesToPython(pybind11::module &rModule)
{
    namespace py = pybind11;

    py::class_<PlaceholderPattern, PlaceholderPattern::Pointer>(
        rModule, "PlaceholderPattern")
        .def(py::init<const std::string&,const PlaceholderPattern::PlaceholderMap&>())
        .def("IsAMatch", &PlaceholderPattern::IsAMatch)
        .def("Match", &PlaceholderPattern::Match)
        .def("Apply", &PlaceholderPattern::Apply)
        .def("Glob", &Glob)
        .def("GetRegexString", &PlaceholderPattern::GetRegexString)
        ;

    py::class_<ModelPartPattern, ModelPartPattern::Pointer, PlaceholderPattern>(
        rModule, "ModelPartPattern")
        .def(py::init<const std::string&>())
        ;
}

} // namespace Python
} // namespace Kratos