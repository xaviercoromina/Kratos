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
#include "pybind11/detail/common.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"
#include "pybind11/stl/filesystem.h"

// HDF5 includes
#include "custom_utilities/vertex.h"
#include "custom_utilities/vertex_utilities.h"
#include "custom_utilities/pattern_utility.h"
#include "custom_utilities/registry_file.h"

// Internal includes
#include "add_custom_utilities_to_python.h"

// STL includes
#include <iterator>
#include <algorithm>


namespace Kratos
{
namespace Python
{


namespace
{
class PointLocatorAdaptorTrampoline : public HDF5::PointLocatorAdaptor
{
public:
    using HDF5::PointLocatorAdaptor::PointLocatorAdaptor;

    const Element::WeakPointer FindElement(const Point& rPoint) const override
    {
        using ReturnType = const Element::WeakPointer;
        using BaseType = HDF5::PointLocatorAdaptor;

        PYBIND11_OVERRIDE_PURE(
            ReturnType,
            BaseType,
            FindElement,
            rPoint);
    }
}; // class PointLocatorAdaptorTrampoline

/// Convert an array of paths to an array of strings
std::vector<std::string> Glob (const ModelPartPattern& rInstance) {
    std::vector<std::string> output;
    auto result = rInstance.Glob();
    std::transform(result.begin(),
                   result.end(),
                   std::back_inserter(output),
                   [](const ModelPartPattern::PathType& rItem) {return rItem.string();});
    return output;
}
} // namespace


void AddCustomUtilitiesToPython(pybind11::module& rModule)
{
    pybind11::class_<HDF5::PointLocatorAdaptor, HDF5::PointLocatorAdaptor::Pointer, PointLocatorAdaptorTrampoline>(rModule, "PointLocatorAdaptor")
        .def(pybind11::init<>())
        .def("FindElement", &HDF5::PointLocatorAdaptor::FindElement)
        ;

    pybind11::class_<HDF5::BruteForcePointLocatorAdaptor, HDF5::BruteForcePointLocatorAdaptor::Pointer, HDF5::PointLocatorAdaptor>(rModule, "BruteForcePointLocatorAdaptor")
        .def(pybind11::init<ModelPart&, const Globals::Configuration, const double>())
        ;

    #define KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(TValue) \
        .def("GetValue", [](const HDF5::Detail::Vertex& rVertex, const Variable<TValue>& rVariable) {return rVertex.GetValue(rVariable);})

    using Array3 = array_1d<double,3>;
    using Array4 = array_1d<double,4>;
    using Array6 = array_1d<double,6>;
    using Array9 = array_1d<double,9>;

    pybind11::class_<HDF5::Detail::Vertex, HDF5::Detail::Vertex::Pointer, Point>(rModule, "Vertex")
        .def(pybind11::init<const array_1d<double,3>&, const HDF5::PointLocatorAdaptor&, bool>())
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(bool)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(int)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(double)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Vector)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array3)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array4)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array6)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array9)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Matrix)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(DenseVector<int>)
        .def_static("MakeShared", &HDF5::Detail::Vertex::MakeShared)
        .def("IsLocated", &HDF5::Detail::Vertex::IsLocated)
        .def("GetID", &HDF5::Detail::Vertex::GetID)
        ;
    #undef KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING

    pybind11::class_<HDF5::Detail::VertexContainerType, HDF5::Detail::VertexContainerType::Pointer>(rModule, "VertexContainer")
        .def(pybind11::init<>())
        .def("push_back", &HDF5::Detail::VertexContainerType::push_back)
        ;

    pybind11::class_<PlaceholderPattern, PlaceholderPattern::Pointer>(rModule, "PlaceholderPattern")
        .def(pybind11::init<const std::string&,const PlaceholderPattern::PlaceholderMap&>())
        .def("IsAMatch",
             &PlaceholderPattern::IsAMatch,
             "Check whether a string satisfies the pattern")
        .def("Match",
             &PlaceholderPattern::Match,
             "Find all placeholders' values in the input string.")
        .def("Apply",
             &PlaceholderPattern::Apply,
             "Substitute values from the input map into the stored pattern.")
        .def("GetRegexString",
             &PlaceholderPattern::GetRegexString,
             "Get the string representation of the regex.")
        ;

    pybind11::class_<ModelPartPattern, ModelPartPattern::Pointer, PlaceholderPattern>(rModule, "ModelPartPattern")
        .def(pybind11::init<const std::string&>())
        .def("Glob",
             &Glob,
             "Collect all file/directory paths that match the pattern.")
        .def("Apply",
             static_cast<std::string(ModelPartPattern::*)(const ModelPartPattern::PlaceholderMap&)const>(&ModelPartPattern::Apply),
             "Substitute values from the input map into the stored pattern.")
        .def("Apply",
             static_cast<std::string(ModelPartPattern::*)(const ModelPart&)const>(&ModelPartPattern::Apply),
             "Substitute values from the model part into the stored pattern.")
        ;

    pybind11::class_<CheckpointPattern, CheckpointPattern::Pointer, ModelPartPattern>(rModule, "CheckpointPattern")
        .def(pybind11::init<const std::string&>())
        .def("Apply",
             static_cast<std::string(CheckpointPattern::*)(const CheckpointPattern::PlaceholderMap&)const>(&CheckpointPattern::Apply),
             "Substitute values from the input map into the stored pattern.")
        .def("Apply",
             static_cast<std::string(CheckpointPattern::*)(const ModelPart&,std::size_t)const>(&CheckpointPattern::Apply),
             "Substitute values from the provided model part and path ID into the stored pattern.")
        ;

    pybind11::class_<RegistryFile, RegistryFile::Pointer>(rModule, "RegistryFile")
        .def(pybind11::init<const std::filesystem::path&>())
        .def(pybind11::init<const std::filesystem::path&,const RegistryFile::Extractor&>())
        .def("GetFilePath", &RegistryFile::GetFilePath, "Get the path to the underlying file.")
        .def("SetExtractor", pybind11::overload_cast<const RegistryFile::Extractor&>(&RegistryFile::SetExtractor))
        .def("Push", &RegistryFile::Push, "Insert a new entry at the end, extracted from the input model.")
        .def("Clear", &RegistryFile::Clear, "Delete the registry file")
        .def("__len__", &RegistryFile::size)
        .def("__iter__", [](const RegistryFile& rRegistryFile){return pybind11::make_iterator(rRegistryFile.begin(), rRegistryFile.end());})
        ;
}


} // namespace Python
} // namespace Kratos