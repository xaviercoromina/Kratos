//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_utilities/response_utilities/mass_response_utilities.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MassResponseUtilities >(m, "MassResponseUtilities")
        .def_static("CalculateMass", &MassResponseUtilities::CalculateMass)
        .def_static("CalculateMassShapeSensitivity", &MassResponseUtilities::CalculateMassShapeSensitivity)
        .def_static("CalculateMassDensitySensitivity", &MassResponseUtilities::CalculateMassDensitySensitivity)
        .def_static("CalculateMassThicknessSensitivity", &MassResponseUtilities::CalculateMassThicknessSensitivity)
        .def_static("CalculateMassCrossAreaSensitivity", &MassResponseUtilities::CalculateMassCrossAreaSensitivity)
        ;

}

}  // namespace Python.
} // Namespace Kratos

