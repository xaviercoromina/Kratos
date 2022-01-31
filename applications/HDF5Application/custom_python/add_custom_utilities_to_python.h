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

#ifndef KRATOS_HDF5_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED
#define KRATOS_HDF5_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED


// External includes
#include "pybind11/pybind11.h"

// Core includes
#include "includes/define.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& rModule);

} // namespace Python
} // namespace Kratos


#endif