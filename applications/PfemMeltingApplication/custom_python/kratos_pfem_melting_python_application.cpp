// KRATOS
// _____   __               __  __      _ _   _
//|  __ \ / _|             |  \/  |    | | | (_)
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "pfem_melting_application.h"
//#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
//#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosPfemMeltingApplication,m)
{

    py::class_<KratosPfemMeltingApplication,
           KratosPfemMeltingApplication::Pointer,
           KratosApplication >(m,"KratosPfemMeltingApplication")
           .def(py::init<>())
           ;
    //AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    //AddCustomResponseFunctionsToPython(m);

    // Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVATION_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ARRHENIUS_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RADIOUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HEAT_OF_VAPORIZATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ARRHENIUS_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALARVELOCITY_X)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALARVELOCITY_Y)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALARVELOCITY_Z)


    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,INITIAL_POSITION)


}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined