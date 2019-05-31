//
//   Project Name:        KratosPlasmaDynamicsApplication $
//   Last Modified by:    $Author:    Marc Chung To Sang $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes

// Project includes
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "../custom_constitutive/DEM_electromagnetic_CL.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    py::class_<DEM_electromagnetic, DEM_electromagnetic::Pointer>(m, "DEM_electromagnetic")
    .def(py::init<>())
    .def("Clone", &DEM_electromagnetic::Clone)
    .def("SetConstitutiveLawInProperties", &DEM_electromagnetic::SetConstitutiveLawInProperties)
    .def("GetTypeOfLaw", &DEM_electromagnetic::GetTypeOfLaw)
    ;
}

}  // namespace Python.
}  // namespace Kratos.
