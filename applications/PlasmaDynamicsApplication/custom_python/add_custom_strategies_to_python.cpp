//
//   Project Name:        KratosPlasmaDynamicsApplication $
//   Last Modified by:    $Author:    Marc Chung To Sang $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//strategies

//builders and solvers

//schemes

//linear solvers


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

}

}  // namespace Python.
} // Namespace Kratos 




