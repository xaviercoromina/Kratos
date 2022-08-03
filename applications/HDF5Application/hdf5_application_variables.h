//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_HDF5_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_HDF5_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{

/// Additional identifier for solution steps to help keeping track of them during checkpointing.
KRATOS_DEFINE_APPLICATION_VARIABLE(HDF5_APPLICATION, int, ANALYSIS_PATH);

} // namespace Kratos

#endif	/* KRATOS_HDF5_APPLICATION_VARIABLES_H_INCLUDED */
