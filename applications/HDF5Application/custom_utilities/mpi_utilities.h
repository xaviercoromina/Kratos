//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#pragma once

// --- Core Includes ---
#include "includes/data_communicator.h"


namespace Kratos {


std::vector<std::string> KRATOS_API(HDF5Application) MPIAllGatherVStrings(const std::vector<std::string>& rLocalStrings,
                                                                          DataCommunicator& rCommunicator);


} // namespace Kratos
