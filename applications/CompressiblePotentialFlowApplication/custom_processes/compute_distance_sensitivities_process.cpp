//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez
//

#include "compute_distance_sensitivities_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
// Constructor for ComputeDistanceSensitivities Process
ComputeDistanceSensitivitiesProcess::ComputeDistanceSensitivitiesProcess(ModelPart& rModelPart,
                    Parameters ThisParameters
                ):
    Process(),
    mrModelPart(rModelPart)
{
}

void ComputeDistanceSensitivitiesProcess::Execute()
{
    KRATOS_TRY;


    KRATOS_CATCH("");
}
}// Namespace Kratos
