//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "utilities/search_utilities.h"

namespace Kratos::Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(SearchUtilities_ComputeGlobalBoundingBox_distributed, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    ModelPartCommunicatorUtilities::SetMPICommunicator(model_part);

    KRATOS_ERROR_IF_NOT(model_part.IsDistributed()) << "ModelPart setup failed!" << std::endl;

    const int my_pid = model_part.GetCommunicator().MyPID();
    const int total_procs = model_part.GetCommunicator().TotalProcesses();

    model_part.CreateNewNode(1, 0.2, 5.3, -8.3);
    model_part.CreateNewNode(2, 8.2+my_pid, 25.3, 16.4);
    model_part.CreateNewNode(3, -9.2, -17.13, 1.5);
    model_part.CreateNewNode(4, 12.6+my_pid, 5.3, -8.3-my_pid);

    const auto bbox = SearchUtilities::ComputeGlobalBoundingBox(model_part);

    // std::cout << SearchUtilities::BoundingBoxStringStream(bbox) << std::endl;

    const SearchUtilities::BoundingBoxType exp_bbox = {
        12.6 + total_procs -1,
        -9.2,
        25.3,
        -17.13,
        16.4,
        -8.3-total_procs+1
    };

    KRATOS_CHECK_EQUAL(bbox.size(), 6);
    for (std::size_t i=0; i<6; ++i){
        KRATOS_CHECK_NEAR(bbox[i], exp_bbox[i], 1e-12);
    }
}

}  // namespace Kratos::Testing