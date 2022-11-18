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

// --- HDF5 Includes ---
#include "custom_utilities/mpi_utilities.h"


namespace Kratos {


std::vector<std::string> MPIAllGatherVStrings(const std::vector<std::string>& rLocalStrings,
                                              DataCommunicator& rCommunicator)
{
    std::vector<std::string> output;

    if (rCommunicator.IsDistributed()) {
        const auto rank = rCommunicator.Rank();
        const auto number_of_ranks = rCommunicator.Size();

        // Send strings from slave ranks to the main one
        if (rank) {
            rCommunicator.Send(rLocalStrings, 0, rank);
        } else {
            output = rLocalStrings;
            std::vector<std::string> buffer;

            for (int i_rank=1; i_rank<number_of_ranks; ++i_rank) {
                rCommunicator.Recv(buffer, i_rank, i_rank);
                output.reserve(output.size() + buffer.size());
                for (auto& r_item : buffer) {
                    output.emplace_back(std::move(r_item));
                }

                buffer.clear();
            }
        }

        // Broadcast collected strings from the main rank
        rCommunicator.Broadcast(output, 0);
    }

    return output;
}


} // namespace Kratos
