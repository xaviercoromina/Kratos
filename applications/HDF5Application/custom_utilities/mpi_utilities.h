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
#include "includes/debug_helpers.h"

// --- STL Includes ---
#include <vector>


namespace Kratos::HDF5 {


struct MPIUtilities
{
    /**
     *  @brief Synchronize a union of all items on every rank.
     *  @warning This function is only meant to be used for trivially serializable
     *           value types.
     */
    template <class TInputIterator, class TOutputIterator>
    static void AllGatherV(TInputIterator itBegin,
                           TInputIterator itEnd,
                           TOutputIterator itOutput,
                           DataCommunicator& rCommunicator)
    {
        using Value = typename std::iterator_traits<TInputIterator>::value_type;

        const int this_rank = rCommunicator.Rank();
        const int number_of_ranks = rCommunicator.Size();

        const int send_to = (this_rank + 1) % number_of_ranks;
        const int receive_from = (number_of_ranks + this_rank - 1) % number_of_ranks;

        std::vector<Value> output_buffer(itBegin, itEnd);
        std::vector<Value> communication_buffer = output_buffer;

        for (int i_rank=0; i_rank<number_of_ranks-1; ++i_rank) {
            // Forward new data
            // - if this is the first iteration, the rank-local data is sent
            // - if this is not the first iteration, the last received data is sent
            communication_buffer = rCommunicator.SendRecv(communication_buffer, send_to, receive_from);
            output_buffer.reserve(output_buffer.size() + communication_buffer.size());

            // Copy received data to the output buffer
            std::copy(communication_buffer.begin(), communication_buffer.end(), std::back_inserter(output_buffer));
        }

        // Move items from the output buffer to the output range
        for (Value& r_item : output_buffer) {
            *itOutput++ = std::move(r_item);
        }
    }
}; // struct MPIUtilities


} // namespace Kratos::HDF5
