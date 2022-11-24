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
        //std::cout << "Data on " << rCommunicator.Rank() << ": ";
        //for (auto it=itBegin; it!=itEnd; ++it) {
        //    std::cout << *it << " ";
        //}
        //std::cout << std::endl;
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
        rCommunicator.Barrier();
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");

        using Value = typename std::iterator_traits<TInputIterator>::value_type;
        std::vector<Value> output_buffer;

        const int master_rank = 0;
        const int this_rank = rCommunicator.Rank();
        const int number_of_ranks = rCommunicator.Size();

        if (this_rank == master_rank) {
            // Don't bother sendind data to ourselves, and just
            // copy the input to the output buffer
            output_buffer.reserve(std::distance(itBegin, itEnd));
            std::copy(itBegin, itEnd, std::back_inserter(output_buffer));

            std::vector<Value> receive_buffer;
            receive_buffer.reserve(1); // <== we're only ever going to store one received vector here

            for (int i_rank=1; i_rank<number_of_ranks; ++i_rank) {
                // Receive objects from a rank
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
                rCommunicator.Recv(receive_buffer, i_rank, i_rank);
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");

                // Move received objects from the buffer to the output
                output_buffer.reserve(output_buffer.size() + receive_buffer.back().size());
                for (Value& r_item : receive_buffer) {
                    output_buffer.emplace_back(std::move(r_item));
                }
                receive_buffer.clear();
            }
        } else {
            // DataCommunicator operates on objects, or vectors of objects,
            // so that's what we need to pack the input data into
            std::vector<Value> local_objects(itBegin, itEnd);
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
            rCommunicator.Send(local_objects, master_rank, this_rank);
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
        }

        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
        rCommunicator.Broadcast(output_buffer, master_rank);
        KRATOS_LINE_WATCH(rCommunicator.Rank() << ": ");
        for (Value& r_item : output_buffer) {
            *itOutput++ = std::move(r_item);
        }

        rCommunicator.Barrier();
    }
}; // struct MPIUtilities


} // namespace Kratos::HDF5
