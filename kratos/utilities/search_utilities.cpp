//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/search_utilities.h"

namespace Kratos {

SearchUtilities::BoundingBoxType SearchUtilities::ComputeGlobalBoundingBox(const ModelPart& rModelPart)
{
    BoundingBoxType local_bounding_box = ComputeLocalBoundingBox(rModelPart);

    array_1d<double,3> max_vals, min_vals;

    // Fill buffers for MPI
    for (std::size_t i=0; i<3; ++i) {
        max_vals[i] = local_bounding_box[i*2];
        min_vals[i] = local_bounding_box[i*2+1];
    }

    // Compute global values
    const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
    if (r_data_comm.IsDefinedOnThisRank()) {
        max_vals = r_data_comm.MaxAll(max_vals);
        min_vals = r_data_comm.MinAll(min_vals);
    }

    BoundingBoxType global_bounding_box;
    // Extract information from buffers
    for (std::size_t i=0; i<3; ++i) {
        global_bounding_box[i*2] = max_vals[i];
        global_bounding_box[i*2+1] = min_vals[i];
    }

    return global_bounding_box;
}

/***********************************************************************************/
/***********************************************************************************/

SearchUtilities::BoundingBoxType SearchUtilities::ComputeLocalBoundingBox(const ModelPart& rModelPart)
{
    BoundingBoxType local_bounding_box {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    // loop over all nodes (local and ghost(necessary if conditions have only ghost nodes) )
    for (auto &r_node : rModelPart.Nodes()) {
        local_bounding_box[0] = std::max(r_node.X(), local_bounding_box[0]);
        local_bounding_box[1] = std::min(r_node.X(), local_bounding_box[1]);
        local_bounding_box[2] = std::max(r_node.Y(), local_bounding_box[2]);
        local_bounding_box[3] = std::min(r_node.Y(), local_bounding_box[3]);
        local_bounding_box[4] = std::max(r_node.Z(), local_bounding_box[4]);
        local_bounding_box[5] = std::min(r_node.Z(), local_bounding_box[5]);
    }
    return local_bounding_box;
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::ComputeBoundingBoxesWithTolerance(
    const std::vector<double>& rBoundingBoxes,
    const double Tolerance,
    std::vector<double>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances
    for (IndexType i=0; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;
    }

    for (IndexType i=1; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::string SearchUtilities::BoundingBoxStringStream(const BoundingBoxType& rBoundingBox)
{
    // xmax, xmin,  ymax, ymin,  zmax, zmin
    std::stringstream buffer;
    buffer << "[" << rBoundingBox[1] << " "    // xmin
                  << rBoundingBox[3] << " "    // ymin
                  << rBoundingBox[5] << "]|["  // zmin
                  << rBoundingBox[0] << " "    // xmax
                  << rBoundingBox[2] << " "    // ymax
                  << rBoundingBox[4] << "]";   // zmax
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

bool SearchUtilities::PointIsInsideBoundingBox(
    const BoundingBoxType& rBoundingBox,
    const array_1d<double, 3>& rCoords
    )
{   // The Bounding Box should have some tolerance already!
    if (rCoords[0] < rBoundingBox[0] && rCoords[0] > rBoundingBox[1])   // check x-direction
        if (rCoords[1] < rBoundingBox[2] && rCoords[1] > rBoundingBox[3])   // check y-direction
            if (rCoords[2] < rBoundingBox[4] && rCoords[2] > rBoundingBox[5])   // check z-direction
                return true;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::FillBufferBeforeLocalSearch(
    const SearchLocalSystemPointerVector& rSearchLocalSystems,
    const std::vector<double>& rBoundingBoxes,
    const SizeType BufferSizeEstimate,
    std::vector<std::vector<double>>& rSendBuffer,
    std::vector<int>& rSendSizes
    )
{
    const SizeType comm_size = rSendBuffer.size();

    BoundingBoxType bounding_box;
    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(rBoundingBoxes.size(), 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        auto& r_rank_buffer = rSendBuffer[i_rank];
        r_rank_buffer.clear();
        r_rank_buffer.reserve(BufferSizeEstimate);
        rSendSizes[i_rank] = 0;

        for (IndexType i_local_sys=0; i_local_sys<rSearchLocalSystems.size(); ++i_local_sys) {
            for (IndexType j=0; j<6; ++j) {
                bounding_box[j] = rBoundingBoxes[(i_rank*6) + j]; // retrieve bounding box of partition
            }

            const auto& rp_local_sys = rSearchLocalSystems[i_local_sys];

            if (!rp_local_sys->IsDoneSearching()) {
                const auto& r_coords = rp_local_sys->Coordinates();
                if (SearchUtilities::PointIsInsideBoundingBox(bounding_box, r_coords)) {
                    // These push_backs are threadsafe bcs only one vector is accessed per thread!
                    r_rank_buffer.push_back(static_cast<double>(i_local_sys)); // this it the "mSourceLocalSystemIndex" of the SearchInterfaceInfo
                    r_rank_buffer.push_back(r_coords[0]);
                    r_rank_buffer.push_back(r_coords[1]);
                    r_rank_buffer.push_back(r_coords[2]);

                    rSendSizes[i_rank] += 4;
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::CreateSearchInterfaceInfosFromBuffer(
    const std::vector<std::vector<double>>& rRecvBuffer,
    const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer
    )
{
    const SizeType comm_size = rSearchInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rRecvBuffer.size() == comm_size) << "Buffer-size mismatch!" << std::endl;

    array_1d<double, 3> coords;
    // Loop the ranks and construct the SearchInterfaceInfos
    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        const SizeType recv_buffer_size_rank = rRecvBuffer[i_rank].size();
        KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(recv_buffer_size_rank, 4) == 0) << "Rank " << CommRank
            << " received a wrong buffer-size from rank " << i_rank+1 << "!" << std::endl;
        const SizeType num_objs = recv_buffer_size_rank / 4; // 1 index and 3 coordinates

        const auto& r_rank_buffer = rRecvBuffer[i_rank];
        auto& r_interface_infos_rank = rSearchInterfaceInfosContainer[i_rank];
        r_interface_infos_rank.clear();

        if (r_interface_infos_rank.size() != num_objs) {
            r_interface_infos_rank.resize(num_objs);
        }

        for (IndexType j=0; j<num_objs; ++j) {
#ifdef KRATOS_DEBUG
            // with this check we make sure that this field
            // only contains doubles converted from ints
            // e.g. 4.5 is not allowed!
            double int_part;
            double fract_part = std::modf((r_rank_buffer[j*4]+0.1), &int_part);

            KRATOS_ERROR_IF(std::abs(fract_part-0.1) > 1e-10)
                << "Buffer contains a double (" << r_rank_buffer[j*4]
                << ") that was not casted from an int, i.e. it contains a "
                << "fractional part of " << std::abs(fract_part-0.1) << "!" << std::endl;
#endif
            // retrive data from buffer
            const int local_sys_idx = static_cast<IndexType>(r_rank_buffer[j*4]+0.1);
            // 0.1 is added to prevent truncation errors like (int)1.9999 = 1
            coords[0] = r_rank_buffer[j*4 + 1];
            coords[1] = r_rank_buffer[j*4 + 2];
            coords[2] = r_rank_buffer[j*4 + 3];
            r_interface_infos_rank[j] = rpRefInterfaceInfo->Create(coords, local_sys_idx, i_rank);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::FillBufferAfterLocalSearch(
    SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer,
    const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    std::vector<std::vector<char>>& rSendBuffer,
    std::vector<int>& rSendSizes
    )
{
    const SizeType comm_size = rSearchInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rSendSizes.size() == comm_size) << "Buffer-size mismatch!" << std::endl;

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        if (i_rank != static_cast<IndexType>(CommRank)) {
            SearchInterfaceInfoSerializer interface_infos_serializer(rSearchInterfaceInfosContainer[i_rank], rpRefInterfaceInfo );

            Kratos::StreamSerializer serializer;
            serializer.save("interface_infos", interface_infos_serializer);

            const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
            const std::string& stream_str = p_serializer_buffer->str();

            const SizeType send_size = sizeof(char) * (stream_str.size()+1); // +1 fof Null-terminated string

            rSendSizes[i_rank] = send_size;

            auto& r_rank_buffer = rSendBuffer[i_rank];
            r_rank_buffer.clear();

            if (r_rank_buffer.size() != send_size) {
                r_rank_buffer.resize(send_size);
            }

            std::memcpy(r_rank_buffer.data(), stream_str.c_str(), send_size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::DeserializeSearchInterfaceInfosFromBuffer(
    const std::vector<std::vector<char>>& rRecvBuffer,
    const SearchInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    SearchInterfaceInfoPointerVectorType& rSearchInterfaceInfosContainer
    )
{
    const SizeType comm_size = rSearchInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rRecvBuffer.size() == comm_size)
        << "Buffer-size mismatch!" << std::endl;

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        if (i_rank != static_cast<IndexType>(CommRank)) {
            Kratos::StreamSerializer serializer;

            const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
            p_serializer_buffer->write(rRecvBuffer[i_rank].data(), rRecvBuffer[i_rank].size());

            SearchInterfaceInfoSerializer interface_infos_serializer(rSearchInterfaceInfosContainer[i_rank], rpRefInterfaceInfo );

            serializer.load("interface_infos", interface_infos_serializer);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchInterfaceInfoSerializer::save(Kratos::Serializer& rSerializer) const
{
    const SizeType num_infos = mrInterfaceInfos.size();
    rSerializer.save("size", num_infos);

    for (IndexType i=0; i<num_infos; ++i) {
        rSerializer.save("E", *(mrInterfaceInfos[i])); // NOT serializing the shared_ptr!
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchInterfaceInfoSerializer::load(Kratos::Serializer& rSerializer)
{
    mrInterfaceInfos.clear(); // make sure it has no leftovers

    SizeType num_infos;
    rSerializer.load("size", num_infos);

    if (mrInterfaceInfos.size() != num_infos) {
        mrInterfaceInfos.resize(num_infos);
    }

    for (IndexType i=0; i<num_infos; ++i) {
        // first we create a new object, then we load its data
        // this is needed bcs of the polymorphic behavior of the InterfaceInfos
        // i.e. in order to create the correct type
        // => the vector contains baseclass-pointers!
        // Jordi I am quite sure that I could get around it by registering it, what do you think... TODO
        // I think doing it manually is more efficient, which I want so I would probably leave it ...
        // The serializer does some nasty casting when pointers are serialized...
        // I could do a benchmark at some point but I highly doubt that the serializer is faster ...
        mrInterfaceInfos[i] = mrpRefInterfaceInfo->Create();
        rSerializer.load("E", *(mrInterfaceInfos[i])); // NOT serializing the shared_ptr!
    }
}

}