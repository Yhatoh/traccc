/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
// #include "../utils/utils.hpp"
#include "traccc/kokkos/seeding/seed_finding.hpp"

#include "traccc/kokkos/utils/definitions.hpp"

// Project include(s).
#include "traccc/device/fill_prefix_sum.hpp"
#include "traccc/device/make_prefix_sum_buffer.hpp"
#include "traccc/edm/device/device_doublet.hpp"
#include "traccc/edm/device/device_triplet.hpp"
#include "traccc/edm/device/doublet_counter.hpp"
#include "traccc/edm/device/seeding_global_counter.hpp"
#include "traccc/edm/device/triplet_counter.hpp"
#include "traccc/kokkos/utils/make_prefix_sum_buff.hpp"
#include "traccc/seeding/device/count_doublets.hpp"
#include "traccc/seeding/device/count_triplets.hpp"
#include "traccc/seeding/device/find_doublets.hpp"
#include "traccc/seeding/device/find_triplets.hpp"
#include "traccc/seeding/device/reduce_triplet_counts.hpp"
#include "traccc/seeding/device/select_seeds.hpp"
#include "traccc/seeding/device/update_triplet_weights.hpp"

// System include(s).
#include <algorithm>
#include <vector>

namespace traccc::kokkos {

seed_finding::seed_finding(const seedfinder_config& config,
                           const seedfilter_config& filter_config,
                           const traccc::memory_resource& mr,
                           vecmem::copy& copy)
    : m_seedfinder_config(config.toInternalUnits()),
      m_seedfilter_config(filter_config.toInternalUnits()),
      m_mr(mr),
      m_copy(copy) {}

seed_finding::output_type seed_finding::operator()(
    const spacepoint_collection_types::const_view& spacepoints_view,
    const sp_grid_const_view& g2_view) const {

    // Get the sizes from the grid view
    auto grid_sizes = m_copy.get_sizes(g2_view._data_view);

    // Create prefix sum buffer
    vecmem::data::vector_buffer sp_grid_prefix_sum_buff =
        make_prefix_sum_buff(grid_sizes, m_copy, m_mr);

    // Set up the doublet counter buffer.
    device::doublet_counter_collection_types::buffer doublet_counter_buffer = {
        m_copy.get_size(sp_grid_prefix_sum_buff), m_mr.main,
        vecmem::data::buffer_type::resizable};
    m_copy.setup(doublet_counter_buffer);

    // Calculate the number of threads and thread blocks to run the doublet
    // counting kernel for.

    // Counter for the total number of doublets and triplets
    vecmem::unique_alloc_ptr<device::seeding_global_counter>
        globalCounter_device =
            vecmem::make_unique_alloc<device::seeding_global_counter>(
                m_mr.main);
    auto globalCounter_device_ptr = globalCounter_device.get();
    Kokkos::parallel_for(
        "memset_globalcounter", 1, KOKKOS_LAMBDA(const uint64_t i) {
            (globalCounter_device_ptr + i)->m_nMidBot = 0;
            (globalCounter_device_ptr + i)->m_nMidTop = 0;
            (globalCounter_device_ptr + i)->m_nTriplets = 0;
        });

    // Count the number of doublets that we need to produce.
    {
        sp_grid_const_view sp_grid = g2_view;
        vecmem::data::vector_view<const device::prefix_sum_element_t>
            sp_prefix_sum = sp_grid_prefix_sum_buff;
        device::doublet_counter_collection_types::view doublet_counter =
            doublet_counter_buffer;
        unsigned int* nMidBot = &((*globalCounter_device).m_nMidBot);
        unsigned int* nMidTop = &((*globalCounter_device).m_nMidTop);
        Kokkos::parallel_for(
            "count_doublets", m_copy.get_size(sp_grid_prefix_sum_buff),
            KOKKOS_LAMBDA(const uint64_t i) {
                device::count_doublets(i, m_seedfinder_config, sp_grid,
                                       sp_prefix_sum, doublet_counter, *nMidBot,
                                       *nMidTop);
            });
    }

    // Get the summary values.
    vecmem::unique_alloc_ptr<device::seeding_global_counter>
        globalCounter_host =
            vecmem::make_unique_alloc<device::seeding_global_counter>(
                (m_mr.host != nullptr) ? *(m_mr.host) : m_mr.main);

    // Set up the doublet counter buffers.
    device::device_doublet_collection_types::buffer doublet_buffer_mb = {
        globalCounter_host->m_nMidBot, m_mr.main};
    m_copy.setup(doublet_buffer_mb);
    device::device_doublet_collection_types::buffer doublet_buffer_mt = {
        globalCounter_host->m_nMidTop, m_mr.main};
    m_copy.setup(doublet_buffer_mt);

    // Calculate the number of threads and thread blocks to run the doublet
    // finding kernel for.
    const unsigned int doublet_counter_buffer_size =
        m_copy.get_size(doublet_counter_buffer);

    // Find all of the spacepoint doublets.
    {
        sp_grid_const_view sp_grid = g2_view;
        device::doublet_counter_collection_types::const_view doublet_counter =
            doublet_counter_buffer;
        device::device_doublet_collection_types::view mb_doublets =
            doublet_buffer_mb;
        device::device_doublet_collection_types::view mt_doublets =
            doublet_buffer_mt;
        Kokkos::parallel_for(
            "find_doublets", doublet_counter_buffer_size,
            KOKKOS_LAMBDA(const uint64_t i) {
                std::cout << "try" << std::endl;
                device::find_doublets(i, m_seedfinder_config, sp_grid,
                                      doublet_counter, mb_doublets,
                                      mt_doublets);
            });
    }
    // Set up the triplet counter buffers
    device::triplet_counter_spM_collection_types::buffer
        triplet_counter_spM_buffer = {doublet_counter_buffer_size, m_mr.main};
    m_copy.setup(triplet_counter_spM_buffer);
    m_copy.memset(triplet_counter_spM_buffer, 0);
    device::triplet_counter_collection_types::buffer
        triplet_counter_midBot_buffer = {globalCounter_host->m_nMidBot,
                                         m_mr.main,
                                         vecmem::data::buffer_type::resizable};
    m_copy.setup(triplet_counter_midBot_buffer);

    // Calculate the number of threads and thread blocks to run the doublet
    // counting kernel for.

    // Count the number of triplets that we need to produce.
    {
        sp_grid_const_view sp_grid = g2_view;
        device::doublet_counter_collection_types::const_view doublet_counter =
            doublet_counter_buffer;
        device::device_doublet_collection_types::const_view mb_doublets =
            doublet_buffer_mb;
        device::device_doublet_collection_types::const_view mt_doublets =
            doublet_buffer_mt;
        device::triplet_counter_spM_collection_types::view spM_counter =
            triplet_counter_spM_buffer;
        device::triplet_counter_collection_types::view midBot_counter =
            triplet_counter_midBot_buffer;
        Kokkos::parallel_for(
            "count_triplets", globalCounter_host->m_nMidBot,
            KOKKOS_LAMBDA(const uint64_t i) {
                device::count_triplets(
                    i, m_seedfinder_config, sp_grid, doublet_counter,
                    mb_doublets, mt_doublets, spM_counter, midBot_counter);
            });
    }
    // Calculate the number of threads and thread blocks to run the triplet
    // count reduction kernel for.

    // Reduce the triplet counts per spM.
    {
        device::doublet_counter_collection_types::const_view doublet_counter =
            doublet_counter_buffer;
        device::triplet_counter_spM_collection_types::view spM_counter =
            triplet_counter_spM_buffer;
        unsigned int num_triplets = (*globalCounter_device).m_nTriplets;
        Kokkos::parallel_for(
            "reduce_triplet_counts", doublet_counter_buffer_size,
            KOKKOS_LAMBDA(const uint64_t i) {
                unsigned int aux_num_triplets = num_triplets;
                device::reduce_triplet_counts(i, doublet_counter, spM_counter,
                                              aux_num_triplets);
            });
        (*globalCounter_device).m_nTriplets = num_triplets;
    }

    // Set up the triplet buffer.
    device::device_triplet_collection_types::buffer triplet_buffer = {
        globalCounter_host->m_nTriplets, m_mr.main};
    m_copy.setup(triplet_buffer);

    // Calculate the number of threads and thread blocks to run the triplet
    // Find all of the spacepoint triplets.
    {
        sp_grid_const_view sp_grid = g2_view;
        device::doublet_counter_collection_types::const_view doublet_counter =
            doublet_counter_buffer;
        device::device_doublet_collection_types::const_view mt_doublets =
            doublet_buffer_mt;
        device::triplet_counter_spM_collection_types::const_view spM_tc =
            triplet_counter_spM_buffer;
        device::triplet_counter_collection_types::const_view midBot_tc =
            triplet_counter_midBot_buffer;
        device::device_triplet_collection_types::view triplet_view =
            triplet_buffer;

        Kokkos::parallel_for(
            "find_triplets", m_copy.get_size(triplet_counter_midBot_buffer),
            KOKKOS_LAMBDA(const uint64_t i) {
                device::find_triplets(i, m_seedfinder_config,
                                      m_seedfilter_config, sp_grid,
                                      doublet_counter, mt_doublets, spM_tc,
                                      midBot_tc, triplet_view);
            });
    }
    // Calculate the number of threads and thread blocks to run the weight
    // updating kernel for.
    // Array for temporary storage of quality parameters for comparing triplets
    // within weight updating kernel
    scalar* data = (scalar*)m_mr.main.allocate(
        sizeof(scalar) * m_seedfilter_config.compatSeedLimit * WARP_SIZE * 2);

    // Update the weights of all spacepoint triplets.
    {
        sp_grid_const_view sp_grid = g2_view;
        device::triplet_counter_spM_collection_types::const_view spM_tc =
            triplet_counter_spM_buffer;
        device::triplet_counter_collection_types::const_view midBot_tc =
            triplet_counter_midBot_buffer;
        device::device_triplet_collection_types::view triplet_view =
            triplet_buffer;
        Kokkos::parallel_for(
            "update_triplet_weights", globalCounter_host->m_nTriplets,
            KOKKOS_LAMBDA(const uint64_t i) {
                scalar* dataPos = &data[(i % (WARP_SIZE * 2)) *
                                        m_seedfilter_config.compatSeedLimit];
                device::update_triplet_weights(i, m_seedfilter_config, sp_grid,
                                               spM_tc, midBot_tc, dataPos,
                                               triplet_view);
            });
    }
    // Create result object: collection of seeds
    seed_collection_types::buffer seed_buffer(
        globalCounter_host->m_nTriplets, m_mr.main,
        vecmem::data::buffer_type::resizable);
    m_copy.setup(seed_buffer);

    // Calculate the number of threads and thread blocks to run the seed
    // selecting kernel for.
    // Array for temporary storage of triplets for comparing within seed
    // selecting triplet
    triplet* data2 = (triplet*)m_mr.main.allocate(
        sizeof(triplet) * m_seedfilter_config.max_triplets_per_spM * WARP_SIZE *
        2);

    // Each thread uses max_triplets_per_spM elements of the array'
    {
        sp_grid_const_view internal_sp_view = g2_view;
        device::triplet_counter_spM_collection_types::const_view spM_tc =
            triplet_counter_spM_buffer;
        device::triplet_counter_collection_types::const_view midBot_tc =
            triplet_counter_midBot_buffer;
        device::device_triplet_collection_types::view triplet_view =
            triplet_buffer;
        seed_collection_types::view seed_view = seed_buffer;

        Kokkos::parallel_for(
            "select_seeds", doublet_counter_buffer_size,
            KOKKOS_LAMBDA(const uint64_t i) {
                // Each thread uses compatSeedLimit elements of the
                // array
                triplet* dataPos =
                    &data2[(i % (WARP_SIZE * 2)) *
                           m_seedfilter_config.max_triplets_per_spM];
                device::select_seeds(i, m_seedfilter_config, spacepoints_view,
                                     internal_sp_view, spM_tc, midBot_tc,
                                     triplet_view, dataPos, seed_view);
            });
    }

    return seed_buffer;
}

}  // namespace traccc::kokkos
