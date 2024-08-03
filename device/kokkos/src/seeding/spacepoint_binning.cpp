/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/kokkos/seeding/spacepoint_binning.hpp"

#include "traccc/kokkos/utils/definitions.hpp"

// Project include(s).
#include "traccc/seeding/device/count_grid_capacities.hpp"
#include "traccc/seeding/device/populate_grid.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

namespace traccc::kokkos {

spacepoint_binning::spacepoint_binning(
    const seedfinder_config& config, const spacepoint_grid_config& grid_config,
    const traccc::memory_resource& mr, vecmem::copy& copy)
    : m_config(config),
      m_axes(get_axes(grid_config, (mr.host ? *(mr.host) : mr.main))),
      m_mr(mr),
      m_copy(copy) {}

spacepoint_binning::output_type spacepoint_binning::operator()(
    const spacepoint_collection_types::const_view& spacepoints_view) const {

    // Get the spacepoint sizes from the view
    auto sp_size = m_copy.get_size(spacepoints_view);

    if (sp_size == 0) {
        return {m_axes.first, m_axes.second, {}, m_mr.main, m_mr.host};
    }

    // Set up the container that will be filled with the required capacities for
    // the spacepoint grid.
    const std::size_t grid_bins = m_axes.first.n_bins * m_axes.second.n_bins;
    vecmem::data::vector_buffer<unsigned int> grid_capacities_buff(grid_bins,
                                                                   m_mr.main);
    m_copy.setup(grid_capacities_buff);
    m_copy.memset(grid_capacities_buff, 0);
    vecmem::data::vector_view<unsigned int> grid_capacities_view =
        grid_capacities_buff;

    // Calculate the number of threads and thread blocks to run the kernels for.
    const unsigned int num_threads = 32 * 8;
    const unsigned int num_blocks = (sp_size + num_threads - 1) / num_threads;

    // Hack to avoid warnings thrown by C++20
    seedfinder_config config = m_config;
    auto axes = m_axes;

    Kokkos::parallel_for(
        "count_grid_capacities", team_policy(num_blocks, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member_type& team_member) {
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team_member, num_threads),
                [&](const int& thr) {
                    device::count_grid_capacities(
                        team_member.league_rank() * team_member.team_size() +
                            thr,
                        config, axes.first, axes.second, spacepoints_view,
                        grid_capacities_view);
                });
        });

    // Copy grid capacities back to the host
    vecmem::vector<unsigned int> grid_capacities_host(m_mr.host ? m_mr.host
                                                                : &(m_mr.main));
    m_copy(grid_capacities_buff, grid_capacities_host);

    // Create the grid buffer.
    sp_grid_buffer grid_buffer(
        m_axes.first, m_axes.second,
        std::vector<std::size_t>(grid_capacities_host.begin(),
                                 grid_capacities_host.end()),
        m_mr.main, m_mr.host, vecmem::data::buffer_type::resizable);

    m_copy.setup(grid_buffer._buffer);
    sp_grid_view grid_view = grid_buffer;

    // Populate the grid.
    Kokkos::parallel_for(
        "populate_grid", team_policy(num_blocks, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member_type& team_member) {
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team_member, num_threads),
                [&](const int& thr) {
                    device::populate_grid(
                        team_member.league_rank() * team_member.team_size() +
                            thr,
                        config, spacepoints_view, grid_view);
                });
        });
    
    // Return the freshly filled buffer.
    return grid_buffer;
}

}  // namespace traccc::kokkos
