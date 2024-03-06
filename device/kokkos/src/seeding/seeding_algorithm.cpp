/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/kokkos/seeding/seeding_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/detail/seeding_config.hpp"

namespace traccc::kokkos {

seeding_algorithm::seeding_algorithm(const seedfinder_config& finder_config,
                                     const spacepoint_grid_config& grid_config,
                                     const seedfilter_config& filter_config,
                                     const traccc::memory_resource& mr,
                                     vecmem::copy& copy)
    : m_spacepoint_binning(finder_config, grid_config, mr, copy),
      m_seed_finding(finder_config, filter_config, mr, copy) {}

seeding_algorithm::output_type seeding_algorithm::operator()(
    const spacepoint_collection_types::const_view& spacepoints_view) const {

    return m_seed_finding(spacepoints_view,
                          m_spacepoint_binning(spacepoints_view));
}

}  // namespace traccc::kokkos
