/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_parameters.hpp"

namespace traccc::device {

/// Function for combinatorial finding.
/// If the chi2 of the measurement < chi2_max, its measurement index and the
/// index of the link from the previous step are added to the link container.
///
/// @param[in] globalIndex        The index of the current thread
/// @param[in] cfg                Track finding config object
/// @param[in] det_data           Detector view object
/// @param[in] measurements_view  Measurements container view
<<<<<<< HEAD
/// @param[in] barcodes_view      Barcode sequence view object
=======
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
/// @param[in] upper_bounds_view  Upper bounds of measurements unique w.r.t
/// barcode
/// @param[in] in_params_view     Input parameters
/// @param[in] n_measurements_prefix_sum_view  Prefix sum of the number of
/// measurements per parameter
/// @param[in] ref_meas_idx_view  The first index of measurements per parameter
/// @param[in] step               Step index
/// @param[in] n_max_candidates   Number of maximum candidates
/// @param[out] out_params_view   Output parameters
/// @param[out] links_view        link container for the current step
/// @param[out] n_candidates      The number of candidates for the current step
///
template <typename detector_t, typename config_t>
TRACCC_DEVICE inline void find_tracks(
    std::size_t globalIndex, const config_t cfg,
<<<<<<< HEAD
    typename detector_t::detector_view_type det_data,
    measurement_collection_types::const_view measurements_view,
    vecmem::data::vector_view<const detray::geometry::barcode> barcodes_view,
    vecmem::data::vector_view<const unsigned int> upper_bounds_view,
=======
    typename detector_t::view_type det_data,
    measurement_collection_types::const_view measurements_view,
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
    bound_track_parameters_collection_types::const_view in_params_view,
    vecmem::data::vector_view<const unsigned int>
        n_measurements_prefix_sum_view,
    vecmem::data::vector_view<const unsigned int> ref_meas_idx_view,
    const unsigned int step, const unsigned int& n_max_candidates,
    bound_track_parameters_collection_types::view out_params_view,
    vecmem::data::vector_view<candidate_link> links_view,
    unsigned int& n_candidates);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/find_tracks.ipp"