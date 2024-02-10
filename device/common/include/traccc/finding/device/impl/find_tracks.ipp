/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"

// System include(s).
#include <limits>

namespace traccc::device {

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
    unsigned int& n_candidates) {

    // Detector
    detector_t det(det_data);

    // Measurement
    measurement_collection_types::const_device measurements(measurements_view);
<<<<<<< HEAD
    vecmem::device_vector<const detray::geometry::barcode> barcodes(
        barcodes_view);
    vecmem::device_vector<const unsigned int> upper_bounds(upper_bounds_view);
=======
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846

    // Input parameters
    bound_track_parameters_collection_types::const_device in_params(
        in_params_view);

    // Output parameters
    bound_track_parameters_collection_types::device out_params(out_params_view);

    // Links
    vecmem::device_vector<candidate_link> links(links_view);

    // Prefix sum of the number of measurements per parameter
    vecmem::device_vector<const unsigned int> n_measurements_prefix_sum(
        n_measurements_prefix_sum_view);

<<<<<<< HEAD
    // Search for out_param index
    const auto lo1 = thrust::lower_bound(thrust::seq, n_threads.begin(),
                                         n_threads.end(), globalIndex + 1);
    const auto in_param_id = std::distance(n_threads.begin(), lo1);

    // Get measurements on surface
    unsigned int ref;
    if (lo1 == n_threads.begin()) {
        ref = 0;
    } else {
        ref = *(lo1 - 1);
    }

    // Get barcode
    const auto bcd = in_params.at(in_param_id).surface_link();

    const auto lo2 =
        thrust::lower_bound(thrust::seq, barcodes.begin(), barcodes.end(), bcd);
    const auto bcd_id = std::distance(barcodes.begin(), lo2);

    // Find the range for measurements
    std::pair<unsigned int, unsigned int> range;
    if (lo2 == barcodes.begin()) {
        range.first = 0u;
        range.second = upper_bounds[bcd_id];
    } else {
        range.first = upper_bounds[bcd_id - 1];
        range.second = upper_bounds[bcd_id];
    }

    const unsigned int offset = globalIndex - ref;
    const unsigned int stride = offset * n_measurements_per_thread;
    const unsigned int n_meas_on_surface = range.second - range.first;

    // Iterate over the measurements
    const detray::surface<detector_t> sf{det, bcd};
=======
    // Reference (first) measurement index per parameter
    vecmem::device_vector<const unsigned int> ref_meas_idx(ref_meas_idx_view);
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846

    // Last step ID
    const unsigned int previous_step =
        (step == 0) ? std::numeric_limits<unsigned int>::max() : step - 1;

    const unsigned int n_measurements_sum = n_measurements_prefix_sum.back();
    const unsigned int stride = globalIndex * cfg.n_measurements_per_thread;

    vecmem::device_vector<const unsigned int>::iterator lo1;

    for (unsigned int i_meas = 0; i_meas < cfg.n_measurements_per_thread;
         i_meas++) {
        const unsigned int idx = stride + i_meas;

        if (idx >= n_measurements_sum) {
            break;
        }
        const unsigned int meas_idx = i + stride + range.first;

        if (i_meas == 0 || idx == *lo1) {
            lo1 = thrust::lower_bound(thrust::seq,
                                      n_measurements_prefix_sum.begin(),
                                      n_measurements_prefix_sum.end(), idx + 1);
        }

        const unsigned int in_param_id =
            std::distance(n_measurements_prefix_sum.begin(), lo1);
        const detray::geometry::barcode bcd =
            in_params.at(in_param_id).surface_link();
        const unsigned int offset =
            lo1 == n_measurements_prefix_sum.begin() ? idx : idx - *(lo1 - 1);
        const unsigned int meas_idx = ref_meas_idx.at(in_param_id) + offset;
        bound_track_parameters in_par = in_params.at(in_param_id);
<<<<<<< HEAD
        const auto meas = measurements.at(meas_idx);
=======

        const auto& meas = measurements.at(meas_idx);
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
        track_state<typename detector_t::transform3> trk_state(meas);
        const detray::surface<detector_t> sf{det, bcd};

        // Run the Kalman update
        sf.template visit_mask<
            gain_matrix_updater<typename detector_t::transform3>>(trk_state,
                                                                  in_par);
        // Get the chi-square
        const auto chi2 = trk_state.filtered_chi2();

        if (chi2 < cfg.chi2_max) {

            // Add measurement candidates to link
            vecmem::device_atomic_ref<unsigned int> num_candidates(
                n_candidates);

            const unsigned int l_pos = num_candidates.fetch_add(1);

<<<<<<< HEAD
            // @TODO; Consider max_num_branches_per_surface
=======
            if (l_pos >= n_max_candidates) {
                n_candidates = n_max_candidates;
                return;
            }

>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
            links[l_pos] = {{previous_step, in_param_id}, meas_idx};

            out_params[l_pos] = trk_state.filtered();
        }
    }
}

}  // namespace traccc::device