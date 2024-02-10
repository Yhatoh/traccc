/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/seeding/experimental/spacepoint_formation.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/edm/spacepoint.hpp"

// Detray include(s).
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/intersection/detail/trajectories.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

TEST(spacepoint_formation, cuda) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;
    traccc::memory_resource mr{mng_mr};

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    // Use rectangle surfaces
    detray::mask<detray::rectangle2D<>> rectangle{
        0u, 10000.f * detray::unit<scalar>::mm,
        10000.f * detray::unit<scalar>::mm};

    // Plane alignment direction (aligned to x-axis)
    detray::detail::ray<transform3> traj{{0, 0, 0}, 0, {1, 0, 0}, -1};

    // Position of planes (in mm unit)
    std::vector<scalar> plane_positions = {20.f,  40.f,  60.f,  80.f, 100.f,
                                           120.f, 140.f, 160.f, 180.f};

    detray::tel_det_config<> tel_cfg{rectangle};
    tel_cfg.positions(plane_positions);
    tel_cfg.pilot_track(traj);

    // Create telescope geometry
    auto [det, name_map] = create_telescope_detector(mng_mr, tel_cfg);
    using device_detector_type =
        detray::detector<detray::telescope_metadata<detray::rectangle2D<>>,
<<<<<<< HEAD
                         covfie::field_view, detray::device_container_types>;

    // Surface lookup
    auto surfaces = det.surface_lookup();
=======
                         detray::device_container_types>;

    // Surface lookup
    auto surfaces = det.surfaces();
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846

    // Prepare measurement collection
    measurement_collection_types::host measurements{&mng_mr};

    // Add a measurement at the first plane
    measurements.push_back({{7.f, 2.f}, {0.f, 0.f}, surfaces[0].barcode()});

    // Add a measurement at the last plane
    measurements.push_back({{10.f, 15.f}, {0.f, 0.f}, surfaces[8u].barcode()});

    // Run spacepoint formation
<<<<<<< HEAD
    cuda::experimental::spacepoint_formation<device_detector_type> sp_formation(
        mr, copy, stream);
=======
    traccc::cuda::experimental::spacepoint_formation<device_detector_type>
        sp_formation(mr, copy, stream);
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
    auto spacepoints_buffer =
        sp_formation(detray::get_data(det), vecmem::get_data(measurements));

    spacepoint_collection_types::device spacepoints(spacepoints_buffer);

    // Check the results
    EXPECT_EQ(copy.get_size(spacepoints_buffer), 2u);
    std::set<point3> test;
    test.insert(spacepoints[0].global);
    test.insert(spacepoints[1].global);

    std::set<point3> ref;
    ref.insert({180.f, 10.f, 15.f});
    ref.insert({20.f, 7.f, 2.f});

    EXPECT_EQ(test, ref);
<<<<<<< HEAD
}
=======
}
>>>>>>> f2918520ddb7c6e26d80f74d95a69f87a90be846
