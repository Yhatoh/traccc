/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/kokkos/seeding/seeding_algorithm.hpp"
#include "traccc/options/common_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"
#include "traccc/options/seeding_input_options.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>

// Kokkos include(s).
#include <Kokkos_Core.hpp>

namespace po = boost::program_options;

int seq_run(const traccc::seeding_input_options& /*i_cfg*/,
            const traccc::common_options& common_opts, bool run_cpu) {

    // Read the surface transforms
    auto surface_transforms =
        traccc::io::read_geometry(common_opts.detector_file);

    // Output stats
    uint64_t n_modules = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_kokkos = 0;

    // Configs
    traccc::seedfinder_config finder_config;
    traccc::spacepoint_grid_config grid_config(finder_config);
    traccc::seedfilter_config filter_config;

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    traccc::memory_resource mr{host_mr, &host_mr};

    traccc::seeding_algorithm sa(finder_config, grid_config, filter_config,
                                 host_mr);
    traccc::track_params_estimation tp(host_mr);

    vecmem::copy copy;

    // Kokkos Seeding algorithm
    traccc::kokkos::seeding_algorithm sa_kokkos{finder_config, grid_config,
                                                filter_config, mr, copy};

    // performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{});

    traccc::performance::timing_info elapsedTimes;

    // Loop over events
    for (unsigned int event = common_opts.skip;
         event < common_opts.events + common_opts.skip; ++event) {

        traccc::io::spacepoint_reader_output reader_output(&host_mr);
        traccc::seeding_algorithm::output_type seeds;
        traccc::track_params_estimation::output_type params;

        // Instantiate Kokkos containers/collections
        traccc::seed_collection_types::buffer seeds_kokkos_buffer(0,
                                                                  *(mr.host));

        {  // Start measuring wall time
            traccc::performance::timer wall_t("Wall time", elapsedTimes);

            /*-----------------
            hit file reading
            -----------------*/
            {
                traccc::performance::timer t("Hit reading  (cpu)",
                                             elapsedTimes);
                // Read the hits from the relevant event file
                traccc::io::read_spacepoints(
                    reader_output, event, common_opts.input_directory,
                    surface_transforms, common_opts.input_data_format);
            }  // stop measuring hit reading timer

            auto& spacepoints_per_event = reader_output.spacepoints;
            auto& modules_per_event = reader_output.modules;

            /*----------------------------
                Seeding algorithm
            ----------------------------*/

            traccc::spacepoint_collection_types::buffer
                spacepoints_kokkos_buffer(spacepoints_per_event.size(),
                                          mr.main);

            copy(vecmem::get_data(spacepoints_per_event),
                 spacepoints_kokkos_buffer);

            {
                traccc::performance::timer t("Seeding (kokkos)", elapsedTimes);

                // Reconstruct the spacepoints into seeds.
                seeds_kokkos_buffer = sa_kokkos(spacepoints_kokkos_buffer);
            }  // stop measuring seedin kokkos timer

            // CPU
            if (run_cpu) {
                traccc::performance::timer t("Seeding  (cpu)", elapsedTimes);
                seeds = sa(spacepoints_per_event);
            }  // stop measuring seeding cpu timer

            /*----------------------------
            Track params estimation
            ----------------------------*/

            // CPU

            if (run_cpu) {
                traccc::performance::timer t("Track params  (cpu)",
                                             elapsedTimes);
                params =
                    tp(std::move(spacepoints_per_event), seeds,
                       modules_per_event, {0.f, 0.f, finder_config.bFieldInZ});
            }  // stop measuring track params cpu timer

        }  // Stop measuring wall time
    }

    if (common_opts.check_performance) {
        sd_performance_writer.finalize();
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_spacepoints << " spacepoints from "
              << n_modules << " modules" << std::endl;
    std::cout << "- created (cpu)  " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (kokkos) " << n_seeds_kokkos << " seeds"
              << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    // Initialise both Kokkos and GoogleTest.
    Kokkos::initialize(argc, argv);

    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::common_options common_opts(desc);
    traccc::seeding_input_options seeding_input_cfg(desc);
    desc.add_options()("run_cpu", po::value<bool>()->default_value(false),
                       "run cpu tracking as well");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    seeding_input_cfg.read(vm);
    auto run_cpu = vm["run_cpu"].as<bool>();

    std::cout << "Running " << argv[0] << " " << common_opts.detector_file
              << " " << common_opts.input_directory << " " << common_opts.events
              << std::endl;

    int ret = seq_run(seeding_input_cfg, common_opts, run_cpu);

    // Finalise Kokkos.
    Kokkos::finalize();

    return ret;
}
