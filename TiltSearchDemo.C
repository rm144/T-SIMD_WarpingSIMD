// ===========================================================================
// 
// TiltSearchDemo.C --
// Christoph Berganski
// 
// This source code file is part of the following software:
// 
//    - the low-level C++ template SIMD library
//    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
//      for local visual homing.
// 
// The software is provided based on the accompanying license agreement
// in the file LICENSE or LICENSE.doc. The software is provided "as is"
// without any warranty by the licensor and without any liability of the
// licensor, and the software may not be distributed by the licensee; see
// the license agreement for details.
// 
// (C) Ralf Möller
//     Computer Engineering
//     Faculty of Technology
//     Bielefeld University
//     www.ti.uni-bielefeld.de
// 
// ===========================================================================

// C standard library header for C++ (contains EXIT_SUCCESS)
#include <cstdlib>
// C++ stream input/output
#include <iostream>

// Strong type angular units
#include "Angle.H"

// SIMD Image type
#include "SIMDImage.H"
// Functions operating on SIMDImages
#include "SIMDImageFunctions.H"
// Bundle of entire warping/min-warping functionality
#include "WarpingBundle.H"

// Grid search strategy
#include "SearchGrid.H"
// Pattern search strategy
#include "SearchPattern.H"
// Nelder-Mead simplex search strategy
#include "SearchNelderMead.H"
// Add logging to the generic search
#include "SearchLogger.H"

// Tilt correction solutions and image transformations
#include "TiltCorrection.H"
// Warping tilt search frontend
#include "TiltSearch.H"

// Environment variable configuration options
#include "Env.H"
// Experiment setup (select types, widths, warping bundle, ... to use)
#include "Setup.H"
// Experiment image data loading
#include "Image.H"

// Run program in verbose mode, configured by environment
ENV_BOOL(verbose, false);
// Dump default options if verbose?
ENV_BOOL(dump_default, true);

// Search strategy selectors
enum Strategy : int {
    GRID = 0, PATTERN = 1, NELDER_MEAD = 2
};

/**
 * @brief Composes a warping tilt search wrapping a warping bundle at runtime
 * based on configuration from YAML::Node
 * @tparam Warping Type of WarpingBundle to wrap
 * @tparam Log Type of search logging to use (either std::list<std::string> or
 * search::NoLog)
 * @param config YAML::Node holding the configuration
 * @param warping Reference to warping bundle to wrap
 * @param log Reference to search log to fill (or search::NoLog{} to disable
 * logging)
 * @return The composed tilt::Search instance
 */
template<class Warping, class Log = search::NoLog>
    auto make_search(Warping &warping, Log &&log = Log{}) {
        // Environment option to select the search strategy
        DERIVE_ENV(Strategy, int, search_strategy, PATTERN);
        // Environment option to select the tilt correction solution
        ENV(tilt::Selector, tilt_solution, tilt::EXACT);

        // Hardcoded search range
        auto X = search::Range{-0.14_rad, +0.14_rad, +0.02_rad};
        auto Y = search::Range{-0.14_rad, +0.14_rad, +0.02_rad};

        // Construct box-constraint function given the ranges
        //  NOTE: This is used by Pattern and Nelder-Mead search to constrain
        //  those to the same range as the grid search
        auto limit = [X, Y](auto x, auto y) -> bool {
            return x >= X.min && x <= X.max && y >= Y.min && y <= Y.max;
        };

        // Construct tilt search based on strategy
        switch (search_strategy) {
            // Construct grid search
            case GRID: {
                // Setup search grid from config
                search::Grid searcher{X, Y};
                // Compose the warping tilt search
                return tilt::Search{
                    tilt::RollPitch<Radian>{},
                    tilt_solution,
                    search::maybe_log(searcher, log),
                    warping
                };
            }
                // Construct pattern search
            case PATTERN: {
                // Hardcoded initial search pattern
                auto searcher = search::Pattern{
                    std::tuple{0.0_rad, 0.0_rad}, 0.14, 0.02, 0.5
                };
                // Compose the warping tilt search
                return tilt::Search{
                    tilt::RollPitch<Radian>{},
                    tilt_solution,
                    search::maybe_log(searcher, log),
                    warping,
                    limit
                };
            }
                // Construct Nelder-Mead searcher
            case NELDER_MEAD: {
                // Hardcoded initial search simplex
                auto searcher = search::NelderMead{
                    std::array{
                        std::tuple{-0.07_rad, -0.07_rad},
                        std::tuple{+0.14_rad, +0.00_rad},
                        std::tuple{+0.00_rad, +0.14_rad}
                    },
                    1.0, 2.0, 0.5, 0.5, 0.04, 50
                };
                // Compose the warping tilt search
                return tilt::Search{
                    tilt::RollPitch<Radian>{},
                    tilt_solution,
                    search::maybe_log(searcher, log),
                    warping,
                    limit
                };
            }
                // Unsupported / Not implemented search strategy
            default: {
                // Adhoc exception class
                class Exception : public std::exception {
                    // Returns hardcoded message
                    const char *what() const noexcept override { // NOLINT
                        return "Search Strategy Not Implemented";
                    }
                };
                // Throw the adhoc exception
                throw Exception();
            }
        }

    }

/**
 * @brief Inserts argument pack into output stream adding space between args.
 * @tparam Args Argument pack to insert
 * @param stream Output stream to insert arguments into
 * @param args Arguments to insert into the stream
 */
template<class... Args>
    void print_args(std::ostream &stream, Args &&... args) {
        // Fold arguments into the stream, each followed by space
        ((stream << args << ' '), ...);
    }

/**
* @brief Generates usage string referring to the program name
* @param name Name of the program (argument) to use.
* @return Usage instructions as std::string
*/
std::string usage(const std::string &name) {
    return "Usage: " + name + " ss_root ss cv_root cv [rotate-ss] [rotate-cv]";
}

// Experiment tool entrypoint
int main(int argc, char **argv) {
    // Two arguments necessary: Directory root and filename
    if (argc != 5 && argc != 6 && argc != 7) {
        // Write error/usage message to standard error stream
        std::cerr << argv[0] << ": Missing Arguments. " << usage(argv[0])
                  << std::endl;
        // Stop here with error code
        exit(EXIT_FAILURE);
    }
    // Image directory (database) and image filename (basename without suffix)
    std::string ss_root = argv[1], ss_filename = argv[2]; // Snapshot image
    std::string cv_root = argv[3], cv_filename = argv[4]; // Current view image

    // Optional last arguments specify to rotate ss and cv images
    Angle<Radian> rotate_ss{(argc >= 6) ? std::stod(argv[5]) : 0.0};
    Angle<Radian> rotate_cv{(argc >= 7) ? std::stod(argv[6]) : 0.0};

    // Load the snapshot and current view images
    auto[ss, ss_descriptor] = load_image(ss_root, ss_filename, rotate_ss);
    auto[cv, cv_descriptor] = load_image(cv_root, cv_filename, rotate_cv);

    // Warping requires the ss and cv images to have the same size
    assert(ss.w == cv.w && ss.h == cv.h);

    // Construct warping bundle instance to use for the experiment
    auto warping = make_warping(ss.w);
    // Construct parameters for warping run
    auto params = make_warping_params();

    // Image cropping to remove invalid regions after tilt transforms
    ENV(Angle<Radian>, crop_upper, 35.0_deg);
    ENV(Angle<Radian>, crop_lower, 00.0_deg);

    // Setup image cropping function used as preprocessing step during tilt
    // search (to remove lower/upper mask regions).
    auto crop = [crop_upper, crop_lower](auto &image) {
        // Image to fill with view of input image
        std::remove_reference_t<decltype(image)> cropped;
        // Short name for vertical resolution of panoramic image
        double vr = image.addOn.verticalResolution;
        // Preprocessing: crop parts of the image
        ns_simd::croppedView(
            image, (int) (crop_upper / vr), (int) (crop_lower / vr), cropped
        );
        // Return cropped image (view)
        return cropped;
    };

    // Interpolation method to use
    DERIVE_ENV(Interpolation, int, interpolation, Interpolation::NEAREST);

    // Environment variable for enabling logging of the search process
    ENV_BOOL(logging, false);
    // Logging of search (states) string representation
    std::list<std::string> log;
    // Set up the tilt search wrapping the warping instance
    //  NOTE: Directly call lambda expression without arguments to wrap
    //  conditional construction.
    tilt::Search tilt_search = [&] {
        // Make search with or without logging
        return logging ? make_search(warping, log)
                       : make_search(warping);
    }();

    // If configured verbose, print some more info
    if (verbose) {
        // Dump the environment configuration to standard output
        env_config::dump(std::cout, dump_default);
    }

    // Start measuring time
    auto t0 = ns_simd::getTimeSpec();
    // Run the (constrained) tilt search
    auto[min, argmin] = tilt_search.run(ss, cv, params, crop, interpolation);
    // Stop measuring time
    auto t1 = ns_simd::getTimeSpec();
    // Accumulate time
    auto total_time = (long) ns_simd::timeSpecDiffUsec(t1, t0);

    // Print snapshot and current view image support data
    std::cout << ss_descriptor << ' ' << cv_descriptor << ' ';
    // Print search results (warping, tilt parameters and time measurements)
    print_args(
        std::cout, argmin, min, warping.count(), warping.avg_time(), total_time
    );
    // Add newline and flush the stream
    std::cout << std::endl;

    // Print the log of collected search states (might be empty in case of
    // logging == false).
    for (const auto &state: log) {
        // Prefix log strings with '>' to allow for filtering of outputs
        std::cout << "> " << state << std::endl;
    }

    // Compute true home directions (beta angle)
    double beta_true = atan2(ss_descriptor.pos_Y - cv_descriptor.pos_Y,
                             ss_descriptor.pos_X - cv_descriptor.pos_X);
    // Compute true orientation angle
    double psi_true = cv_descriptor.delta_R - ss_descriptor.delta_R;
    // Compute true alpha parameter
    double alpha_true = psi_true - beta_true - cv_descriptor.delta_R + M_PI;
    // Compute estimated beta
    double beta_est = min.psi - cv_descriptor.delta_R - min.alpha + M_PI;

    // Compute error of home direction estimates
    double beta_err = acos(cos(beta_true - beta_est));
    double alpha_err = acos(cos(alpha_true - min.alpha));
    double psi_err = acos(cos(psi_true - min.psi));

    // Extract true and estimated tilt angles
    double dx_true = cv_descriptor.delta_X, dx_est = (double) argmin.delta_X;
    double dy_true = cv_descriptor.delta_Y, dy_est = (double) argmin.delta_Y;

    // Compute tilt estimate error
    double tilt_error = acos(
        cos(dx_true - dx_est) * cos(dy_true) * cos(dy_est) + sin(
            dy_true) * sin(dy_est));

    // Print  home direction angle estimates
    printf(
        "# EST:   beta=%f, alpha=%f, psi=%f\n", beta_est, min.alpha, min.psi
    );
    // Print true (expected) home direction angles
    printf(
        "# TRUE:  beta=%f, alpha=%f, psi=%f\n", beta_true, alpha_true, psi_true
    );
    // Print home direction errors
    printf(
        "# ERROR: beta=%f, alpha=%f, psi=%f\n", beta_err, alpha_err, psi_err
    );
    // Print home direction angle estimates
    printf(
        "# TILT EST:   delta_X=%f, delta_Y=%f\n", dx_est, dy_est
    );
    // Print true (expected) home direction angles
    printf(
        "# TILT TRUE:  delta_X=%f, delta_Y=%f\n", dx_true, dy_true
    );
    // Print home direction errors
    printf(
        "# TILT ERROR: %f\n", tilt_error
    );
    // Print time measurements
    printf(
        "# TIME:  warping=%f, total=%ld\n", warping.avg_time(), total_time
    );
    // Print count of warping runs
    printf(
        "# COUNT: warping=%zu\n", warping.count()
    );

    // If no errors so far, experiment is done - exit with success code
    return EXIT_SUCCESS;
}
