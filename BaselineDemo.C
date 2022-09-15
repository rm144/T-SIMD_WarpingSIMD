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

// Tilt correction solutions and image transformations
#include "TiltCorrection.H"

// Environment variable configuration options
#include "Env.H"
// Experiment setup (select types, widths, warping bundle, ... to use)
#include "Setup.H"
// Experiment image data loading
#include "Image.H"

// Run program in verbose mode, configured by environment
BOOL(verbose, false);
// Dump default options if verbose?
BOOL(dump_default, true);

/**
 * @brief Generates usage string referring to the program name
 * @param name Name of the program (argument) to use.
 * @return Usage instructions as std::string
 */
std::string usage(const std::string &name) {
    return "Usage: " + name +
           " ss_root ss cv_root cv [rotate-ss] [rotate-cv] [no-tilt]";
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

// Experiment tool entrypoint
int main(int argc, char **argv) {
    // Two arguments necessary: Directory root and filename
    if (argc != 5 && argc != 6 && argc != 7 && argc != 8) {
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
    // Optional argument specifies to do no tilt correction
    bool no_tilt = (argc >= 8) && std::string(argv[7]) == "no-tilt";

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

    // Environment option to select the tilt correction solution
    ENV(tilt::Selector, tilt_solution, tilt::EXACT);
    // Tilt correction solution selector
    auto correction = tilt::make_correction<decltype(cv)>(tilt_solution);
    // Interpolation method to use
    DERIVE_ENV(Interpolation, int, interpolation, Interpolation::NEAREST);

    // Warping run results output variable
    WarpingBundle::Result min;
    // Tilt correction parameters result
    tilt::RollPitch<Radian> argmin{0.0_rad, 0.0_rad};

    // If configured verbose, print some more info
    if (verbose) {
        // Dump the environment configuration to standard output
        env_config::dump(std::cout, dump_default);
    }

    // Start measuring time
    auto t0 = ns_simd::getTimeSpec();
    // Select which baseline test to run: tilt correction or not tilt correction
    if (no_tilt) {
        // Temporary for tilt corrected and preprocessed images
        auto pre_ss = crop(ss), pre_cv = crop(cv);
        // Run warping on cropped images without tilt correction
        warping.run(pre_ss, pre_cv, params, min);
    } else {
        // Get the ground truth tilt parameters from current view descriptor
        argmin = tilt::RollPitch<Radian>{
            cv_descriptor.delta_X, cv_descriptor.delta_Y
        };
        // Temporary for tilt corrected and preprocessed images
        decltype(cv) tmp;
        // Apply tilt correction to current view image
        correction(cv, tmp, argmin, interpolation);
        // Crop invalid regions from images
        auto pre_ss = crop(ss), pre_cv = crop(tmp);
        // Run warping on cropped and tilt corrected images
        warping.run(pre_ss, pre_cv, params, min);
    }
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

    // Compute tilt magnitude (error does not make sense here - it is either
    // zero if corrected or the magnitude if not)
    double tilt_magnitude = acos(cos(dx_true) * cos(dy_true));

    // Print home direction angle estimates
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
        "# TILT MAGNITUDE: %f\n", tilt_magnitude
    );
    // Print time measurements
    printf(
        "# TIME:  warping=%f, total=%ld\n", warping.avg_time(), total_time
    );
    // Print count of warping runs
    printf(
        "# COUNT: warping=%ld\n", warping.count()
    );

    // If no errors so far, experiment is done - exit with success code
    return EXIT_SUCCESS;
}
