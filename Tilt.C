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
// Execution time measurement via struct timespec
#include "TimeMeasurement.H"

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
 * @param program_name Name of the program (argument) to use.
 * @return Usage instructions as std::string
 */
std::string usage(const std::string &name) {
    return "Usage: " + name
           + " root filename output_root solution [interpolation=NEAREST]";
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
    if (argc != 5 && argc != 6) {
        // Write error/usage message to standard error stream
        std::cerr << argv[0] << ": Missing Arguments. " << usage(argv[0])
                  << std::endl;
        // Stop here with error code
        exit(EXIT_FAILURE);
    }
    // Image directory (database) and image filename (basename without suffix)
    std::string root = argv[1], filename = argv[2], output_root = argv[3];
    // Tilt solution selector (number)
    tilt::Selector solution = tilt::Selector{std::stoi(argv[4])};
    // Optional interpolation mode selector (if not given defaults to nearest
    // neighbor interpolation)
    Interpolation interpolation = (argc == 6) ? Interpolation{
        std::stoi(argv[5])} : Interpolation::NEAREST;

    // Load the image
    auto[image, descriptor] = load_image(root, filename);

    // Image cropping to remove invalid regions after tilt transforms
    ENV(Angle<Radian>, crop_upper, 00.0_deg);
    ENV(Angle<Radian>, crop_lower, 00.0_deg);

    // Setup image cropping function used as preprocessing step during tilt
    // search (to remove lower/upper mask regions).
    auto crop = [crop_upper, crop_lower](auto &img) {
        // Image to fill with view of input image
        std::remove_reference_t<decltype(img)> cropped;
        // Short name for vertical resolution of panoramic image
        double vr = img.addOn.verticalResolution;
        // Preprocessing: crop parts of the image
        ns_simd::croppedView(
            img, (int) (crop_upper / vr), (int) (crop_lower / vr), cropped
        );
        // Return cropped image (view)
        return cropped;
    };

    // If configured verbose, print some more info
    if (verbose) {
        // Dump the environment configuration to standard output
        env_config::dump(std::cout, dump_default);
    }

    // Tilt correction solution selector
    auto correction = tilt::make_correction<decltype(image)>(solution);

    // Get the ground truth tilt parameters from the image descriptor
    auto tilt = tilt::RollPitch<Radian>{
        descriptor.delta_X, descriptor.delta_Y
    };

    // Start measuring time
    auto t0 = ns_simd::getTimeSpec();
    // Apply tilt correction to the image. Write to temporary as source and
    // destination of tilt operation may not overlap.
    // @formatter:off
    decltype(image) tmp; correction(image, tmp, tilt, interpolation);
    // @formatter:on
    // Stop measuring time
    auto t1 = ns_simd::getTimeSpec();
    // Accumulate time
    auto total_time = (long) ns_simd::timeSpecDiffUsec(t1, t0);

    // Build path to image file including directory and suffix
    std::string image_path = output_root + "/" + descriptor.filename + ".pgm";
    // Save the tilted image as configured
    if (!ns_simd::savePGM(image_path, crop(tmp))) {
        // Saving the image failed (probably the directory does not exist)
        throw ns_simd::SIMDException{
            __FUNCTION__, "Failed to save the image file: "
                          + image_path + ": " + strerror(errno)
        };
    }

    // Build path to image file including directory and suffix
    image_path = output_root + "/" + descriptor.filename + ".orig.pgm";
    // Move the original image to the output directory for comparison
    if (!ns_simd::savePGM(image_path, crop(image))) {
        // Saving the image failed (probably the directory does not exist)
        throw ns_simd::SIMDException{
            __FUNCTION__, "Failed to save the image file: "
                          + image_path + ": " + strerror(errno)
        };
    }

    // Print the image parameters and measured time
    std::cout << descriptor << ' ' << total_time << std::endl;

    // If no errors so far, experiment is done - exit with success code
    return EXIT_SUCCESS;
}
