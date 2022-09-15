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
 * @brief Output stream insertion of generic std::vector of Type
 * @tparam Type Element type of the vector
 * @param lhs Reference to the output stream to insert values into
 * @param rhs Reference to the vector to get the values from
 * @return Reference to the modified output stream
 */
template<class Type>
    std::ostream &operator<<(std::ostream &lhs, const std::vector<Type> &rhs) {
        // Iterate over all elements of the vector
        for (const Type &value: rhs) {
            // Insert element into output stream
            lhs << value << ' ';
        }
        // Return reference to modified stream
        return lhs;
    }

/**
 * @brief Generates usage string referring to the program name
 * @param program_name Name of the program (argument) to use.
 * @return Usage instructions as std::string
 */
std::string usage(const std::string &name) {
    return "Usage: " + name + " root filename";
}

// Preprocessing tool entrypoint
int main(int argc, char **argv) {
    // Two arguments necessary: Directory root and filename
    if (argc != 3) {
        // Write error/usage message to standard error stream
        std::cerr << argv[0] << ": Missing Arguments. " << usage(argv[0])
                  << std::endl;
        // Stop here with error code
        exit(EXIT_FAILURE);
    }
    // Image directory (database) and image filename (basename without suffix)
    std::string root = argv[1], filename = argv[2];

    // Load the image
    // auto[image, descriptor] = load_image(root, filename);
    auto image = std::get<0>(load_image(root, filename));
    
    // Image cropping to remove invalid regions after tilt transforms
    ENV(Angle<Radian>, crop_upper, 35.0_deg);
    ENV(Angle<Radian>, crop_lower, 00.0_deg);

    // Target to write the cropped image into
    decltype(image) cropped;
    // Short name for vertical resolution of panoramic image
    double vr = image.addOn.verticalResolution;
    // Preprocessing: crop parts of the image
    ns_simd::croppedView(
        image, (int) (crop_upper / vr), (int) (crop_lower / vr), cropped
    );

    // Construct warping bundle instance to use for the experiment
    auto warping = make_warping(cropped.w);

    // If configured verbose, print some more info
    if (verbose) {
        // Dump the environment configuration to standard output
        env_config::dump(std::cout, dump_default);
    }

    // Compute the ideal pixel scales (vector) using the scale plane stack
    // computation of the warping bundle
    std::vector<double> pixel = warping.spsComp->idealPixelScale(
        warping.spsComp->maxDenom(cropped)
    );
    // Get match type of constructed warping bundle instance
    using MATCH_TYPE = decltype(warping)::MATCH_TYPE;
    // Compute the post scale using the scale plane stack computation depending
    // on image size and type used by the matching step
    double post = warping.spsComp->idealPostScale(
        double(ns_simd::SIMDTypeInfo<MATCH_TYPE>::max()) / cropped.w
    );

    // Print the scaling parameters
    std::cout << pixel << " " << post << std::endl;

    // If no errors so far, computing scales is done - exit with success code
    return EXIT_SUCCESS;
}
