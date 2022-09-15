#ifndef EXP_IMAGE_H
#define EXP_IMAGE_H

// Math functions for tilt parameter adaptation
#include <cmath>
// std::pair
#include <utility>
// (Image) Filenames are std::strings
#include <string>
// Input/Output streams for loading parameters
#include <iostream>
// C++ file stream input/output
#include <fstream>

// SIMD library Exception class
#include "SIMDException.H"
// SIMD Image type and Panorama addon
#include "SIMDImage.H"
// Functions operating on SIMDImage (e.g. loadPGM, rotateHor)
#include "SIMDImageFunctions.H"

// Strong type angular units
#include "Angle.H"

// Tilt parameter models
#include "TiltModel.H"

/**
 * @brief Image descriptor parameters
 */
struct ImageParams {
    // Image filename (without path/prefix/suffix)
    std::string filename{};
    // Grid position and tilt/direction angles
    double pos_X = 0, pos_Y = 0, delta_X = 0, delta_Y = 0, delta_R = 0;

    /**
     * @brief Extracts parameters from input stream
     */
    friend std::istream &operator>>(std::istream &lhs, ImageParams &rhs) {
        // Extract individual parameters from the input stream
        lhs >> rhs.filename >> rhs.pos_X >> rhs.pos_Y
            >> rhs.delta_X >> rhs.delta_Y >> rhs.delta_R;
        // Return reference to the modified stream (left-hand side).
        return lhs;
    }

    /**
     * @brief Inserts parameters into output stream
     */
    friend std::ostream &operator<<(std::ostream &lhs, ImageParams &rhs) {
        // Insert individual parameters into the output stream
        lhs << rhs.filename << " " << rhs.pos_X << " " << rhs.pos_Y << " "
            << rhs.delta_X << " " << rhs.delta_Y << " " << rhs.delta_R;
        // Return reference to the modified stream (left-hand side).
        return lhs;
    }
};

/**
 * @brief Extracts Panorama addon parameters from input stream
 * @param lhs Input stream to extract parameters from
 * @param rhs Panorama addon to fill
 * @return Reference to input stream after extraction (lhs)
 */
std::istream &operator>>(std::istream &lhs, ns_simd::Panorama &rhs) {
    // Extract parameter in order (extraction returns ref)
    return lhs >> rhs.verticalResolution >> rhs.horizon;
}

/**
 * @brief Loads a SIMDImage panorama given the root directory and filename.
 * @param root Root directory of the image database
 * @param filename Filename (basename without .txt or .pgm suffix) of the image
 * to load
 * @return Loaded panoramic SIMDImage
 */
auto load_image(const std::string &root, const std::string &filename,
                const Angle<Radian> rotate = 0.0_rad) {
    // Open the image descriptor file (containing support data, e.g. panorama
    // parameters)
    std::ifstream descriptor_file(root + "/" + filename + ".txt");
    // Image parameter descriptor to read file into
    ImageParams descriptor;
    // (Uninitialized) Panorama addon to fill
    ns_simd::Panorama panorama;
    // Extract parameters from descriptor file
    descriptor_file >> descriptor >> panorama;
    // Uninitialized SIMDImage to read image file into
    ns_simd::SIMDImage<ns_simd::SIMDByte, SW, SA, ns_simd::Panorama> image;
    // Build the full path to the image file including directory and suffix
    std::string image_path = root + "/" + filename + ".pgm";
    // Load the image from file and insert the panorama addon
    if(!ns_simd::loadPGM(image_path, image, panorama)) {
        // Loading the image failed (probably the file does not exist)
        throw ns_simd::SIMDException{
            __FUNCTION__, "Failed to load the image file: " + image_path
        };
    }
    // If specified, rotate the image
    if (rotate != 0.0_rad) {
        // Compute shift in pixels from angle and panorama image resolution
        int shift = (double) (rotate / image.addOn.verticalResolution);
        // Image to write the rotated into
        ns_simd::SIMDImage<ns_simd::SIMDByte, SW, SA, ns_simd::Panorama> rot;
        // Rotate image horizontally by shift pixels
        ns_simd::rotateHor(image, shift, rot);
        // Add rotation to the image descriptor
        descriptor.delta_R += (double)rotate;
        // Convert tilt parameters to axis-angle representation
        //  NOTE: It is easier to adapt tilt to orientation change in axis angle
        //  representation. Axis-angle decouples magnitude and direction of
        //  tilt with only the latter changing.
        tilt::AxisAngle<Radian> tilt = tilt::RollPitch<Radian> {
            descriptor.delta_X, descriptor.delta_Y
        };
        // Add changed rotation to tilt in axis-angle representation
        //  NOTE: Adapt tilt in opposite direction of rotation
        tilt.theta_R -= rotate;
        // Extract modified roll-pitch parameters
        auto [delta_X, delta_Y] = tilt::RollPitch<Radian>{tilt};
        // Adapt the tilt parameters to change of orientation
        descriptor.delta_X = (double)delta_X;
        descriptor.delta_Y = (double)delta_Y;
        // Return the rotated image
        return std::pair{rot, descriptor};
    }
    // Return (normal) loaded image
    return std::pair{image, descriptor};
}

#endif /*EXP_IMAGE_H*/
