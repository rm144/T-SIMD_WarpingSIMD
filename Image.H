// ===========================================================================
// 
// Image.H --
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
    if (!ns_simd::loadPGM(image_path, image, panorama)) {
        // Loading the image failed (probably the file does not exist)
        throw ns_simd::SIMDException{
            __FUNCTION__, "Failed to load the image file: " + image_path
        };
    }
    // If specified, rotate the image
    if (rotate != 0.0_rad) {
        // Compute shift in pixels from angle and panorama image resolution
        int shift = (int) ((double) rotate / image.addOn.verticalResolution);
        // Image to write the rotated into
        ns_simd::SIMDImage<ns_simd::SIMDByte, SW, SA, ns_simd::Panorama> rot;
        // Rotate image horizontally by shift pixels
        ns_simd::rotateHor(image, shift, rot);
        // Add rotation to the image descriptor
        descriptor.delta_R += (double) rotate;
        // Extract roll-pitch parameters
        double delta_X = descriptor.delta_X, delta_Y = descriptor.delta_Y;
        // Adapt the tilt parameters to change of orientation:
        //  NOTE: This approximates the effect of the similarity transform
        //  describing the change of coordinate system due to rotating the
        //  view (reference) direction:
        //      Rz(-rotate) * Rx(delta_X) * Ry(delta_Y) * Rz(+rotate)
        //  Note: This essentially is a rotation of the tilt parameters
        //  represented as a vector (delta_X, delta_Y) in 2d:
        //      (delta_X', delta_Y') = Rz(-rotate) * (delta_X, delta_Y)
        descriptor.delta_X = cos(rotate) * delta_X + sin(rotate) * delta_Y;
        descriptor.delta_Y = cos(rotate) * delta_Y - sin(rotate) * delta_X;
        // Return the rotated image
        return std::pair{rot, descriptor};
    }
    // Return (normal) loaded image
    return std::pair{image, descriptor};
}

#endif /*EXP_IMAGE_H*/
