#ifndef OCAM_MODEL_H
#define OCAM_MODEL_H

// Math functions: root, powers, trig. to compute coordinate mappings
#include <cmath>
// Dynamic size generic vector (array). Used to hold polynomial coefficients
#include <vector>
// Standard library string. Used to specify filename
#include <string>
// Standard library exception base class (to derive custom exception)
#include <exception>

// Transformation matrices
#include "TransformMatrix.H"

// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"
// SIMD Image type
#include "SIMDImage.H"
// Functions operating on SIMDImages
#include "SIMDImageFunctions.H"


/**
 * @brief Camera calibration model interface (loading+transforming)
 * @see OCamCalib-Toolbox by Davide Scaramuzza
 */
class OCamModel {
    // Coefficients of the camera to world mapping polynomial
    std::vector<double> cam2world_coeff;
    // Coefficients of the world to camera mapping polynomial
    std::vector<double> world2cam_coeff;
    // Center coordinates
    double center_x, center_y;
    // OCamCalib (c,d,e) parameters (?)
    double c, d, e;
    // Image width and height
    double width, height;

public:
    /**
     * @brief Loads camera calibration from text file.
     * @see OCamCalib-Toolbox by Davide Scaramuzza
     * @param filename Path to calibration file (e.g. calib_results.txt)
     * @param verbose Enables verbose loading (prints calibration file comments)
     */
    explicit OCamModel(const char *filename, const bool verbose = false) {
        // Adhoc exception class
        class Exception : public std::exception {
            // Returns hardcoded message
            const char *what() const noexcept override {
                return "Error while reading the calibration file";
            }
        };

        // TODO: This "parsing" should really be refactored...

        // Open the calibration file for reading
        FILE *calib = fopen(filename, "r");

        if (!calib)
            // Throw the adhoc exception
            throw Exception();

        // Read comments and empty lines into this buffer
        char buffer[1024];

        // Drop first line (comment)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        if (!fgets(buffer, sizeof(buffer), calib))
            // Throw the adhoc exception
            throw Exception();
        if (verbose) printf("%s", buffer);

        // Drop next line (empty)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        // Read number of polynomial coefficients for the DIRECT mapping
        // function and resize the coefficients vector
        std::size_t n_coeff;

        if (fscanf(calib, "%lu ", &n_coeff) != 1)
            // Throw the adhoc exception
            throw Exception();

        cam2world_coeff.resize(n_coeff);

        for (std::size_t idx = 0; idx < n_coeff; ++idx)
            if (fscanf(calib, "%lf ", &cam2world_coeff.at(idx)) != 1)
                // Throw the adhoc exception
                throw Exception();

        // Drop next line (comment)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        if (!fgets(buffer, sizeof(buffer), calib))
            // Throw the adhoc exception
            throw Exception();
        if (verbose) printf("%s", buffer);

        // Drop next line (empty)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        // Read number of polynomial coefficients for the inverse mapping
        // function and resize the coefficients vector
        if (fscanf(calib, "%lu ", &n_coeff) != 1)
            // Throw the adhoc exception
            throw Exception();

        world2cam_coeff.resize(n_coeff);

        for (std::size_t idx = 0; idx < n_coeff; ++idx)
            if (fscanf(calib, "%lf ", &world2cam_coeff.at(idx)) != 1)
                // Throw the adhoc exception
                throw Exception();

        // Drop next line (comment)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        if (!fgets(buffer, sizeof(buffer), calib))
            // Throw the adhoc exception
            throw Exception();
        if (verbose) printf("%s", buffer);

        // Drop next line (empty)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        // Read image center coordinates
        if (fscanf(calib, "%lf %lf\n", &center_x, &center_y) != 2)
            // Throw the adhoc exception
            throw Exception();

        // Drop next line (comment)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        if (!fgets(buffer, sizeof(buffer), calib))
            // Throw the adhoc exception
            throw Exception();
        if (verbose) printf("%s", buffer);

        // Drop next line (empty)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        // Read affine parameters "c", "d", "e"
        if (fscanf(calib, "%lf %lf %lf\n", &c, &d, &e) != 3)
            // Throw the adhoc exception
            throw Exception();

        // Drop next line (comment)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        if (!fgets(buffer, sizeof(buffer), calib))
            // Throw the adhoc exception
            throw Exception();
        if (verbose) printf("%s", buffer);

        // Drop next line (empty)
        if (fscanf(calib, "\n") != 0)
            // Throw the adhoc exception
            throw Exception();

        // Read image size: "height" and "width"
        if (fscanf(calib, "%lf %lf\n", &width, &height) != 2)
            // Throw the adhoc exception
            throw Exception();

        if (fclose(calib))
            // Throw the adhoc exception
            throw Exception();
    }

    /**
     * @brief Loads camera calibration from text file.
     * @see OCamCalib-Toolbox by Davide Scaramuzza
     * @param filename Path to calibration file (e.g. calib_results.txt)
     * @param verbose Enables verbose loading (prints calibration file comments)
     */
    explicit OCamModel(const std::string &filename, const bool verbose = false)
        : OCamModel(filename.c_str(), verbose) {}

    /**
     * @brief World to camera coordinate mapping (inverse mapping function)
     * @see OCamCalib-Toolbox by Davide Scaramuzza
     * @param W 3d world coordinate vector
     * @param C 2d image coordinate vector (output)
     */
    void world2cam(const double W[3], double C[2]) const {
        // Length of the world  coordinate vector
        const double norm = std::sqrt(W[0] * W[0] + W[1] * W[1]);
        // Non-zero length vector needs proper mapping
        if (norm != 0.0) { // hm...
            // Kind of polar coordinate angle of the world coordinates
            const double theta = std::atan(W[2] / norm);
            // Inverse of the norm (length) for normalizing (scaling)
            const double invnorm = 1.0 / norm;
            // Start evaluating the polynomial with the first coefficient
            double rho = world2cam_coeff.at(0);
            // Evaluate the polynomial (starting with the highest power)
            for (std::size_t i = world2cam_coeff.size() - 1; i > 0; --i)
                rho += std::pow(theta, i) * world2cam_coeff.at(i);
            // Coordinate mapping
            const double x = W[0] * invnorm * rho;
            const double y = W[1] * invnorm * rho;
            C[0] = x * c + y * d + center_x;
            C[1] = x * e + y + center_y;
        } else {
            // Zero length vectors are just projected to the center
            C[0] = center_x;
            C[1] = center_y;
        }
    }

    /**
     * @brief Camera to world coordinate mapping (forward mapping function)
     * @see OCamCalib-Toolbox by Davide Scaramuzza
     * @param C 2d image coordinate vector
     * @param W 3d world coordinate vector (output)
     */
    void cam2world(const double C[2], double W[3]) const {
        // 1/det(A), where A = [c,d;e,1] as in the Matlab file
        const double invdet = 1 / (c - d * e);
        // ?
        double xp = invdet * ((C[0] - center_x) - d * (C[1] - center_y));
        double yp = invdet * (-e * (C[0] - center_x) + c * (C[1] - center_y));
        // distance [pixels] of  the point from the image center
        const double r = std::sqrt(xp * xp + yp * yp);
        // Start evaluating the polynomial with the first coefficient
        double zp = cam2world_coeff.at(0);
        // Evaluate the polynomial (starting with the highest power)
        for (std::size_t i = cam2world_coeff.size(); i > 0; --i)
            zp += std::pow(r, i) * cam2world_coeff.at(i);
        // Normalize
        const double invnorm = 1 / std::sqrt(xp * xp + yp * yp + zp * zp);
        // Fill normalized world coordinates (W is out parameter)
        W[0] = invnorm * xp;
        W[1] = invnorm * yp;
        W[2] = invnorm * zp;
    }
};

// WarpingSIMD / SIMDImage namespace
namespace ns_simd {
/**
 * @brief Transforms raw input SIMDImage (camera image) to panoramic image
 * @note This include relevant masking operations to remove invalid regions from
 * the transformed image (e.g. those outside of [depr,elev]).
 * @note Pixels transformed outside of the original source image are replaced
 * with a border value (default max of the image Type, coinciding with the
 * invalid pixel mask but afterwards removed to give way for the actual mask:
 * most of the SIMDImageFunctions do not allow stray invalid pixels in the
 * middle or at the bottom of the image).
 * @tparam Type SIMDImage pixel value type
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @tparam AddOn SIMDImage addon of source image (ignored, not Panorama)
 * @param src Input SIMDImage to unroll
 * @param dst Destination image to put the transformed into (will be resized)
 * @param width Target width of the panoramic image
 * @param start Start angle [rad] of the panoramic image (relative to view
 * direction)
 * @param depr Depression angle [rad] (amount below the horizon to be included)
 * @param elev Elevation angle [rad] (amount above the horizon to be included)
 * @param model Camera calibration model (see OCamModel)
 * @param T Homogenous coordinate transform applied to 3d world coordinates
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template<class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
    void
    unroll(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
           SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst,
           const int width, const double start, const double depr,
           const double elev, const OCamModel &model,
           const transform::Matrix &T = transform::Id,
           const Interpolation interpolation = Interpolation::NEAREST,
           const Type border = SIMDTypeInfo<Type>::max()) {
        // Compute height of mapping
        int height = (int) std::lround(
            (width * (fabs(elev) + fabs(depr))) / (2 * M_PI));
        // Panorama: Compute vertical resolution (angle per pixel)
        double vr = (fabs(elev) + fabs(depr)) / height;
        // Panorama: Compute horizon
        double h = elev / vr;
        // Prepare the panorama addon
        Panorama panorama(vr, h);
        // Insert addon into destination image
        dst.addOn = panorama;
        // Set up coordinate mappings of appropriate size
        Remap::Map mapx(width, height), mapy(width, height);
        // Set up mask image of appropriate size
        SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama>
            mask(width, height, panorama);
        // Set all pixels invalid
        mask.clear();
        // Iterate the target image size in row-major order filling mask and
        // maps
        for (int i = 0; i < height; ++i) {
            // 3d world coordinates in spherical coordinates: First angle
            //  NOTE: Placed in outer loop to not recompute for each pixel
            const double alpha = -std::fabs(elev) + (i + 0.5) * vr;
            const double ca = cos(alpha), sa = sin(alpha);
            // Inner loop over image columns
            for (int j = 0; j < width; ++j) {
                // 3d world coordinates in spherical coordinates: Second angle
                const double beta = (j + 0.5) * vr + start;
                const double cb = cos(beta), sb = sin(beta);
                // Transform spherical to cartesian coordinates (unit vector)
                const double tX = std::fabs(ca) * sb;
                const double tY = std::fabs(ca) * cb;
                const double tZ = sa;
                // Make the transform:: Matrix-Vector * operator available
                using transform::operator*;
                // Apply homogenous coordinate transform to world coordinates
                auto W = T * transform::Vector{{tY, tX, tZ, 1.0}};
                // 2d vector to put image coordinates into
                double C[2] = {0, 0};
                // Project world coordinates to camera (image coordinates) using
                // calibrated camera model
                model.world2cam((double *) &W, C);
                // Compute the polar angle of the transformed coordinates
                const double polar = std::atan2(
                    W[2], std::sqrt(W[0] * W[0] + W[1] * W[1])
                );
                // If the polar angle is within [depr,elev], this is a valid
                // pixel mapping
                if (polar < depr && polar > -elev) {
                    // Fill coordinate mappings with image coordinate
                    mapx(i, j) = C[0];
                    mapy(i, j) = C[1];
                    // Mark pixel as valid
                    mask[i][j] = SIMDTypeInfo<Type>::max();
                }
            }
        }
        // Apply the remapping transformation from source to destination image
        remap(src, dst, mapx, mapy, interpolation, border);
        // Max value of image value type might be interpreted as masked
        Type max = SIMDTypeInfo<Type>::max();
        // Remove potentially masked pixels (coincidentally having
        // the mask value or due to interpolation artifacts)
        substituteValue(dst, max, (Type) (max - 1), dst);
        // Blend in the mask image
        maskImage(dst, mask, dst);
    }
}

#endif /*OCAM_MODEL_H*/
