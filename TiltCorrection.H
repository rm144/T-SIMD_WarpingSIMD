// ===========================================================================
//
// TiltCorrection.H --
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

#ifndef TILT_CORRECTION_H
#define TILT_CORRECTION_H

// Math functions for image transformations
#include <cmath>
// General-purpose polymorphic function wrapper: std::function
#include <functional>
// operator>>
#include <istream>

// Strong type angular units
#include "Angle.H"

// Tilt parameter models
#include "TiltModel.H"

// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"
// SIMD Image type
#include "SIMDImage.H"
// Functions operating on SIMDImages
#include "SIMDImageFunctions.H"

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {}
// Export a short, flat name to this namespace
namespace tilt = ns_simd::tilt;

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {
/**
 * @brief Exact tilt correction solution tag type (compile time selector)
 */
struct Exact
{};

/**
 * @brief Approximate tilt correction solution tag type (compile time selector)
 */
struct Approximate
{};

/**
 * @brief Vertical tilt correction solution tag type (compile time selector)
 */
struct Vertical
{};

/**
 * @brief Tilt correction solution runtime selector number
 * @see Factory function make_correction(Selector)
 */
enum Selector : int { EXACT = 0, APPROXIMATE = 1, VERTICAL = 2 };

/**
 * @brief Stream extracts a tilt selector number
 * @param lhs Input stream to extract the value from
 * @param rhs Tilt selector to insert the value into
 * @return Reference to the modified input stream after extraction
 */
std::istream &operator>>(std::istream &lhs, Selector &rhs)
{
  return lhs >> (int &) rhs;
}

/**
 * @brief Applies tilt correction to panoramic SIMDImage
 * @note This is kind of a hub template selecting the solution overload via tag
 * dispatching.
 * @tparam Solution Tilt correction solution selector specifying the
 * approximation to use
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @param src Input image to tilt correct
 * @param dst Output image to write the tilt corrected into (will be resized)
 * @param tilt Tilt in Axis-Angle representation to apply to the image
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Solution, class Type, int SIMD_WIDTH, int SIMD_ALIGN>
void correct(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &src,
             SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst,
             const tilt::AxisAngle<Radian> &tilt,
             const Interpolation interpolation = Interpolation::NEAREST,
             const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Call solution overload using tag dispatching
  correct(Solution {}, src, dst, tilt, interpolation, border);
}

/**
 * @brief Applies the exact solution tilt correction to panoramic SIMDImage
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @param src Input image to tilt correct
 * @param dst Output image to write the tilt corrected into (will be resized)
 * @param tilt Tilt in Axis-Angle representation to apply to the image
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
void correct(Exact, // Tag dispatch
             const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &src,
             SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst,
             const tilt::AxisAngle<Radian> &tilt,
             const Interpolation interpolation = Interpolation::NEAREST,
             const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Short names for horizon and vertical resolution of th input panoramic
  // image. These do not change and will be inserted into the destination
  // Panorama addon.
  const double h  = src.addOn.horizon;
  const double vr = src.addOn.verticalResolution;
  // Compute horizontal resolution used for converting between column
  // angles and pixel index.
  const double hr = (2 * M_PI) / src.w;
  // Convert (horizontal) rotation angle to pixel
  const double r = (double) tilt.theta_R / hr;
  // Compute sine and cosine of the tilt angle (avoid recomputing this in
  // the inner loop)
  double sp = sin((double) tilt.phi), cp = cos((double) tilt.phi);
  // Table holding sine and cosine of horizontal angle
  std::vector<double> St(src.w), Ct(src.w);
  // Precompute sine and cosine of horizontal angle to avoid
  // recomputing per inner loop iteration
  for (int j = 0; j < src.w; ++j) {
    // Convert column index [pixel] to angle [rad] and shift by
    // horizontal rotation angle
    const double theta = j * hr - (double) tilt.theta_R;
    // Fill table with sine and cosine of angle
    St[j] = sin(theta);
    Ct[j] = cos(theta);
  }
  // Set up row and column coordinate mappings of appropriate size.
  //  NOTE: Output will have the same size as the input
  ns_simd::Remap::Map mapx(src.w, src.h), mapy(src.w, src.h);
  // Process mappings pixel-wise in row-major order
  for (int i = 0; i < src.h; ++i) {
    // Vertical angular coordinate
    const double delta = (h - i) * vr;
    // Compute sine and cosine of the vertical angle (avoid recomputing
    // this in the inner loop)
    const double sd = sin(delta), cd = cos(delta);
    // Inner loop: columns
    for (int j = 0; j < src.w; ++j) {
      // Get precomputed sine and cosine from table
      const double st = St[j], ct = Ct[j];
      // Exact solution row coordinate mapping
      //  NOTE: Combines exact solution transform an angles asin(...)
      //  and transformation back to pixel coordinates (.) / vr + h.
      double mi = -asin(-sp * st * cd + cp * sd) / vr + h;
      // Exact solution columns coordinate mapping
      //  NOTE: Combines exact solution transform an angles atan2(...)
      //  and transformation back to pixel coordinates (.) / hr + r.
      double mj = atan2(cp * st * cd + sp * sd, ct * cd) / hr + r;
      // Remap from [-pi,+pi] to [0,+2pi] ([-w/2,+w/2] to [0,w])
      if (mj < 0) { mj += src.w; }
      // Fill coordinate mapping tables
      mapx(i, j) = mi;
      mapy(i, j) = mj;
    }
  }
  // Insert panorama parameters into destination image addon
  //  NOTE: These are the same as extracted from the input image
  dst.addOn = Panorama(vr, h);
  // Apply the remapping transform to the input image
  ns_simd::remap(src, dst, mapx, mapy, interpolation, border);
}

/**
 * @brief Applies the approximate solution tilt correction to panoramic
 * SIMDImage
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @param src Input image to tilt correct
 * @param dst Output image to write the tilt corrected into (will be resized)
 * @param tilt Tilt in Axis-Angle representation to apply to the image
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
void correct(Approximate, // Tag dispatch
             const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &src,
             SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst,
             const tilt::AxisAngle<Radian> &tilt,
             const Interpolation interpolation = Interpolation::NEAREST,
             const Type border                 = SIMDTypeInfo<Type>::max())
{
  // NOTE: The following is almost identical to the exact solution (see
  // above). The difference is in the inner loop computing the mapping
  // based on the approximate solution. However, outsourcing common parts
  // does not really make sense as this is primarily variable/constant and
  // table setup code.

  // Short names for horizon and vertical resolution of th input panoramic
  // image. These do not change and will be inserted into the destination
  // Panorama addon.
  const double h  = src.addOn.horizon;
  const double vr = src.addOn.verticalResolution;
  // Compute horizontal resolution used for converting between column
  // angles and pixel index.
  const double hr = (2 * M_PI) / src.w;
  // Compute sine and cosine of the tilt angle (avoid recomputing this in
  // the inner loop)
  double sp = sin((double) tilt.phi);
  // Table holding sine and cosine of horizontal angle
  std::vector<double> St(src.w), Ct(src.w);
  // Precompute sine and cosine of horizontal angle to avoid
  // recomputing per inner loop iteration
  for (int j = 0; j < src.w; ++j) {
    // Convert column index [pixel] to angle [rad] and shift by
    // horizontal rotation angle
    const double theta = j * hr - (double) tilt.theta_R;
    // Fill table with sine and cosine of angle
    St[j] = sin(theta);
    Ct[j] = cos(theta);
  }
  // Set up row and column coordinate mappings of appropriate size.
  //  NOTE: Output will have the same size as the input
  ns_simd::Remap::Map mapx(src.w, src.h), mapy(src.w, src.h);
  // Process mappings pixel-wise in row-major order
  for (int i = 0; i < src.h; ++i) {
    // Vertical angular coordinate
    const double delta = (h - i) * vr;
    // Inner loop: columns
    for (int j = 0; j < src.w; ++j) {
      // Get precomputed sine and cosine from table
      const double st = St[j], ct = Ct[j];
      // Approximate solution row coordinate mapping
      //  NOTE: sp * st / vr is the angular transform converted back
      //  to pixels, the effect on depending on 'i' is linearized and
      //  thus does not even be converted to angle and back.
      double mi = i + sp * st / vr;
      // Exact solution columns coordinate mapping
      //  NOTE: sp * delta * ct / hr is the angular transform
      //  converted back to pixels, the effect on depending on 'j' is
      //  linearized and thus does not even be converted to angle and
      //  back.
      double mj = j + sp * delta * ct / hr;
      // Remap from [-pi,+pi] to [0,+2pi] ([-w/2,+w/2] to [0,w])
      if (mj < 0) { mj += src.w; }
      // Fill coordinate mapping tables
      mapx(i, j) = mi;
      mapy(i, j) = mj;
    }
  }
  // Insert panorama parameters into destination image addon
  //  NOTE: These are the same as extracted from the input image
  dst.addOn = Panorama(vr, h);
  // Apply the remapping transform to the input image
  ns_simd::remap(src, dst, mapx, mapy, interpolation, border);
}

/**
 * @brief Applies the vertical solution tilt correction to panoramic SIMDImage
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @param src Input image to tilt correct
 * @param dst Output image to write the tilt corrected into (will be resized)
 * @param tilt Tilt in Axis-Angle representation to apply to the image
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
void correct(Vertical, // Tag dispatch
             const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &src,
             SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst,
             const tilt::AxisAngle<Radian> &tilt,
             const Interpolation interpolation = Interpolation::NEAREST,
             const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Short names for horizon and vertical resolution of th input panoramic
  // image. These do not change and will be inserted into the destination
  // Panorama addon.
  const double h  = src.addOn.horizon;
  const double vr = src.addOn.verticalResolution;
  // Compute horizontal resolution used for converting between column
  // angles and pixel index.
  const double hr = (2 * M_PI) / src.w;
  // Compute sine of the tilt angle (avoid recomputing this in the inner
  // loop)
  double sp = sin((double) tilt.phi);
  // Set up column shift offset table of image size
  ns_simd::Shift::ShiftCol shift(src.w);
  // Iterate the columns to fill the mapping
  for (int j = 0; j < src.w; ++j) {
    // Convert column index [pixel] to angle [rad] and shift by
    // horizontal rotation angle, then compute mapping, convert back
    // to [pixel] and insert into offset table
    shift(j) = sp * sin((double) j * hr - (double) tilt.theta_R) / vr;
  }
  // Insert panorama parameters into destination image addon
  //  NOTE: These are the same as extracted from the input image
  dst.addOn = Panorama(vr, h);
  // Apply the shift remapping transform to the input image
  ns_simd::shift(src, dst, shift, interpolation, border);
}

/**
 * @brief Type erased, runtime polymorphic tilt correction interface
 * @tparam Image SIMDImage type to be processed
 */
template <class Image>
class Correction;

/**
 * @brief Type erased, runtime polymorphic tilt correction interface
 * @tparam Type Value type of SIMDImage to be processed
 * @tparam SIMD_WIDTH Width of SIMD vectors of SIMDImage to be processed
 * @tparam SIMD_ALIGN Align of SIMD vectors of SIMDImage to be processed
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
class Correction<SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> >
{
  // Shorter type names for images
  using Image = SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama>;

  // Tilt correction identity function, transfers image from src to dst.
  static inline void id(const Image &src, Image &dst, AxisAngle<Radian>,
                        Interpolation, Type)
  {
    dst = src;
  }

public:
  /**
   * @brief Default initializes tilt correction with identity function.
   */
  Correction() : wrapped(id) {}

  /**
   * @brief Constructs tilt correction by wrapping arbitrary object
   * realising the tilt-correction interface
   * @tparam Wrapped Type of the wrapped implementation
   * @param wrapped The implementation to wrap
   */
  template <class Wrapped>
  Correction(Wrapped wrapped) // NOLINT: Intentional implicit
    : wrapped {std::move(wrapped)}
  {}

  /**
   * @brief Applies tilt correction to panoramic SIMDImage
   * @param src Input image to tilt correct
   * @param dst Output image to write the tilt corrected into
   * @param tilt Tilt in Axis-Angle representation to apply to the image
   * @param interpolation Interpolation method to use for remapping
   * non-integer coordinates
   * @param border Border value to insert if a pixel is outside the image
   */
  void operator()(const Image &src, Image &dst,
                  const tilt::AxisAngle<Radian> tilt,
                  const Interpolation interpolation = Interpolation::NEAREST,
                  const Type border = SIMDTypeInfo<Type>::max()) const
  {
    // Delegate to the wrapped, polymorphic implementation
    wrapped(src, dst, tilt, interpolation, border);
  }

private:
  // Wrapped function realizing the tilt correction interface
  std::function<void(const Image &, Image &, AxisAngle<Radian>, Interpolation,
                     Type)>
    wrapped;
};

/**
 * @brief Constructs tilt correction function appropriate to process the given
 * Image type
 * @tparam Type Value type of SIMDImage to be processed
 * @tparam SIMD_WIDTH Width of SIMD vectors of SIMDImage to be processed
 * @tparam SIMD_ALIGN Align of SIMD vectors of SIMDImage to be processed
 * @param selector Tilt correction solution selector
 * @return Tilt correction function (polymorphic wrapper) operating on Image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
Correction<SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> > make_correction(
  const Selector selector,
  SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> = {})
{
  switch (selector) {
  case EXACT: return correct<Exact, Type, SIMD_WIDTH, SIMD_ALIGN>;
  case APPROXIMATE: return correct<Approximate, Type, SIMD_WIDTH, SIMD_ALIGN>;
  case VERTICAL: return correct<Vertical, Type, SIMD_WIDTH, SIMD_ALIGN>;
  default:
    // Adhoc exception class
    class Exception : public std::exception
    {
      // Returns hardcoded message
      const char *what() const noexcept override
      { // NOLINT
        return "Tilt Correction Method Not Implemented";
      }
    };
    // Throw the adhoc exception
    throw Exception();
  }
}

/**
 * @brief Constructs tilt correction function appropriate to process the given
 * Image type
 * @tparam Image Type of SIMDImage to be processed
 * @param selector Tilt correction solution selector
 * @return Tilt correction function (polymorphic wrapper) operating on Image
 */
template <class Image>
Correction<Image> make_correction(const Selector selector)
{
  return make_correction(selector, Image {});
}
} // namespace ns_simd::tilt

#endif /*TILT_CORRECTION_H*/
