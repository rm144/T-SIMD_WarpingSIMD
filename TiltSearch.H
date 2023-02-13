// ===========================================================================
//
// TiltSearch.H --
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

#ifndef TILT_SEARCH_H
#define TILT_SEARCH_H

// Dynamic memory management
#include <memory>
// Pair of two objects: std::pair
#include <utility>

// Strong type angular units
#include "Angle.H"

// Generic search (minimizer)
#include "Search.H"

// Tilt parameter models
#include "TiltModel.H"
// Tilt correction methods
#include "TiltCorrection.H"

// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"
// SIMD Image type
#include "SIMDImage.H"
// Functions operating on SIMDImages
#include "SIMDImageFunctions.H"
// Bundle of entire warping/min-warping functionality
#include "WarpingBundle.H"

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {}
// Export a short, flat name to this namespace
namespace tilt = ns_simd::tilt;

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {
/**
 * @brief Generic tilt search interface: Combines warping, tilt correction and
 * search strategy
 * @tparam Model Tilt parameter model to use
 * @tparam Warping WarpingBundle to wrap with the tilt search (this yields the
 * search objective)
 */
template <class Model, class Warping>
class Search
{
  // Warping parameter bundle
  using Parameters = typename Warping::Parameters;
  // Warping result bundle
  using Result = typename Warping::Result;
  // Warping Image type
  using ImageType = typename Warping::ImageType;
  // Short typename of the SIMDImage to use with the warping instance
  using Image = typename Warping::Image;
  // Signature of image preprocessing
  using Preprocess = std::function<Image(Image &)>;
  // Fix image type for tilt correction
  using Correction = ns_simd::tilt::Correction<Image>;

  // Image processing identity function (does nothing, hopefully
  // optimized)
  static inline Image &id(Image &image) { return image; }

  // Pure virtual tilt correction interface definition. This is a hidden
  // implementation detail used to create the runtime polymorphic
  // behavior.
  struct Interface
  {
    // Pure virtual tilt correction method
    virtual std::pair<Result, Model> run(Image &, Image &, const Parameters &,
                                         Preprocess, Interpolation,
                                         ImageType) = 0;
  };

  // Implementation of (optionally constrained) tilt search given a
  // generic search method and optional constraint function
  template <class Searcher, class Constraint = search::Unconstrained>
  struct Wrapper final : public Interface
  {
    // Tilt correction solution to use
    Correction correction;
    // Search strategy to use
    Searcher searcher;
    // Wrapped warping instance
    Warping &warping;
    // Constraint function object
    Constraint constraint;

    // Wraps tilt correction, search strategy, warping bundle (by
    // reference!) and (optional) constraint
    Wrapper(Correction c, Searcher s, Warping &w, Constraint g)
      : correction {c}, searcher {s}, warping {w}, constraint {g}
    {}

    // Core implementation of the warping tilt correction search
    std::pair<Result, Model> run(Image &ss, Image &cv, const Parameters &params,
                                 Preprocess pre,
                                 const Interpolation interpolation,
                                 const ImageType border) override
    {
      // Construct the warping tilt search objective
      auto objective = [&](auto &&...xs) {
        // Temporary for tilt corrected and preprocessed images
        Image tmp, pre_ss, pre_cv;
        // Apply tilt correction to current view image
        correction(cv, tmp, Model(xs...), interpolation, border);
        // Preprocess both images
        pre_ss = pre(ss), pre_cv = pre(tmp);
        // Warping results target object (output parameters
        // of run)
        Result result;
        // Run warping on the image pair
        warping.run(pre_ss, pre_cv, params, result);
        // Yield warping result as objective function value
        return result;
      };
      // Minimize the warping tilt search objective with
      // constraint
      auto [min, argmin] = search::minimize(objective, searcher, constraint);
      // Convert argmin tuple to requested tilt Model
      return std::pair {min, Model(argmin)};
    }
  };

  /**
   * @brief Constructs a shared pointer to a tilt search implementation
   * wrapper
   * @tparam Searcher Search method type
   * @tparam Constraint Constraint function type
   * @param c Tilt correction selector
   * @param s Search strategy
   * @param w Warping bundle (by reference!)
   * @param g Constraint function (or ::Search::Unconstrained{})
   * @return Shared pointer to implementation wrapper
   */
  template <class Searcher, class Constraint>
  auto wrap(Correction c, Searcher s, Warping &w, Constraint g)
  {
    return std::make_shared<Wrapper<Searcher, Constraint>>(c, s, w, g);
  }

  /**
   * @brief Constructs a shared pointer to a tilt search implementation
   * wrapper
   * @tparam Searcher Search method type
   * @tparam Constraint Constraint function type
   * @param c Tilt correction selector
   * @param s Search strategy
   * @param w Warping bundle (by reference!)
   * @param g Constraint function (or ::Search::Unconstrained{})
   * @return Shared pointer to implementation wrapper
   */
  template <class Searcher, class Constraint>
  auto wrap(Selector c, Searcher s, Warping &w, Constraint g)
  {
    return wrap(make_correction<Image>(c), s, w, g);
  }

public:
  /**
   * @brief Initializes the warping tilt search by wrapping a search
   * strategy and warping bundle
   * @note The first (unnamed) parameter allows for template argument
   * deduction on constructor call: Either all or none might be deduced.
   * @param correction Tilt correction solution to use (function object)
   * @param searcher Search strategy to use for tilt search
   * @param warping Warping bundle to use as part of the search objective
   * @param constraint Optional constraint function passed to the searcher
   */
  template <class Searcher, class Constraint = search::Unconstrained>
  Search(Model, Correction correction, Searcher searcher, Warping &warping,
         Constraint constraint = {})
    : ptr {wrap(correction, searcher, warping, constraint)}
  {}

  /**
   * @brief Initializes the warping tilt search by wrapping a search
   * strategy and warping bundle
   * @note The first (unnamed) parameter allows for template argument
   * deduction on constructor call: Either all or none might be deduced.
   * @note This constructs the tilt correction given the selector number
   * using the make_correction\<Image\> factory.
   * @param correction Tilt correction solution selector
   * @param searcher Search strategy to use for tilt search
   * @param warping Warping bundle to use as part of the search objective
   * @param constraint Optional constraint function passed to the searcher
   */
  template <class Searcher, class Constraint = search::Unconstrained>
  Search(Model, Selector correction, Searcher searcher, Warping &warping,
         Constraint constraint = {})
    : ptr {wrap(correction, searcher, warping, constraint)}
  {}

  /**
   * @brief Runs a warping tilt search for ss-cv image pair
   * @param ss Snapshot image (SIMDImage appropriate for the wrapped
   * warping bundle)
   * @param cv Current view image (SIMDImage appropriate for the wrapped
   * warping bundle)
   * @param params Warping parameters (bundled)
   * @param pre Image preprocessing function applied AFTER the tilt
   * @param interpolation Interpolation method to use for remapping
   * non-integer coordinates
   * @param border Border value to insert if a pixel is outside the image
   * @return A std::pair of the warping result and the minimizing tilt
   * parameters in terms of the requested Model
   */
  std::pair<Result, Model> run(
    Image &ss, Image &cv, const Parameters &params, Preprocess pre = id,
    const Interpolation interpolation = Interpolation::NEAREST,
    const ImageType border            = SIMDTypeInfo<ImageType>::max())
  {
    // Delegate to runtime polymorphic implementation
    return ptr->run(ss, cv, params, pre, interpolation, border);
  }

  // Pointer to polymorphic object implementing the tilt correction
  // interface
  std::shared_ptr<Interface> ptr;
};
} // namespace ns_simd::tilt

#endif /*TILT_SEARCH_H*/
