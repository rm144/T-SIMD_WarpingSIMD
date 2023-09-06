// ===========================================================================
//
// Setup.H --
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

#pragma once
#ifndef EXP_SETUP_H
#define EXP_SETUP_H

// All SIMD vectors and corresponding definitions
#include "SIMDVecAll.H"
// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"
// Bundle of entire warping/min-warping functionality
#include "WarpingBundle.H"

// Select width and alignment of SIMD vectors (and images)
inline constexpr int SW = NATIVE_SIMD_WIDTH;
inline constexpr int SA = NATIVE_SIMD_WIDTH;

// Type alias to ns_simd::WarpingBundle with preconfigured Types and width
using WarpingBundle = ns_simd::WarpingBundle<
  // Type of pixels in input images
  ns_simd::SIMDByte,
  // Type used in preprocessing step
  ns_simd::SIMDFloat,
  // Type used in distance measure step
  ns_simd::SIMDShort,
  // Type used in scale plane stack
  ns_simd::SIMDByte,
  // Type used for matches
  ns_simd::SIMDWord,
  // Type used for compass estimate / compass acceleration
  ns_simd::SIMDInt,
  // Width of SIMD vectors
  SW,
  // Alignment of SIMD vectors
  SA>;

// Standard C++ time strong types (e.g. microseconds)
#include <chrono>
// Execution time measurement via struct timespec
#include "TimeMeasurement.H"

// Augments Type with time and count measuring the calls of the run method
template <class Type>
class RunStatistics : public Type
{
  // Counter of Type::run calls
  std::size_t counter = 0;
  // Total (summed) runtime of Type::run calls
  std::chrono::microseconds time {0};

public:
  /**
   * @brief Returns the number of run calls so far
   * @return Number of Type::run calls
   */
  [[nodiscard]] std::size_t count() const { return counter; }

  /**
   * @brief Returns the total runtime accumulated so far
   * @return Total runtime in microseconds
   */
  [[nodiscard]] double total_time() const { return time.count(); }

  /**
   * @brief Returns the average runtime accumulated so far
   * @return Average runtime in microseconds
   */
  [[nodiscard]] double avg_time() const { return (time / counter).count(); }

  // Inherit the Type constructor
  using Type::Type;

  // Wraps the Type::run method forwarding the arguments after
  // incrementing the counter and measuring the execution time.
  template <class... Args>
  void run(Args &&...args)
  {
    // Start measuring time
    auto t0 = ns_simd::getTimeSpec();
    // Count and run
    counter++, Type::run(std::forward<Args>(args)...);
    // Stop measuring time
    auto t1 = ns_simd::getTimeSpec();
    // Accumulate time
    time +=
      std::chrono::microseconds {(long) ns_simd::timeSpecDiffUsec(t1, t0)};
  }
};

// Configuration parameters from environment variables
#include "Env.H"

/**
 * @brief Constructs a WarpingBundle augmented by RunStatistics with parameters
 * from the env_config.
 * @param width Width of SIMDImages to be processed by the warping bundle
 * @return Configured instance of RunStatistics\<WarpingBundle\>
 */
auto make_warping(const int width)
{
  // Number of alpha and psi search steps
  ENV_INT(nAlpha, 128);
  ENV_INT(nPsi, 128);
  // Number of scale planes
  ENV_INT(nScalePlanes, 9);
  // Maximal scale factor
  ENV_DOUBLE(maxScaleFactor, 2.0);
  // Maximal threshold
  ENV_DOUBLE(maxThresholdMinWarping, 2.5);
  // Number of rho steps
  ENV_INT(nRhoWarping, 20);
  // Restriction of min-warping template to certain values of rho; (0.0, 100.0)
  // means that there is no restriction.
  ENV_DOUBLE(rhoMinMinWarping, 0.0);
  ENV_DOUBLE(rhoMaxMinWarping, 100.0);
  ENV_DOUBLE(rhoMinWarping, 0.0);
  ENV_DOUBLE(rhoMaxWarping, 1.0);
  // First phase method selector
  ENV_INT(firstPhase, 23);
  // Harris Parameters
  ENV_DOUBLE(harrisK, 0.04);
  ENV_INT(harrisBinomialFilterCount, 1);
  // Sigmoid weights
  ENV_DOUBLE_VEC(sigmoidW, );
  ENV_DOUBLE_VEC(sigmoidW0, );
  // MinWarping search method (search loop order)
  ENV_INT(searchMethodMinWarping, 0);
  // Warping search method (search loop order)
  ENV_INT(searchMethodWarping, 0);
  // Construct an instance of WarpingBundle
  return RunStatistics<WarpingBundle> {
    width,
    nAlpha,
    nPsi,
    nScalePlanes,
    maxScaleFactor,
    maxThresholdMinWarping,
    rhoMinMinWarping,
    rhoMaxMinWarping,
    rhoMinWarping,
    rhoMaxWarping,
    nRhoWarping,
    firstPhase,
    searchMethodWarping,
    searchMethodMinWarping,
    ns_simd::HarrisParameter(harrisK, harrisBinomialFilterCount),
    sigmoidW,
    sigmoidW0};
}

/**
 * @brief Constructs a WarpingBundle parameter set with parameters from the
 * env_config.
 * @return A WarpingBundle::Parameters object ready to use with the warping
 * bundle run method.
 */
auto make_warping_params()
{
  // First phase interpolation method
  ENV_INT(firstPhaseInterpolation, 0);
  // Pixel scale parameters
  ENV_DOUBLE_VEC(pixelScale, 15.427);
  // Post scale parameter
  ENV_DOUBLE(postScale, 170.664);
  // Warping (1) or MinWarping (0)?
  ENV_INT(searcher, 0);
  // Switch for compass acceleration (1=sum, 2=quantil)
  ENV_INT(compassAcceleration, 1);
  // Quantil fraction for quantil-based compass acceleration
  ENV_DOUBLE(quantilCompassFraction, 0.2);
  // Parameter for compass acceleration
  ENV_DOUBLE(psiFraction, 0.3);
  // Switch for partial search (otherwise full search)
  ENV_INT(partialSearch, 0);
  // Partial searcher, but full search range (tests only)
  ENV_INT(partialFullSearch, 0);
  // Search radii for partial search
  ENV_DOUBLE(alphaRad, M_PI / 4);
  ENV_DOUBLE(psiRad, M_PI / 4);
  // Switch for double search (second run with exchanged snapshot/current view)
  // affects both partial and full search
  ENV_INT(doubleSearch, 1);
  // Fine search is executed
  ENV_INT(fineSearch, 0);
  // Fine search uses double search
  ENV_INT(doubleSearchFine, 1);
  // Fine search only uses immediate surroundings
  ENV_INT(fineSearchSurround, 0);
  // Construct an instance of WarpingBundle::Parameters
  return WarpingBundle::Parameters {firstPhaseInterpolation,
                                    pixelScale,
                                    postScale,
                                    0.0,
                                    0.0,
                                    searcher,
                                    compassAcceleration,
                                    psiFraction,
                                    quantilCompassFraction,
                                    partialSearch,
                                    alphaRad,
                                    psiRad,
                                    partialFullSearch,
                                    doubleSearch,
                                    fineSearch,
                                    doubleSearchFine,
                                    fineSearchSurround};
}

#endif /*EXP_SETUP_H*/
