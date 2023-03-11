// ===========================================================================
//
// WarpingBundle.H --
// bundle of entire warping/min-warping functionality
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

// https://
// en.wikibooks.org/wiki/C%2B%2B_Programming/Code/Design_Patterns#Factory

#ifndef WARPING_BUNDLE_H_
#define WARPING_BUNDLE_H_

#include "TimeMeasurement.H"
#include "WarpingCompound.H"
#include "WarpingFactories.H"

#include <ostream>
#include <vector>

namespace simd {

template <typename ImgType, typename ProcType, typename MeasType,
          typename SPSType, typename MatchType, typename CompassType,
          int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingBundle
  : public WarpingSearcherSelector
  , public ComplexSearchSelector
{
public:
  WarpingCompound<SPSType, MatchType, CompassType, SIMD_WIDTH, SIMD_ALIGN> *wsc;
  WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                        SIMD_ALIGN> *spsComp;
  MinWarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
    *fullSearcherMinWarping;
  MinWarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
    *partialSearcherMinWarping;
  WarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN> *fullSearcherWarping;
  WarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
    *partialSearcherWarping;

  WarpingBundle(int w, int nAlpha, int nPsi, int nScalePlanes,
                double maxScaleFactor, double maxThresholdMinWarping,
                double rhoMinMinWarping, double rhoMaxMinWarping,
                double rhoMinWarping, double rhoMaxWarping, int nRhoWarping,
                int firstPhase, int searchMethodWarping,
                int searchMethodMinWarping, const HarrisParameter &harrisPar,
                const std::vector<double> &sigmoidW,
                const std::vector<double> &sigmoidW0)
  {
    // create first phase computation
    spsComp =
      WarpingSPSComputationFactory<ImgType, ProcType, MeasType, SPSType,
                                   SIMD_WIDTH, SIMD_ALIGN>::create(firstPhase,
                                                                   harrisPar,
                                                                   sigmoidW,
                                                                   sigmoidW0);
    // the same type of searcher (e.g. PAXY) is created for both the
    // full and the partial search
    // Warping
    fullSearcherWarping =
      WarpingFullFactory<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>::create(
        searchMethodWarping);
    partialSearcherWarping =
      WarpingPartialFactory<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>::create(
        searchMethodWarping);
    // MinWarping
    fullSearcherMinWarping =
      MinWarpingFullFactory<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>::create(
        searchMethodMinWarping);
    partialSearcherMinWarping =
      MinWarpingPartialFactory<SPSType, MatchType, SIMD_WIDTH,
                               SIMD_ALIGN>::create(searchMethodMinWarping);
    // create compound
    wsc = new WarpingCompound<SPSType, MatchType, CompassType, SIMD_WIDTH,
                              SIMD_ALIGN>(
      w, nAlpha, nPsi, nScalePlanes, maxScaleFactor, maxThresholdMinWarping,
      rhoMinMinWarping, rhoMaxMinWarping, rhoMinWarping, rhoMaxWarping,
      nRhoWarping);
  }

  // TODO: Not copiable due to resource management unclear
  WarpingBundle(const WarpingBundle &) = delete;

  ~WarpingBundle()
  {
    delete wsc;
    delete spsComp;
    delete fullSearcherMinWarping;
    delete partialSearcherMinWarping;
    delete fullSearcherWarping;
    delete partialSearcherWarping;
  }

  // compute best match between ss and cv under azimuthal rotation
  // (compass), can be used for localization
  // adapted from code provided by Michael Horst
  // 27. Feb 18 (rm): vert.res. and horizon taken from images
  void compassMatch(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
                    SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
                    int interpolation, const std::vector<double> &pixelScale,
                    double postScale, double &dMin, double &t1)
  {
    struct timespec ts0, ts1;
    ts0 = getTimeSpec();
    wsc->computeSPS(ss, cv, interpolation, pixelScale, postScale, *spsComp);
    SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> compass;
    wsc->spsArray[SPS_ORIG]->compassEstimate(compass);
    dMin = min(compass);
    ts1  = getTimeSpec();
    t1   = timeSpecDiffUsec(ts1, ts0);
  }

  // run warping / min-warping computation; compute home vector and
  // some additional information
  // 27. Feb 18 (rm): vert.res. and horizon taken from images
  void run(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
           SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
           int interpolation, const std::vector<double> &pixelScale,
           double postScale, double alphaEst, double psiEst, int searcher,
           int compassAcceleration, double psiFraction,
           double quantilCompassFraction, int partialSearch, double alphaRad,
           double psiRad, int partialFullSearch, int doubleSearch,
           int fineSearch, int doubleSearchFine, int fineSearchSurround,
           double &alphaMin, double &psiMin, double &dMin, int &iAlphaMin,
           int &iPsiMin, int &iAlphaMinInv, int &iPsiMinInv, int &jAlphaOff,
           int &jPsiOff, SIMDImage<MatchType, 1, 1> &fineMatch, double &t1,
           double &t2)
  {
    // ----- first phase -----
    // signed long long t0;
    // t0 = getTimeUsec();
    struct timespec ts0, ts1, ts2;
    ts0 = getTimeSpec();
    // <ProcType,MeasType> not needed, defined through spsComp
    switch (searcher) {
    case minWarpingSearcher:
    case warpingSearcher:
      wsc->computeSPS(ss, cv, interpolation, pixelScale, postScale, *spsComp);
      break;

    case minWarpingSearcher1SP:
    case warpingSearcher1SP:
      wsc->computeSPS1SP(ss, cv, interpolation, pixelScale, postScale,
                         *spsComp);
      break;

    default:
      throw SIMDException("WarpingBundle::run", "invalid warping searcher");
    };

    ts1 = getTimeSpec();
    t1  = timeSpecDiffUsec(ts1, ts0);
    // ----- second phase -----
    switch (searcher) {
    case minWarpingSearcher:
      // MinWarping
      wsc->minWarping.search(
        *fullSearcherMinWarping, *partialSearcherMinWarping, alphaEst, psiEst,
        compassAcceleration, psiFraction, quantilCompassFraction, partialSearch,
        alphaRad, psiRad, partialFullSearch, doubleSearch, fineSearch,
        doubleSearchFine, fineSearchSurround, alphaMin, psiMin, dMin, iAlphaMin,
        iPsiMin, iAlphaMinInv, iPsiMinInv, jAlphaOff, jPsiOff, fineMatch);
      break;

    case warpingSearcher:
      // Warping
      wsc->warping.search(
        *fullSearcherWarping, *partialSearcherWarping, alphaEst, psiEst,
        compassAcceleration, psiFraction, quantilCompassFraction, partialSearch,
        alphaRad, psiRad, partialFullSearch, doubleSearch, fineSearch,
        doubleSearchFine, fineSearchSurround, alphaMin, psiMin, dMin, iAlphaMin,
        iPsiMin, iAlphaMinInv, iPsiMinInv, jAlphaOff, jPsiOff, fineMatch);
      break;

    case minWarpingSearcher1SP:
      // MinWarping1SP
      wsc->minWarping1SP.search(
        *fullSearcherMinWarping,    // may change
        *partialSearcherMinWarping, // may change
        alphaEst, psiEst, compassAcceleration, psiFraction,
        quantilCompassFraction, partialSearch, alphaRad, psiRad,
        partialFullSearch, doubleSearch, fineSearch, doubleSearchFine,
        fineSearchSurround, alphaMin, psiMin, dMin, iAlphaMin, iPsiMin,
        iAlphaMinInv, iPsiMinInv, jAlphaOff, jPsiOff, fineMatch);
      break;

    case warpingSearcher1SP:
      // Warping1SP
      wsc->warping1SP.search(
        *fullSearcherWarping,    // may change
        *partialSearcherWarping, // may change
        alphaEst, psiEst, compassAcceleration, psiFraction,
        quantilCompassFraction, partialSearch, alphaRad, psiRad,
        partialFullSearch, doubleSearch, fineSearch, doubleSearchFine,
        fineSearchSurround, alphaMin, psiMin, dMin, iAlphaMin, iPsiMin,
        iAlphaMinInv, iPsiMinInv, jAlphaOff, jPsiOff, fineMatch);
      break;

    default:
      throw SIMDException("WarpingBundle::run", "invalid warping searcher");
    }
    ts2 = getTimeSpec();
    t2  = timeSpecDiffUsec(ts2, ts1);
  }

  // -----------------------------------------------------------------------
  // Condensed interface to get cleaner code for wrapping tilt search
  // -----------------------------------------------------------------------
  // Expose the image type (template argument)
  using ImageType = ImgType;
  // Expose the MatchType (template argument)
  using MATCH_TYPE = MatchType;
  // SIMDImage to use with this bundle instance
  using Image = SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama>;

  /**
   * @brief Bundles WarpingBundle::run parameter structure (with defaults)
   */
  struct Parameters
  {
    // NOTE: Comments taken from warpingSIMDTest.C

    /**
     * @brief First phase (scale-plane stack computation) interpolation
     */
    int interpolation = 0;

    /**
     * @brief Pixel scale parameter for scale plane stack computation
     */
    std::vector<double> pixelScale = {16};

    /**
     * @brief Post scale parameter for scale plane stack computation
     */
    double postScale = 200;

    // TODO: ???
    double alphaEst = 0;
    double psiEst   = 0;

    /**
     * @brief Select search method
     */
    int searcher = ns_simd::WarpingSearcherSelector::minWarpingSearcher;

    /**
     * @brief Second Phase: Kind of compass acceleration to use 1=sum,
     * 2=quantile
     */
    int compassAcceleration = 0;

    /**
     * @brief Parameter for compass acceleration
     */
    double psiFraction = 0.3;

    /**
     * @brief Quantile fraction for quantile-based compass acceleration
     */
    double quantileCompassFraction = 0.2;

    /**
     * @brief Switch for partial search (otherwise full search)
     */
    int partialSearch = 0;

    /**
     * @brief Search radii for partial search
     */
    double alphaRad = M_PI / 4.0;

    /**
     * @brief Search radii for partial search
     */
    double psiRad = M_PI / 4.0;

    /**
     * @brief Partial searcher, but full search range (tests only)
     */
    int partialFullSearch = 0;

    /**
     * @brief Switch for double search (second run with exchanged
     *  snapshot/current view) affects both partial and full search
     */
    int doubleSearch = 0;

    /**
     * @brief Fine search is executed
     */
    int fineSearch = 0;

    /**
     * @brief Fine search uses double search
     */
    int doubleSearchFine = 0;

    /**
     * @brief Fine search only uses immediate surroundings
     */
    int fineSearchSurround = 0;
  };

  /**
   * @brief Bundled WarpingBundle::run results structure
   */
  struct Result
  {
    // Warping angle estimates
    double alpha {}, psi {}; // Maybe strong type?
    // Distance under the warp curve
    // NOTE: Default initializing this to infinity is important to make
    // constrained optimization of the warping objective work:
    // std::numeric_limits<>::infinity() returns a default initialized
    // value in case of custom/non-specialized types.
    double d = std::numeric_limits<double>::infinity();
    // Indices corresponding to warping angle estimates
    int iAlpha {}, iPsi {}, iAlphaInv {}, iPsiInv {}, jAlphaOff {}, jPsiOff {};
    // Fine match
    ns_simd::SIMDImage<MatchType, 1, 1> fineMatch;

    /**
     * @brief Converts the warping result to a scalar value of type
     * double.
     * @note For now, this is sufficient to enable arithmetic and
     * comparison operators via implicit conversion to use the warping
     * result as the objective of the search methods.
     * @return Distance under the warp curve d as double
     */
    operator double() const
    { // NOLINT: Intentional implicit
      return d;
    }

    // TODO: Add proper implementation of arithmetic and comparison
    //  operators. Note however, arithmetic is rather difficult to define
    //  for the angle (alpha, psi) index and fineMatch members...

    /**
     * @brief Inserts warping result into an output stream
     * @note This inserts alpha, psi and d
     * @param rhs Output stream to insert the result representation into
     * @param lhs Result to insert
     * @return Reference to the modified output stream
     */
    friend std::ostream &operator<<(std::ostream &rhs, const Result &lhs)
    {
      return rhs << lhs.alpha << " " << lhs.psi << " " << lhs.d;
    }
  };

  // Simplified run interface bundling parameters and results
  void run(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
           SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
           const Parameters &params, Result &result, double &t1, double &t2)
  {
    run(ss, cv,
        // Unpack the parameter bundle
        params.interpolation, params.pixelScale, params.postScale,
        params.alphaEst, params.psiEst, params.searcher,
        params.compassAcceleration, params.psiFraction,
        params.quantileCompassFraction, params.partialSearch, params.alphaRad,
        params.psiRad, params.partialFullSearch, params.doubleSearch,
        params.fineSearch, params.doubleSearchFine, params.fineSearchSurround,
        // Results written to variables passed by reference
        result.alpha, result.psi, result.d, result.iAlpha, result.iPsi,
        result.iAlphaInv, result.iPsiInv, result.jAlphaOff, result.jPsiOff,
        result.fineMatch, t1, t2);
  }

  // Simplified interface omitting the time measurement
  void run(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
           SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
           const Parameters &params, Result &result)
  {
    // Dummies to write times into
    double t1 = 0, t2 = 0;
    // Delegate to the already simplified interface (hopefully inlined)
    run(ss, cv, params, result, t1, t2);
  }
};

} // namespace simd

#endif
