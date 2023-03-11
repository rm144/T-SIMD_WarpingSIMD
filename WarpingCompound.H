// ===========================================================================
//
// WarpingCompound.H --
// MinWarping / Warping interface
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

#ifndef WARPINGCOMPOUND_H_
#define WARPINGCOMPOUND_H_

#include "ComplexSearch.H"
#include "MinWarpingMatch.H"
#include "MinWarpingSearch.H"
#include "MinWarpingTemplate.H"
#include "SIMDException.H"
#include "SIMDExperimental.H"
#include "SIMDImage.H"
#include "WarpingParameter.H"
#include "WarpingSPS.H"
#include "WarpingSPSComputation.H"
#include "WarpingSearch.H"
#include "WarpingTemplate.H"

#include <cmath>
#include <vector>

namespace simd {

// TODO: - consistency check of parameters between SPS/match/template?

class WarpingSearcherSelector
{
public:
  enum SearcherType {
    minWarpingSearcher,
    warpingSearcher,
    minWarpingSearcher1SP,
    warpingSearcher1SP
  };
};

// =========================================================================
// WarpingCompound
// =========================================================================

// template class WarpingCompound provides minWarping and warping
// variables of type ComplexSearch which are initialized with shared
// scale-plane stack, shared match arrays, and shared parameters

template <typename SPSType, typename MatchType, typename CompassType,
          int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingCompound
  : public WarpingSearcherSelector
  , public ComplexSearchSelector
{
public:
  WarpingParameter<SPSType, SIMD_WIDTH> param;
  SymmetricScaleFactors symmScaleFac;
  SingleScaleFactor singleScaleFac;
  WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> *spsArray[NUM_SPS];
  WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> *spsArray1SP[NUM_SPS];
  MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN> *matchArray[NUM_MATCH];
  ComplexSearch<SPSType, MatchType, CompassType, SIMD_WIDTH, SIMD_ALIGN,
                MinWarpingTemplate>
    minWarping;
  ComplexSearch<SPSType, MatchType, CompassType, SIMD_WIDTH, SIMD_ALIGN,
                WarpingTemplate>
    warping;
  ComplexSearch<SPSType, MatchType, CompassType, SIMD_WIDTH, SIMD_ALIGN,
                MinWarpingTemplate>
    minWarping1SP;
  ComplexSearch<SPSType, MatchType, CompassType, SIMD_WIDTH, SIMD_ALIGN,
                WarpingTemplate>
    warping1SP;

  // TODO: Not copiable due to resource management unclear
  WarpingCompound(const WarpingCompound &) = delete;

  WarpingCompound(int width, int nAlpha, int nPsi, int nScalePlanes,
                  double maxScaleFactor, double maxThresholdMinWarping,
                  double rhoMinMinWarping, double rhoMaxMinWarping,
                  double rhoMinWarping, double rhoMaxWarping, int nRhoWarping)
    : param(width, nAlpha, nPsi), symmScaleFac(nScalePlanes, maxScaleFactor),
      singleScaleFac(),
      // search templates of all ComplexSearch members are deleted
      // in ~ComplexSearch
      minWarping(param, spsArray,
                 new MinWarpingTemplate(
                   width, param.wSPS, symmScaleFac.scaleFactors,
                   maxThresholdMinWarping, rhoMinMinWarping, rhoMaxMinWarping),
                 matchArray),
      warping(param, spsArray,
              new WarpingTemplate(width, param.wSPS, symmScaleFac.scaleFactors,
                                  rhoMinWarping, rhoMaxWarping, nRhoWarping),
              matchArray),
      minWarping1SP(param, spsArray1SP,
                    new MinWarpingTemplate(width, param.wSPS,
                                           singleScaleFac.scaleFactors,
                                           maxThresholdMinWarping,
                                           rhoMinMinWarping, rhoMaxMinWarping),
                    matchArray),
      warping1SP(param, spsArray1SP,
                 new WarpingTemplate(width, param.wSPS,
                                     singleScaleFac.scaleFactors, rhoMinWarping,
                                     rhoMaxWarping, nRhoWarping),
                 matchArray)
  {
    // TODO: are there any size/parameter checks needed in WarpingCompound?
    // TODO: - multi-in/multi-out cases in processing considered?
    // TODO: - checks from functions in SIMDImageFunctions.H now in
    // TODO:   WarpingSPSComputation::validForWidth
    // TODO: - which objects have to be created to do a full check?

    // create multipe scale-plane stacks and multiple match objects
    for (int i = 0; i < NUM_SPS; i++) {
      // printf("creating spsArray[%d]\n", i);
      spsArray[i] = new WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN>(
        param, symmScaleFac.scaleFactors.size());
      // printf("creating spsArrayMin[%d]\n", i);
      spsArray1SP[i] = new WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN>(
        param, singleScaleFac.scaleFactors.size());
    }
    for (int i = 0; i < NUM_MATCH; i++) {
      // printf("creating matchArray[%d]\n", i);
      matchArray[i] =
        new MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN>(nAlpha, nPsi);
      matchArray[i]->clear();
    }
  }

  virtual ~WarpingCompound()
  {
    for (int i = 0; i < NUM_SPS; i++) delete spsArray[i];
    for (int i = 0; i < NUM_SPS; i++) delete spsArray1SP[i];
    for (int i = 0; i < NUM_MATCH; i++) delete matchArray[i];
  }

  // =========================================================================
  // access to spsArray of different searchers
  // =========================================================================

  WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> **getSPSArray(int searcher)
  {
    switch (searcher) {
    case minWarpingSearcher: return minWarping.spsArray;
    case warpingSearcher: return warping.spsArray;
    case minWarpingSearcher1SP: return minWarping1SP.spsArray;
    case warpingSearcher1SP: return warping1SP.spsArray;
    default:
      throw SIMDException("WarpingCompound::getSPSArray",
                          "invalid warping searcher");
    }
  }

  // =========================================================================
  // SPS computation
  // =========================================================================

  // compute spsArray[SPS_ORIG]
  // (ImgType is last template argument since it can be deduced from ss/cv)
  // 27. Feb 18 (rm): vert.res. and hor. taken from images
  template <typename ProcType, typename MeasType, typename ImgType>
  void computeSPS(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    int interpolation, const std::vector<double> &pixelScale, double postScale,
    WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                          SIMD_ALIGN> &spsComp)
  {
    // printf("computeSPS(i=%d,pi=%g,po=%g\n",
    //        interpolation, pixelScale, postScale);
    spsComp.computeSPS(ss, cv, symmScaleFac.scaleFactors, interpolation,
                       pixelScale, postScale, *(spsArray[SPS_ORIG]));
  }

  // compute spsArray1SP[SPS_ORIG]
  // (ImgType is last template argument since it can be deduced from ss/cv)
  // 27. Feb 18 (rm): vert.res. and hor. taken from images
  template <typename ProcType, typename MeasType, typename ImgType>
  void computeSPS1SP(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    int interpolation, const std::vector<double> &pixelScale, double postScale,
    WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                          SIMD_ALIGN> &spsComp)
  {
    // printf("computeSPS1SP(i=%d,pi=%g,po=%g\n",
    //        interpolation, pixelScale, postScale);
    spsComp.computeSPS(ss, cv, symmScaleFac.scaleFactors, interpolation,
                       pixelScale, postScale, *(spsArray[SPS_ORIG]));
    spsArray[SPS_ORIG]->planeMinimum(*(spsArray1SP[SPS_ORIG]));
  }
};

} // namespace simd

#endif
