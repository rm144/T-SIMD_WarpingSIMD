// ===========================================================================
//
// MinWarpingSearchCollection.H --
// collection of search algorithms for minWarping
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
#ifndef MINWARPING_SEARCH_COLLECTION_H_
#define MINWARPING_SEARCH_COLLECTION_H_

#include "MinWarpingSearch.H"
#include "SIMDException.H"
#include "SIMDTypes.H"
#include "SIMDVecAll.H"

#include <limits>
#include <vector>

namespace simd {

// TODO: test if _mm_lddqu_si128 is faster than loadu
// TODO: faster with different arrangement of SPS (more cache-efficient?)
// TODO: test gcc 4.9 runtime (M. Pflanzer: 50% slower!?)

// ###########################################################################
// full matching
// ###########################################################################

// ===========================================================================
// XPAY
// ===========================================================================

// XPAY denotes the search order: x, psi, alpha, y
template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingFull_XPAY
  : public MinWarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
{
  virtual void search(const MinWarpingTemplate &minTemplate,
                      const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
                      MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN> &match)
  {
    // we assume that sizeof(SPSType) <= sizeof(MatchType)
    // (since we add up SPSType elements in MatchType sum)
    // (just to provide the user with a clear compiler message)
    SIMD_STATIC_ASSERT(sizeof(SPSType) <= sizeof(MatchType));

    // 30. Sep 14 (rm): local copies (somewhat faster, surprisingly)
    //  1. Jun 15 (rm): now const
    const int w    = minTemplate.w;
    const int nPsi = sps.param.nPsi, nAlpha = sps.param.nAlpha;
    int *const shuffleModulo = sps.param.shuffleModulo;
    int *const modulo2wpw    = sps.param.modulo2wpw;
    int *const jPsiVec       = sps.param.jPsiVec;
    int *const jAlphaVec     = sps.param.jAlphaVec;
    MatchType *const _match  = match.match.data;
    SPSType *const _stack    = sps.stack.data;
    const int simd_sps_elems = sps.param.simd_sps_elems;
    // numInVecs is always 1, see assertion above
    const int numOutVecs = numOutputSIMDVecs<MatchType, SPSType>();
    const int numOutVecElems =
      // numOutputSIMDVecElements<MatchType,SPSType,SIMD_WIDTH>();
      numSIMDVecElements<MatchType, SIMD_WIDTH>();
    // numElems is always simd_sps_elems, see assertion above
    const SPSType maxSPS = SIMDTypeInfo<SPSType>::max();
    // SIMDVec<SPSType,SIMD_WIDTH> zero = setzero<SPSType,SIMD_WIDTH>();
    // faster if not defined inside loops?
    SIMDVec<SPSType, SIMD_WIDTH> minv, v;
    SIMDVec<MatchType, SIMD_WIDTH> minvM[numOutVecs];

    // we write simd_sps_elems results to match array in one go
    if (nAlpha % simd_sps_elems != 0)
      throw SIMDException("MinWarpingFull_XPAY::search",
                          "nAlpha has to be a multiple of simd_sps_elems");

    // here the match array is accumulated
    // TODO: entire clear in the beginning or split (when needed) - which
    // TODO: is faster?
    match.clear();
    // go trough all template columns
#if 0
      for (int jX = 0; jX < w; jX++) {
#else
    // slightly faster, but only on AVX
    // idea: SPS is accessed in stepAlpha blocks
    const int stepAlpha = sps.param.stepAlpha;
    for (int jXbase = 0; jXbase < stepAlpha; jXbase++)
      for (int jX = jXbase; jX < w; jX += stepAlpha) {
#endif
    // get the template column
    lp_min_templt_col_t tCol = minTemplate.minTemplateColumns[jX];
    // we only do min search for non-empty template columns
    // (otherwise we would obtain initial value of min!)
    if (tCol->len != 0) {
      // precompute part of shuffled (jX + jAlpha) % w
      int *shuffleModuloPlusX = shuffleModulo + jX;
      // pointer into match array (arranged in nPsi rows @ nAlpha cols)
      MatchType *matchp = _match;
      // go through all psi parameters
      for (int iPsi = 0; iPsi < nPsi; iPsi++) {
        // precompute part of (jY - jPsi) % w (times wSPS)
        int *modulo2wpwMinusPsi = modulo2wpw - jPsiVec[iPsi];
        // go through all alpha parameters
        for (int iAlpha = 0; iAlpha < nAlpha; iAlpha += simd_sps_elems) {
          // precompute scale-plane stack with column offset
          // select word in scale-plane row, shuffled:
          // jTheta = (jX + jAlpha) % w
          // 30. Sep 14 (rm): array access instead of dereferentiation
          SPSType *spsCol = _stack + shuffleModuloPlusX[jAlphaVec[iAlpha]];
          // we search for a minimum, so we start with the highest
          // positive value
          minv = set1<SPSType, SIMD_WIDTH>(maxSPS);
          // go through all elements in this template column
          // 19. Jun 15 (rm): pointer-based loop is faster
          const min_templt_col_entry_t *const start = tCol->entries();
          const min_templt_col_entry_t *const end   = start + tCol->len;
          for (const min_templt_col_entry_t *entry = start; entry < end;
               entry++) {
            // compute address of sps word
            // 30. Sep 14 (rm): array access instead of dereferentiation
            SPSType *vp =
              // select scale plane
              spsCol + entry->sp +
              // select row in scale plane: jDelta=wSPS * [(jY - jPsi) % w]
              modulo2wpwMinusPsi[entry->iy];
            // we need to have an unaligned access
            v = loadu<SIMD_WIDTH>(vp);
            // update minimum (vpminub)
            minv = min(minv, v);
          }
          // accumulate minima
          extend(minv, minvM);
          for (int vo = 0; vo < numOutVecs; vo++, matchp += numOutVecElems)
            store(matchp, adds(load<SIMD_WIDTH>(matchp), minvM[vo]));
        }
      }
    }
  }
}
};

// ###########################################################################
// min-warping curve (structure taken from full-XPAY, but no SIMD code)
// ###########################################################################

// 23. Feb 18 (rm): added jYVec
template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
void minWarpingCurve(const MinWarpingTemplate &minTemplate,
                     const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
                     int iAlpha, int iPsi,
                     std::vector<int> &jDeltaVec,     // over jTheta
                     std::vector<int> &planeIndexVec, // over jTheta
                     std::vector<SPSType> &minVec,    // over jTheta
                     std::vector<int> &jYVec,         // over jX
                     std::vector<int> &jXVec,         // over jTheta
                     std::vector<int> &jThetaVec,     // over jX
                     MatchType &matchSum)
{
  // 30. Sep 14 (rm): local copies (somewhat faster, surprisingly)
  //  1. Jun 15 (rm): now const
  const int w = minTemplate.w, wSPS = minTemplate.wSPS;
  int *const shuffle    = sps.param.shuffle;
  int *const modulo     = sps.param.modulo;
  int *const modulopw   = sps.param.modulopw;
  SPSType *const _stack = sps.stack.data;
  const int numPlanes   = sps.stack.numPlanes;
  const SPSType maxSPS  = SIMDTypeInfo<SPSType>::max();
  // resize result vectors, gaps in min-template stay negative!
  jDeltaVec.resize(w, -1);
  planeIndexVec.resize(w, -1);
  minVec.resize(w, SIMDTypeInfo<SPSType>::max());
  jYVec.resize(w, -1);
  jXVec.resize(w, -1);
  jThetaVec.resize(w, -1);
  matchSum = 0;
  // jAlpha, jPsi:
  const int jAlpha = sps.param.jAlphaVec[iAlpha];
  const int jPsi   = sps.param.jPsiVec[iPsi];
  // go trough all template columns
  for (int jX = 0; jX < w; jX++) {
    // get the template column
    lp_min_templt_col_t tCol = minTemplate.minTemplateColumns[jX];
    // we only do min search for non-empty template columns
    // (otherwise we would obtain initial value of min!)
    if (tCol->len != 0) {
      // precompute scale-plane stack with column offset
      // select word in scale-plane row, shuffled:
      // jTheta = (jX + jAlpha) % w
      int jTheta = modulo[jX + jAlpha];
      // 30. Sep 14 (rm): array access instead of dereferentiation
      SPSType *spsCol = _stack + shuffle[jTheta];
      // we search for a minimum, so we start with the highest
      // positive value
      SPSType minv = maxSPS;
      // initialize with invalid values, should be updated in any case
      int jDeltaMin = w, planeIndexMin = numPlanes, jYMin = w;
      // go through all elements in this template column
      // 19. Jun 15 (rm): pointer-based loop is faster
      const min_templt_col_entry_t *const start = tCol->entries();
      const min_templt_col_entry_t *const end   = start + tCol->len;
      for (const min_templt_col_entry_t *entry = start; entry < end; entry++) {
        int jY = (int) entry->iy;
        // jDelta = (jY - jPsi) % w
        int jDelta = modulopw[jY - jPsi];
        // entry->sp is an offset, not an index, convert back to index
        int planeIndex = entry->sp / (wSPS * w);
        // compute address of sps word
        // 30. Sep 14 (rm): array access instead of dereferentiation
        SPSType *vp =
          // select scale plane
          spsCol + entry->sp +
          // select row in scale plane: jDelta = wSPS * [(jY - jPsi) % w]
          wSPS * jDelta;
        // we need to have an unaligned access
        SPSType v = *vp;
        // update minimum
        // 6. Oct 15 (rm): <= instead of < (if all values are equal,
        //                 jDeltaMin was not updated)
        if (v <= minv) {
          minv          = v;
          jDeltaMin     = jDelta;
          planeIndexMin = planeIndex;
          jYMin         = jY;
        }
      }
      // these vectors run over jTheta
      jDeltaVec[jTheta]     = jDeltaMin;
      planeIndexVec[jTheta] = planeIndexMin;
      minVec[jTheta]        = minv;
      // jYVec runs over jX!
      jYVec[jX] = jYMin;
      //  4. Jul 18 (rm): association jX <-> jTheta
      jXVec[jTheta] = jX;
      jThetaVec[jX] = jTheta;
      matchSum += minv;
    }
  }
}

// ###########################################################################
// partial matching
// ###########################################################################

// TODO: check: is partial search really working on arbitrary range sets
// TODO: with holes?

// ===========================================================================
// XPAY
// ===========================================================================

// XPAY denotes the search order: x, psi, alpha, y
template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingPartial_XPAY
  : public MinWarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
{
public:
  virtual void search(
    const MinWarpingTemplate &minTemplate,
    const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
    const MinWarpingSearchRange<SPSType, SIMD_WIDTH, SIMD_ALIGN> &range,
    MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN> &match, int jAlphaOff,
    int jPsiOff)
  {
    // we assume that sizeof(SPSType) <= sizeof(MatchType)
    // (since we add up SPSType elements in MatchType sum)
    // (just to provide the user with a clear compiler message)
    SIMD_STATIC_ASSERT(sizeof(SPSType) <= sizeof(MatchType));

    //  1. Oct 14 (rm): local copies (sometimes faster, not here though)
    //  3. Jun 15 (rm): now const
    const int w                = minTemplate.w;
    const int nAlpha           = sps.param.nAlpha;
    int *const shuffleModuloPW = sps.param.shuffleModuloPW;
    int *const modulo2wpw      = sps.param.modulo2wpw;
    int *const jAlphaVec       = sps.param.jAlphaVec;
    int *const jPsiVec         = sps.param.jPsiVec;
    MatchType *const _match    = match.match.data;
    SPSType *const _stack      = sps.stack.data;
    const int simd_sps_elems   = sps.param.simd_sps_elems;
    // numInVecs is always 1, see assertion above
    const int numOutVecs = numOutputSIMDVecs<MatchType, SPSType>();
    const int numOutVecElems =
      // numOutputSIMDVecElements<MatchType,SPSType,SIMD_WIDTH>();
      numSIMDVecElements<MatchType, SIMD_WIDTH>();
    // numElems is always simd_sps_elems, see assertion above
    const SPSType maxSPS = SIMDTypeInfo<SPSType>::max();
    // faster if not defined inside loops
    SIMDVec<MatchType, SIMD_WIDTH> zero = setzero<MatchType, SIMD_WIDTH>();
    SIMDVec<SPSType, SIMD_WIDTH> minv, v;
    SIMDVec<MatchType, SIMD_WIDTH> minvM[numOutVecs];

    // we write simd_sps_elems results to match array in one go
    if (nAlpha % simd_sps_elems != 0)
      throw SIMDException("MinWarpingPartial_XPAY::search",
                          "nAlpha has to be a multiple of simd_sps_elems");

    // this is a bit more complicated than in the PAXY version: since we
    // accumulate in the match array, we need to initialize all match
    // elements which will not be touched later with invalid (e.g. 0xffff)
    // and all elements which will be touched with 0x0000,
    // all untouched elements in alpha-psi array (match) are invalid
    // (e.g. 0xffff)
    match.fillInvalid();
    // the elements which will be touched are zeroed out
    // go through all psi parameters
    for (unsigned int kPsi = 0; kPsi < range.psiIndices.size(); kPsi++) {
      int iPsi            = range.psiIndices[kPsi];
      MatchType *matchRow = _match + iPsi * nAlpha;
      // go through all alpha parameters
      // note that alphaIndices contain indices with a spacing of
      // simd_sps_elems!
      for (unsigned int kAlpha = 0; kAlpha < range.alphaIndices.size();
           kAlpha++) {
        int iAlpha = range.alphaIndices[kAlpha];
        // set consecutive values to 0
        store(matchRow + iAlpha, zero, numOutVecs);
      }
    }
    // go trough all template columns
    // TODO: same x-loop spliting as in full search?
    for (int jX = 0; jX < w; jX++) {
      // get the template column
      lp_min_templt_col_t tCol = minTemplate.minTemplateColumns[jX];
      // we only do min search for non-empty template columns
      // (otherwise we would obtain initial value of min!)
      if (tCol->len != 0) {
        // precompute part of shuffled (jX + jAlpha) % w
        // 20. Nov 15 (rm): shuffleModuloPW instead of shuffleModulo
        // 17. Dec 15 (rm): addition of offset moved here
        int *shuffleModuloPlusXPlusAlphaOff = shuffleModuloPW + jX + jAlphaOff;
        // go through all psi parameters
        for (unsigned int kPsi = 0; kPsi < range.psiIndices.size(); kPsi++) {
          // get index
          int iPsi = range.psiIndices[kPsi];
          // precompute part of (jY - jPsi) % w (times wSPS)
          // 20. Nov 15 (rm): offset added
          int *modulo2wpwMinusPsi = modulo2wpw - (jPsiVec[iPsi] + jPsiOff);
          // pointer into match row
          MatchType *matchRow = _match + iPsi * nAlpha;
          // go through all alpha parameters
          // note that alphaIndices contain indices with a spacing of 8!
          for (unsigned int kAlpha = 0; kAlpha < range.alphaIndices.size();
               kAlpha++) {
            // get index
            int iAlpha = range.alphaIndices[kAlpha];
            // precompute scale-plane stack with column offset
            // select word in scale-plane row, shuffled: (jX + jAlpha) % w
            //  1. Oct 14 (rm): index access
            // 20. Nov 15 (rm): offset added
            SPSType *spsCol =
              _stack + shuffleModuloPlusXPlusAlphaOff[jAlphaVec[iAlpha]];
            // we search for a minimum, so we start with the highest
            // positive value
            minv = set1<SPSType, SIMD_WIDTH>(maxSPS);
            // go through all elements in this template column
            // 19. Jun 15 (rm): pointer-based loop is faster
            const min_templt_col_entry_t *const start = tCol->entries();
            const min_templt_col_entry_t *const end   = start + tCol->len;
            for (const min_templt_col_entry_t *entry = start; entry < end;
                 entry++) {
              // compute address of sps word
              SPSType *vp =
                // select scale plane
                //  1. Oct 14 (rm): index access
                spsCol + entry->sp +
                // select row in scale plane: wSPS * [(jY - jPsi) % w]
                modulo2wpwMinusPsi[entry->iy];
              // we need to have an unaligned access
              v = loadu<SIMD_WIDTH>(vp);
              // update minimum
              minv = min(minv, v);
            }
            // accumulate minima
            // MatchType *matchp = match + iPsi * nAlpha + iAlpha;
            MatchType *matchp = matchRow + iAlpha;
            extend(minv, minvM);
            for (int vo = 0; vo < numOutVecs; vo++, matchp += numOutVecElems)
              store(matchp, adds(load<SIMD_WIDTH>(matchp), minvM[vo]));
          }
        }
      }
    }
  }
};

} // namespace simd

#endif
