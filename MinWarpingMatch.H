// ===========================================================================
//
// MinWarpingMatch.H --
// C++ template for match array of minWarping algorithm
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

#ifndef _MINWARPING_MATCH_H_
#define _MINWARPING_MATCH_H_

#include "SIMDImage.H"

#include <limits>
// 20. Apr 18 (rm): thanks to Benedikt Volkmer
#include "SIMDTypes.H"
#include "SIMDVecAll.H"

namespace simd {

// ===========================================================================
// MinWarpingMatch
// ===========================================================================

template <typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingMatch
{
public:
  const int nAlpha, nPsi;
  const MatchType invalid;
  SIMDImage<MatchType, SIMD_WIDTH, SIMD_ALIGN> match;
  int *modulo_na, *modulo_na_pna;

  MinWarpingMatch(int nAlpha, int nPsi,
                  MatchType invalid = SIMDTypeInfo<MatchType>::max())
    : nAlpha(nAlpha), nPsi(nPsi), invalid(invalid), match(nAlpha, nPsi),
      modulo_na(0), modulo_na_pna(0)
  {
    // match = new SIMDImage<MatchType,SIMD_ALIGN>(nAlpha, nPsi);
    // modulo nAlpha
    // 19. Jun 13 (rm): new alpha handling, factor 3 instead of 2
    modulo_na = new int[3 * nAlpha];
    for (int i = 0; i < 3 * nAlpha; i++) modulo_na[i] = i % nAlpha;
    modulo_na_pna = modulo_na + nAlpha;
  }

  virtual ~MinWarpingMatch()
  {
    if (modulo_na) delete[] modulo_na;
  }

  MatchType *operator[](int iPsi)
  {
    if ((iPsi < 0) || (iPsi >= nPsi))
      throw SIMDException("MinWarpingMatch::operator[]", "invalid index");
    return match[iPsi];
  }

  // check whether size and invalid are the same
  bool sameConfig(const MinWarpingMatch &otherMatch)
  {
    return (nAlpha == otherMatch.nAlpha) && (nPsi == otherMatch.nPsi) &&
           (invalid == otherMatch.invalid);
  }

  bool save(FILE *f) const { return match.save(f); }

  void clear() { match.clear(); }

  void fillInvalid()
  {
    // if MatchType is an unsigned integer and invalid coincides with
    // the maximum value of MatchType, we can fill with 0xff bytes,
    // otherwise we fill in a loop
    if (SIMDTypeInfo<MatchType>::isInteger &&
        !SIMDTypeInfo<MatchType>::isSigned &&
        (invalid == SIMDTypeInfo<MatchType>::max()))
      match.fillByte(0xff);
    else
      match.fill(invalid);
  }

  void invertSolution(int iAlpha0, int iPsi0, int &iAlpha1, int &iPsi1)
  {
    iPsi1   = (nPsi - iPsi0) % nPsi;
    iAlpha1 = modulo_na_pna[nAlpha / 2 - (iPsi0 * nAlpha) / nPsi + iAlpha0];
  }

  // TODO: can we accelerate invertPsiRangeFrom using vector instructions?
  void invertPsiRangeFrom(const MinWarpingMatch &otherMatch,
                          const int iPsi0Start, const int iPsi0End)
  {
    if (!sameConfig(otherMatch))
      throw SIMDException("MinWarpingMatch::invertPsiRangeFrom",
                          "configuration mismatch");
    int nAlpha2 = nAlpha / 2;
    // indices with suffix 0 for match in normal order (SS-CV match)
    // indices with suffix 1 for match in inverse order (CV-SS match)
    for (int iPsi0 = iPsi0Start; iPsi0 <= iPsi0End; iPsi0++) {
      // psi1 = -psi0
      int iPsi1        = (nPsi - iPsi0) % nPsi;
      int iPsi1_nAlpha = iPsi1 * nAlpha;
      int iPsi0_nAlpha = iPsi0 * nAlpha;
      // 29. Aug 15 (rm): iPsi0*nAlpha->iPsi0_nAlpha
      int nAlpha2_iPsi0_nAlpha_nPsi = nAlpha2 - iPsi0_nAlpha / nPsi;
      for (int iAlpha0 = 0; iAlpha0 < nAlpha; iAlpha0++) {
        // alpha1 = pi - psi0 * (...) + alpha0, nAlpha added because of %
        // TODO: if nAlpha!=nPsi, can we have holes in the transformation?
        // int iAlpha1 = (nAlpha + nAlpha/2 -(iPsi0*nAlpha)/nPsi + iAlpha0)
        //               % nAlpha;
        // 22. Dec 13 (rm)
        int iAlpha1 = modulo_na_pna[nAlpha2_iPsi0_nAlpha_nPsi + iAlpha0];
        // int idx0 = iPsi0 * nAlpha + iAlpha0;
        // int idx1 = iPsi1 * nAlpha + iAlpha1;
        int idx0         = iPsi0_nAlpha + iAlpha0;
        int idx1         = iPsi1_nAlpha + iAlpha1;
        match.data[idx0] = otherMatch.match.data[idx1];
      }
    }
  }

  void invertFrom(const MinWarpingMatch &otherMatch)
  {
    invertPsiRangeFrom(otherMatch, 0, nPsi - 1);
  }

  void invertSinglePsiFrom(const MinWarpingMatch &otherMatch, int iPsiSingle)
  {
    match.fill(invalid);
    invertPsiRangeFrom(otherMatch, iPsiSingle, iPsiSingle);
  }

  void averageOf(const MinWarpingMatch &match1, const MinWarpingMatch &match2)
  {
    if (!(sameConfig(match1) && sameConfig(match2)))
      throw SIMDException("MinWarpingMatch::averageOf",
                          "configuration mismatch");
    SIMDVec<MatchType, SIMD_WIDTH> vec1, vec2, result, resultIsInvalid,
      invalidVec = set1<MatchType, SIMD_WIDTH>(invalid);
    for (size_t i = 0; i < match.size;
         i += SIMDVec<MatchType, SIMD_WIDTH>::elements) {
      vec1            = load<SIMD_WIDTH>(match1.match.data + i);
      vec2            = load<SIMD_WIDTH>(match2.match.data + i);
      resultIsInvalid = or_(cmpeq(vec1, invalidVec), cmpeq(vec2, invalidVec));
      result          = ifelse(resultIsInvalid, invalidVec, avgrd(vec1, vec2));
      store(match.data + i, result);
    }
  }

  void averageSinglePsiOf(const MinWarpingMatch &match1,
                          const MinWarpingMatch &match2, int iPsiSingle)
  {
    if (!(sameConfig(match1) && sameConfig(match2)))
      throw SIMDException("MinWarpingMatch::averageSinglePsiOf",
                          "configuration mismatch");
    SIMDVec<MatchType, SIMD_WIDTH> vec1, vec2, result, resultIsInvalid,
      invalidVec    = set1<MatchType, SIMD_WIDTH>(invalid);
    MatchType *row  = match[iPsiSingle];
    MatchType *row1 = match1.match[iPsiSingle];
    MatchType *row2 = match2.match[iPsiSingle];
    match.fill(invalid);
    for (int i = 0; i < nAlpha; i += SIMDVec<MatchType, SIMD_WIDTH>::elements) {
      vec1            = load<SIMD_WIDTH>(row1 + i);
      vec2            = load<SIMD_WIDTH>(row2 + i);
      resultIsInvalid = or_(cmpeq(vec1, invalidVec), cmpeq(vec2, invalidVec));
      result          = ifelse(resultIsInvalid, invalidVec, avgrd(vec1, vec2));
      store(row + i, result);
    }
  }

  void copyFrom(const MinWarpingMatch &otherMatch)
  {
    if (!sameConfig(otherMatch))
      throw SIMDException("MinWarpingMatch::copyFrom", "size mismatch");
    match.copyFrom(otherMatch);
  }

  // looking for best (minimal) match
  void bestMatchFull(int &iAlphaMin, int &iPsiMin, MatchType &dMin)
  {
    MatchType *minp = match.data;
    for (MatchType *mp = match.data + 1; mp < match.data + nAlpha * nPsi; mp++)
      if (*mp < *minp) minp = mp;
    iPsiMin   = (minp - match.data) / nAlpha;
    iAlphaMin = (minp - match.data) % nAlpha;
    // printf("minp = %p, match = %p, minp - match = %d, iPsiMin = %d,
    //         iAlphaMin = %d, nPsi = %d, nAlpha = %d\n",
    //	    minp, match, minp - match, iPsiMin, iAlphaMin, nPsi, nAlpha);
    dMin = *minp;
  }

  // obtain min and max which is not invalid (for visualization, not
  // speed-tuned)
  // fails if no pixel is valid (all pixels are invalid):
  // minv, maxv not changed
  void extremaExceptInvalid(MatchType &min, MatchType &max)
  {
    MatchType minv  = SIMDTypeInfo<MatchType>::min(),
              maxv  = SIMDTypeInfo<MatchType>::max();
    MatchType *minp = &maxv, *maxp = &minv;
    for (MatchType *mp = match.data; mp < match.data + nAlpha * nPsi; mp++)
      if (*mp != invalid) {
        if (*mp < *minp) minp = mp;
        if (*mp > *maxp) maxp = mp;
      }
    min = *minp;
    max = *maxp;
  }

  template <typename ImgType>
  void getImage(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &res,
                double shift = 0.0, double scale = 1.0) const
  {
    ImgType imgMax = SIMDTypeInfo<ImgType>::max();
    copyShiftAndScale(match, shift, scale, res);
    for (size_t i = 0; i < match.size; i++)
      if (match.data[i] == invalid) res.data[i] = imgMax;
  }
};

} // namespace simd

#endif
