// ===========================================================================
//
// MinWarpingSearchRange.H --
// class defining search range for partial minWarping
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

#ifndef _MINWARPING_SEARCHRANGE_H_
#define _MINWARPING_SEARCHRANGE_H_

#include "SIMDImage.H"

#include <algorithm>
#include <cmath>
#include <vector>
// 20. Apr 18 (rm): thanks to Benedikt Volkmer
#include "SIMDVecAll.H"

namespace ns_simd {

// ===========================================================================
// MinWarping search range
// ===========================================================================

template <typename SPSType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingSearchRange
{
public:
  // the number of SPS (scale-plane stack) elements processed
  // in parallel in a SIMD vector
  const int simd_sps_elems;
  // full search range
  const int nAlpha, nPsi;
  // restricted search range: indices
  std::vector<int> alphaIndices, psiIndices;
  // these are currently unused, but are computed
  int iAlphaSIMDBlockStart, iAlphaSIMDBlockStop;

  MinWarpingSearchRange(int nAlpha, int nPsi)
    : // simd_sps_elems(SIMD_WIDTH / sizeof(SPSType)),
      simd_sps_elems(NUM_SIMDVEC_ELEMENTS(SPSType, SIMD_WIDTH)), nAlpha(nAlpha),
      nPsi(nPsi)
  {}

  void setPartialRange(double alphaEst, double alphaRadius, double psiEst,
                       double psiRadius)
  {
    getRangeIndicesSpacingSIMDBlock(alphaEst, alphaRadius, nAlpha, alphaIndices,
                                    iAlphaSIMDBlockStart, iAlphaSIMDBlockStop);
    getRangeIndices(psiEst, psiRadius, nPsi, psiIndices);
  }

  void setFullRange()
  {
    getAllIndicesSpacingSIMDBlock(nAlpha, alphaIndices, iAlphaSIMDBlockStart,
                                  iAlphaSIMDBlockStop);
    getAllIndices(nPsi, psiIndices);
  }

  void setSingleRange(int iAlpha, int iPsi)
  {
    getRangeIndicesSpacingSIMDBlockSingle(
      iAlpha, alphaIndices, iAlphaSIMDBlockStart, iAlphaSIMDBlockStop);
    getRangeIndicesSingle(iPsi, psiIndices);
  }

  template <typename CompassType>
  void setFullCompassRange(
    double psiFraction,
    const SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> &compass)
  {
    getAllIndicesSpacingSIMDBlock(nAlpha, alphaIndices, iAlphaSIMDBlockStart,
                                  iAlphaSIMDBlockStop);
    getCompassIndices(psiFraction, compass, psiIndices);
  }

protected:
  // determine index of an angle
  // same function as angle2index in MathExt.H/Math.C
  // moved here to make it independent of other parts of PROG system
  static int a2i(double a, int n)
  {
    // map angle to interval [0,2pi) and return mapped angle / 2pi
    double ra   = a * (0.5 / M_PI);
    double ra01 = ra - ::floor(ra);
    // determine corresponding index
    int i = (int) (ra01 * n + 0.5);
    return (i == n) ? 0 : i;
  }

  // provides a vector with indices given an angular range
  // [value-radius,value+radius], n specifies the index range [0,n-1]
  // only works for radius in [0,pi) (and may even fail in that case if
  // n is very small)
  static void getRangeIndices(double value, double radius, int n,
                              std::vector<int> &indices)
  {
    int il, ir, i, I;

    // printf("center %u\n", a2i(value, n));
    // index of left border
    il = a2i(value - radius, n);
    // index of right border
    ir = a2i(value + radius, n);
    // printf("il, ir = %u, %u\n", il, ir);
    if (il <= ir) {
      // left border is left of right border
      indices.resize((ir - il) + 1);
      for (i = il, I = 0; i <= ir; i++, I++) indices[I] = i;
    } else {
      // left border is right of right border
      indices.resize((ir + 1) + (n - il));
      for (i = 0, I = 0; i <= ir; i++, I++) indices[I] = i;
      for (i = il; i < n; i++, I++) indices[I] = i;
    }
  }

  static void getRangeIndicesSingle(int index, std::vector<int> &indices)
  {
    indices.resize(1);
    indices[0] = index;
  }

  // vector with all indices
  static void getAllIndices(int n, std::vector<int> &indices)
  {
    indices.resize(n);
    for (int i = 0; i < n; i++) indices[i] = i;
  }

  // indices obtained from a compass estimate ("fraction" smallest values)
  template <typename CompassType>
  static void getCompassIndices(
    double fraction,
    const SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> &compass,
    std::vector<int> &indices)
  {
    int n        = compass.w;
    int iLargest = int(fraction * (n - 1));
    CompassType largest;
    // we need to make a copy since nth_element overwrites the array
    // CompassType compassCopy[n];
    SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> compassCopy = compass;
    // memcpy(compassCopy, compass.data, n * sizeof(CompassType));

    std::nth_element(compassCopy.data, compassCopy.data + iLargest,
                     compassCopy.data + n);
    indices.clear();
    largest = compassCopy.data[iLargest];
    // printf("largest %u @ %u\n", largest, iLargest);
    for (int i = 0; i < n; i++)
      // here we work on the original compass array
      if (compass.data[i] <= largest) {
        // printf("sel %d\n", i);
        indices.push_back(i);
      }
  }

  // (for 16: was (N >> 4) + 1)
  // TODO: fast way to find exponent of a power of 2, and so use >> again?
  // 10. Oct 14 (rm): left this unsigned to possibly avoid case distinctions
  unsigned sizeBlk(unsigned N) { return ((N / simd_sps_elems) + 1); };
  // (for 16: was I & (~15))
  unsigned lowerBlk(unsigned I) { return (I & ~(simd_sps_elems - 1)); };

  // for vector indices with spacing of simd_sps_elems:
  // provides a vector with indices given an angular range
  // [value-radius,value+radius], n specifies the index range [0,n-1]
  // only works for radius in [0,pi) (and may even fail in that case if
  // n is very small)
  void getRangeIndicesSpacingSIMDBlock(double value, double radius, int n,
                                       std::vector<int> &indices,
                                       int &iStartBlk, int &iStopBlk)
  {
    int il, ir, ilBlk, irBlk, i, I;

    // printf("center %u\n", a2i(value, n));
    // index of left border
    il = a2i(value - radius, n);
    // index of right border
    ir = a2i(value + radius, n);
    // with spacing of simd_sps_elems we need the start index of the
    // block of simd_sps_elems which contains il or ir
    iStartBlk = ilBlk = lowerBlk(il);
    irBlk             = lowerBlk(ir);
    // 19. Jun 13 (rm): no wrapping here, and assert iStartBlk < iStopBlk
    // iStopBlk = (irBlk + simd_sps_elems) % n;
    iStopBlk = irBlk + simd_sps_elems;
    if (iStopBlk <= iStartBlk) iStopBlk += n;
    // printf("il, ir = %u, %u\n", il, ir);
    if (il <= ir) {
      // left border is left of right border
      indices.resize(sizeBlk(irBlk - ilBlk));
      for (i = ilBlk, I = 0; i <= irBlk; i += simd_sps_elems, I++)
        indices[I] = i;
    } else {
      // left border is right of right border
      indices.resize(sizeBlk(irBlk) + sizeBlk((n - ilBlk) - 1));
      for (i = 0, I = 0; i <= irBlk; i += simd_sps_elems, I++) indices[I] = i;
      for (i = ilBlk; i < n; i += simd_sps_elems, I++) indices[I] = i;
    }
  }

  void getRangeIndicesSpacingSIMDBlockSingle(int index,
                                             std::vector<int> &indices,
                                             int &iStartBlk, int &iStopBlk)
  {
    indices.resize(1);
    iStartBlk = indices[0] = lowerBlk(index);
    iStopBlk               = iStartBlk + simd_sps_elems;
  }

  // for tests: vector with all indices with spacing of 8
  void getAllIndicesSpacingSIMDBlock(int n, std::vector<int> &indices,
                                     int &iStartBlk, int &iStopBlk)
  {
    // identical to (n-1) / simd_sps_elems + 1:
    indices.resize(sizeBlk(n - 1));
    int j = 0;
    for (int i = 0; i < n; i += simd_sps_elems, j++) indices[j] = i;
    iStartBlk = 0;
    iStopBlk  = n;
  }
};

} // namespace ns_simd

#endif
