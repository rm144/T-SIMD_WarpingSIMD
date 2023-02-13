// ===========================================================================
//
// SIMDImageFunctions.H --
// functions operating on SIMDImage objects
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

#ifndef _SIMD_IMAGE_FUNCTIONS_H_
#define _SIMD_IMAGE_FUNCTIONS_H_

#include "SIMDImage.H"
#include "SIMDVecAll.H"
#include "saturate_cast.H"

#include <algorithm>
#include <math.h>
#include <vector>

/**
 * @brief Interpolation method selectors
 */
enum class Interpolation : int {
  /**
   * @brief Nearest neighbor interpolation method selector
   */
  NEAREST = 0,

  /**
   * @brief (Bi-)Linear interpolation method selector
   */
  LINEAR = 1
};

namespace ns_simd {

// ===========================================================================
// NOTES:
// - copy*() functions allow same image as input and output
// ===========================================================================

// ===========================================================================
// functions operating on SIMDImage
// ===========================================================================

// ---------------------------------------------------------------------------
// sameness / overlap check
// ---------------------------------------------------------------------------

// images with same type could be identical (same data pointer)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
bool same(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &i1,
          const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &i2)
{
  return i1.same(i2);
}

// images with different type can never be identical
template <typename T, typename U, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
bool same(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &,
          const SIMDImage<U, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &)
{
  return false;
}

// stacks with same type could be identical (same data pointer)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
bool same(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &s1,
          const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &s2)
{
  return s1.same(s2);
}

// stacks with different type can never be identical
template <typename T, typename U, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
bool same(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &,
          const SIMDStack<U, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &)
{
  return false;
}

// sharedptrs with same element type can overlap (some data pointer)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
bool overlap(const SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN> &i1,
             const SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN> &i2)
{
  // check if start of i1 falls into range of i2 and vice versa
  return (((i1.data >= i2.data) && (i1.data < i2.data + i2.size)) ||
          ((i2.data >= i1.data) && (i2.data < i1.data + i1.size)));
}

// sharedptrs with different element types can never overlap
template <typename T, typename U, int SIMD_WIDTH, int SIMD_ALIGN>
bool overlap(const SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN> &,
             const SharedSIMDPtr<U, SIMD_WIDTH, SIMD_ALIGN> &)
{
  return false;
}

// ---------------------------------------------------------------------------
// sum of absolute values
// ---------------------------------------------------------------------------

// sum of absolute values in all columns
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void sumOfAbs(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img,
              SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &absSum)
{
  SIMDVec<T, SIMD_WIDTH> as;
  if (overlap(img, absSum))
    throw SIMDException("sumOfAbs", "input and output image overlap");
  // number of elements in vector
  const int simd_elems = as.elements;
  const int w = img.w, h = img.h;
  absSum.resize(w, 1, img.addOn);
  // pointer to input image: we work with two pointers, one (pi0)
  // pointing to the top row, the other (pi1) traversing the column block
  // pi0 traverses top row in jumps of simd_elems values
  const T *pi0 = img.data;
  // pas traverses result (single row) in jumps of simd_elems values
  T *pas = absSum.data;
  // go through top column in jumps of simd_elems
  for (int x = 0; x < w; x += simd_elems) {
    // pi1 traverses block of 8 columns, starting at the top row
    const T *pi1 = pi0;
    // sum of absolute values = 0
    // SIMDVec<T,SIMD_WIDTH>
    as = setzero<T, SIMD_WIDTH>();
    // go through column block
    for (int y = 0; y < h; y++) {
      // add absolute values with saturation
      // 25. Jan 14 (rm): using adds: we have to avoid overflow!
      as = adds(as, abs(load<SIMD_WIDTH>(pi1)));
      // go to next row
      pi1 += w;
    }
    // store sum
    store(pas, as);
    // next block of 8 columns
    pi0 += simd_elems;
    pas += simd_elems;
  }
}

// ---------------------------------------------------------------------------
// match operators for numerator
// ---------------------------------------------------------------------------

// match operator for NSAD
template <typename T, int SIMD_WIDTH>
struct MatchAbsDiff
{
  static SIMD_INLINE SIMDVec<T, SIMD_WIDTH> match(
    const SIMDVec<T, SIMD_WIDTH> &a, const SIMDVec<T, SIMD_WIDTH> &b)
  {
    return abs(sub(a, b));
  }
};

// match operator for ASC2
template <typename T, int SIMD_WIDTH>
struct MatchAsc2
{
  static SIMD_INLINE SIMDVec<T, SIMD_WIDTH> match(
    const SIMDVec<T, SIMD_WIDTH> &a, const SIMDVec<T, SIMD_WIDTH> &b)
  {
    // max(min(a,b), -max(a,b))
    return max(min(a, b), neg(max(a, b)));
  }
};

// match operator for ASC without inversion
// (slow, for tests only, effect: addition of 1 only in denominator)
template <typename T, int SIMD_WIDTH>
struct MatchAscDirect
{
  static SIMD_INLINE SIMDVec<T, SIMD_WIDTH> match(
    const SIMDVec<T, SIMD_WIDTH> &a, const SIMDVec<T, SIMD_WIDTH> &b)
  {
    // |a| + |b| - (|a + b| - |a - b|)
    return sub(add(abs(a), abs(b)), sub(abs(add(a, b)), abs(sub(a, b))));
  }
};

// match operator for NSADA
template <typename T, int SIMD_WIDTH>
struct MatchAbsDiffAbs
{
  static SIMD_INLINE SIMDVec<T, SIMD_WIDTH> match(
    const SIMDVec<T, SIMD_WIDTH> &a, const SIMDVec<T, SIMD_WIDTH> &b)
  {
    return abs(sub(abs(a), abs(b)));
  }
};

// match operator for NSAD-NSADA mix
template <typename T, int SIMD_WIDTH>
struct MatchAvgAbsDiffAbsDiffAbs
{
  static SIMD_INLINE SIMDVec<T, SIMD_WIDTH> match(
    const SIMDVec<T, SIMD_WIDTH> &a, const SIMDVec<T, SIMD_WIDTH> &b)
  {
    return avg(abs(sub(a, b)), abs(sub(abs(a), abs(b))));
  }
};

// ---------------------------------------------------------------------------
// pairing of image columns 1 (without re-arrange)
// ---------------------------------------------------------------------------

// 1 is faster!
#define CPM1_ACCUMULATE_IN_MEMORY 1
// (if CPM1_ACCUMULATE_IN_MEMORY == 1) 1 is slightly faster
#define CPM1_SPLIT_CLEAR 1

// loop structurefrom old scalarProductSIMDWord2: accumulation in memory
// (I tried a version of this where I replaced pointer by index access, but
// it was a bit slower)

// 2x Panorama -> DefaultAddOn
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class MatchOperator>
void columnPairMatching1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  // same horizon required, but vertical resolution can differ
  if (l1.addOn.horizon != l2.addOn.horizon)
    throw SIMDException("columnPairMatching1", "panorama horizons differ");
  if (overlap(l1, dir) || overlap(l2, dir))
    throw SIMDException("columnPairMatching1",
                        "input and output images overlap");
  if (!l1.sameSize(l2))
    throw SIMDException("columnPairMatching1",
                        "input images have different size");
  SIMDVec<T, SIMD_WIDTH> l1_Nx, l2_N;
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  // image size
  const int sz = l1.size, w = l1.w;
  dir.resize(w, w);
  const T *const l1data = l1.data, *const l2data = l2.data;
#if CPM1_ACCUMULATE_IN_MEMORY
  // accumulate in memory: faster!
  // set entire result to zero (we accumulate here)
#if !CPM1_SPLIT_CLEAR
  dir.clear();
#endif
  // start at first row in dir
  T *dirRow = dir.data;
  // slowest loop: single pixel access in l1
  // l1Col points to the top element of each column in l1)
  for (const T *l1Col = l1data; l1Col < l1data + w; l1Col++) {
    // start at first row in l2
    const T *l2El = l2data;
#if CPM1_SPLIT_CLEAR
    memset(dirRow, 0x00, w * sizeof(T));
#endif
    // middle loop: direction term accumulation
    // l1El goes through all elements in the current column of l1)
    for (const T *l1El = l1Col; l1El < l1Col + sz; l1El += w) {
      // this is the common term for all row elements
      // SIMDVec<T,SIMD_WIDTH>
      l1_Nx = set1<T, SIMD_WIDTH>(*l1El);
      // pointer dirEl runs through all row elements
      // fastest loop: step-wise computation of direction term
      for (T *dirEl = dirRow; dirEl < dirRow + w; dirEl += simd_elems) {
        // load
        // SIMDVec<T,SIMD_WIDTH>
        l2_N = load<SIMD_WIDTH>(l2El);
        // sum up match
        store(dirEl,
              adds(load<SIMD_WIDTH>(dirEl), MatchOperator::match(l1_Nx, l2_N)));
        // next element in l2
        l2El += simd_elems;
      }
    }
    // next row in result
    dirRow += w;
  }
#else  // CPM1_ACCUMULATE_IN_MEMORY
  // accumulate in register s_N
  SIMDVec<T, SIMD_WIDTH> s_N;
  // loop structure from scalarProductSIMDWord1 in ImageSSE.C, modified
  T *dirEl = dir.data;
  for (const T *l1Col = l1data; l1Col < l1data + w; l1Col++) {
    for (const T *l2Col = l2data; l2Col < l2data + w;
         l2Col += simd_elems, dirEl += simd_elems) {
      // reset summation variable
      s_N = setzero<T, SIMD_WIDTH>();
      for (const T *l1El = l1Col, *l2El = l2Col; l1El < l1Col + sz;
           l1El += w, l2El += w) {
        // get factor from l1 (this is probably what makes it slow)
        l1_Nx = set1<T, SIMD_WIDTH>(*l1El);
        // get elements from l2
        l2_N = load<SIMD_WIDTH>(l2El);
        s_N  = adds(s_N, MatchOperator::match(l1_Nx, l2_N));
      }
      // store sum
      store(dirEl, s_N);
    }
  }
#endif // CPM1_ACCUMULATE_IN_MEMORY
}

// NSAD term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAbsDiffTerm1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching1<T, SIMD_WIDTH, SIMD_ALIGN, MatchAbsDiff<T, SIMD_WIDTH>>(
    l1, l2, dir);
}

// ASC+ term (in second version using min/max)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void asc2Term1(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
               const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
               SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  SIMDVec<T, SIMD_WIDTH> dirV;
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  columnPairMatching1<T, SIMD_WIDTH, SIMD_ALIGN, MatchAsc2<T, SIMD_WIDTH>>(
    l1, l2, dir);
  // multiply result by 2 since the max-min expression is only half as
  // large as the one in ascDirTerm1
  T *const dirdata = dir.data;
  const int sz     = dir.size;
  for (T *dirp = dirdata; dirp < dirdata + sz; dirp += simd_elems) {
    // SIMDVec<T,SIMD_WIDTH>
    dirV = load<SIMD_WIDTH>(dirp);
    store(dirp, adds(dirV, dirV));
  }
}

// NSADA term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAbsDiffAbsTerm1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching1<T, SIMD_WIDTH, SIMD_ALIGN,
                      MatchAbsDiffAbs<T, SIMD_WIDTH>>(l1, l2, dir);
}

// NSAD-NSADA mix term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAvgAbsDiffAbsDiffAbsTerm1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching1<T, SIMD_WIDTH, SIMD_ALIGN,
                      MatchAvgAbsDiffAbsDiffAbs<T, SIMD_WIDTH>>(l1, l2, dir);
}

// ---------------------------------------------------------------------------
// columnPairMatchingTilt1 (ignoring invalid values)
// contributed by Annika Hoffmann
// ---------------------------------------------------------------------------

// loop structure from old scalarProductSIMDWord2: accumulation in memory
// consider invalid pixels
// additional parameters for precomputed sum of absolute values
// (colMeasRaw1, colMeasRaw2) - here the corresponding value has to be picked
// in colMeasure1 and colMeasure2 the used colMeasRaw-values are saved
// (for normalization step)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class MatchOperator>
void columnPairMatchingTilt1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeasure1,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeasure2,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeasure12)
{
  // same horizon required, but vertical resolution can differ
  // (due to magnification of images)
  if (l1.addOn.horizon != l2.addOn.horizon)
    throw SIMDException("columnPairMatchingTilt1", "panorama horizons differ");
  if (overlap(l1, dir) || overlap(l2, dir))
    throw SIMDException("columnPairMatchingTilt1",
                        "input and output images overlap");
  if (!l1.sameSize(l2))
    throw SIMDException("columnPairMatchingTilt1",
                        "input images have different size");

  // invalidVal: we always use maximal value of type T to indicate
  // invalid pixels
  T invalidVal                      = SIMDTypeInfo<T>::max();
  SIMDVec<T, SIMD_WIDTH> invalidVec = set1<T, SIMD_WIDTH>(invalidVal);
  // image size
  const int sz = l1.size, w = l1.w /*, h=l1.h */;
  dir.resize(w, w);
  // colMeasure1 and colMeasure2: may be required for special
  // distance measures where the individual sums are needed, for
  // NSAD etc. where we divde by the sum, colMeasure12 (their sum)
  // would be sufficient; colMeasure1 and colMeasure2 need to be of
  // size w x w, since individual combinations of columns are
  // required
  colMeasure1.resize(w, w);
  colMeasure2.resize(w, w);
  colMeasure12.resize(w, w);
  SIMDVec<T, SIMD_WIDTH> l1_Nx, l2_N;
  const int simd_elems  = numSIMDVecElements<T, SIMD_WIDTH>();
  const T *const l1data = l1.data, *const l2data = l2.data;
  SIMDVec<T, SIMD_WIDTH> zero, one, summand, s_N, colMeasA, colMeasB, colMeasAB;
  // zero-vector for zero-addition, if element is invalid
  zero = setzero<T, SIMD_WIDTH>();
  // vector with 1 elements for incrementing counter
  one = setunity<T, SIMD_WIDTH>();
  // accumulated absolute values
  const T *pAccuA = colMeasRaw1.data;
  const T *pAccuB = colMeasRaw2.data;

  // individual
  T *pcolMeasA  = colMeasure1.data;
  T *pcolMeasB  = colMeasure2.data;
  T *pcolMeasAB = colMeasure12.data;

  // counter for last valid index in image1/2 columns
  // we use the same type as the image such that we have the same number
  // of elements (should be ok even for small types, since row number is
  // typically small)
  SIMDVec<T, SIMD_WIDTH> firstInvalidRow;
  // boolean vector: invalid elements in a given row
  SIMDVec<T, SIMD_WIDTH> invalidElem;
  // buffer for firstInvalidRow
  T firstInvalidRowBuf[simd_elems] SIMD_ATTR_ALIGNED(SIMD_ALIGN);
  // loop structure from scalarProductSIMDWord1 in ImageSSE.C, modified
  T *dirEl = dir.data;
  // l1_offset holds the index of the currently processed column in l1data
  // l2_offset hold the index of the currently processed column in l2data
  int l1_offset, l2_offset;
  // l1Col runs through the top element of each column in l1
  l1_offset = 0;
  for (const T *l1Col = l1data; l1Col < l1data + w; l1Col++, l1_offset++) {
    // reset l2_offset; we start in column 0 in l2data
    l2_offset = 0;
    // l2Col runs through top element of columns in l2,
    // in steps of simd_elems, which are used in a vector
    for (const T *l2Col = l2data; l2Col < l2data + w; l2Col += simd_elems,
                 dirEl += simd_elems, pcolMeasA += simd_elems,
                 pcolMeasB += simd_elems, pcolMeasAB += simd_elems,
                 l2_offset += simd_elems) {
      // reset summation variable
      s_N = setzero<T, SIMD_WIDTH>();
      // initial assumption: all rows are invalid
      firstInvalidRow = setzero<T, SIMD_WIDTH>();
      // l1El runs through each element in the current column of l1
      // l2El runs through each element in the current column of l2
      for (const T *l1El = l1Col, *l2El = l2Col; l1El < l1Col + sz;
           l1El += w, l2El += w) {
        // if pixel in image1 is invalid, we leave the loop
        if (*l1El == invalidVal) break;

        // get factor from l1
        l1_Nx = set1<T, SIMD_WIDTH>(*l1El);
        // get elements from l2
        l2_N = load<SIMD_WIDTH>(l2El);

        // check whether all pixels in l2_N are invalid, leave the loop
        invalidElem = cmpeq(l2_N, invalidVec);

        if (test_all_ones(invalidElem)) break;

        // add 1 to all columns with valid elements
        firstInvalidRow = add(firstInvalidRow, ifelse(invalidElem, zero, one));

        // here we know that l1El is valid, and some elements of l2N are also
        // valid
        // compute sum, but consider invalid values
        summand = ifelse(invalidElem, zero, MatchOperator::match(l1_Nx, l2_N));
        s_N     = adds(s_N, summand);
      }
      // store sum
      store(dirEl, s_N);
      // store firstInvalidRow
      store(firstInvalidRowBuf, firstInvalidRow);
      // firstInvalidRowBuf contains the index of the first invalid row for
      // each of simd_elems columns
      for (int i = 0; i < simd_elems; i++) {
        int lastValidRow = firstInvalidRowBuf[i] - 1;
        // we don't allow column pairs without any valid rows
        // if colMeasA and colMeasB are not required, we could add
        // here instead of afterwards
        // if there are no valid pixels (in the column pair), we just take
        // row 0 from the accumalated row data, since s_N element is 0,
        // normalized measure will be zero
        if (lastValidRow < 0)
          // should not happen, but just to be on the safe side
          fprintf(stderr,
                  "columnPairMatchingTilt1: lastValidRow (%d) < 0 "
                  "at l1Col %d l2col %d i %d\n",
                  lastValidRow, int(l1Col - l1data), int(l2Col - l2data), i);
        int validRowOff  = (lastValidRow >= 0) ? lastValidRow * w : 0;
        *(pcolMeasA + i) = *(pAccuA + validRowOff + l1_offset);
        *(pcolMeasB + i) = *(pAccuB + validRowOff + l2_offset + i);
      }
      // compute sum of colMeasA and colMeasB -> later used for normalization
      colMeasA  = load<SIMD_WIDTH>(pcolMeasA);
      colMeasB  = load<SIMD_WIDTH>(pcolMeasB);
      colMeasAB = adds(colMeasA, colMeasB);

      store(pcolMeasAB, colMeasAB);
    }
  }
}

// ---------------------------------------------------------------------------
// pairing of image columns 2 (with swizzle)
// ---------------------------------------------------------------------------

// 0 is faster
#define CPM2_ACCUMULATE_IN_MEMORY 0
// (if CPM2_ACCUMULATE_IN_MEMORY == 1) is slightly faster
#define CPM2_SPLIT_CLEAR 1
// not much difference
#define CPM2_ORDER_VH 0

// 2xPanorama -> DefaultAddOn
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class MatchOperator>
void columnPairMatching2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgV,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgH, int w, int wSPS,
  int step, int n, int blkSize, const int *const modn,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  /*
  printf("columnPairMatching2: SIMD_WIDTH=%d, w=%d, wSPS=%d, n=%d, step=%d,
  imgV.w=%d, imgV.h=%d\n", SIMD_WIDTH, w, wSPS, n, step, imgV.w, imgV.h);
  */
  // same horizon required, but vertical resolution can differ
  if (imgV.addOn.horizon != imgH.addOn.horizon)
    throw SIMDException("columnPairMatching1", "panorama horizons differ");
  if (overlap(imgV, dir) || overlap(imgH, dir))
    throw SIMDException("columnPairMatching2",
                        "input and output images overlap");
  if (!imgV.sameSize(imgH))
    throw SIMDException("columnPairMatching2",
                        "input images have different size");
  if (imgV.w != wSPS)
    throw SIMDException("columnPairMatching2", "invalid size of input images");
  // image size (includes complementary columns)
  const int szImg = imgV.size;
  // complementary columns: number of elements
  const int nElems = numSIMDVecElements<T, SIMD_WIDTH>();
  // short-hands to image data
  const T *const dataV = imgV.data, *const dataH = imgH.data;
  // resize result array
  dir.resize(w, w);
  // shorthand to result data
  T *dirData = dir.data;
  // step x w
  const int step_x_w = step * w;
  // vectors
  SIMDVec<T, SIMD_WIDTH> vecH, vecV;
  // 2 accumulation modes
#if CPM2_ACCUMULATE_IN_MEMORY
  // set entire result to zero (we accumulate here)
  // we also tested a version where the clear is split (slightly faster)
#if !CPM2_SPLIT_CLEAR
  // clear in one go
  dir.clear();
#endif // CPM2_SPLIT_CLEAR
#else  // CPM2_ACCUMULATE_IN_MEMORY
  // we accumulate in a register (vecS)
  SIMDVec<T, SIMD_WIDTH> vecS;
#endif // CPM2_ACCUMULATE_IN_MEMORY
       // two different loop orders
#if CPM2_ORDER_VH
  // order: vertical->horizontal
  // vertical block (i.e. block in vertical image imgV)
  const T *imgBlkV = dataV;
  for (int iBlkV = 0; iBlkV < step; iBlkV++, imgBlkV += blkSize) {
    // horizontal block (i.e. block in horizontal image imgH)
    const T *imgBlkH = dataH;
    T *dirBlkH       = dirData;
    for (int iBlkH = 0; iBlkH < step;
         iBlkH++, imgBlkH += blkSize, dirBlkH += n) {
#else
  // order: horizontal->vertical
  // horizontal block (i.e. block in horizontal image imgH)
  const T *imgBlkH = dataH;
  T *dirBlkH       = dirData;
  for (int iBlkH = 0; iBlkH < step; iBlkH++, imgBlkH += blkSize, dirBlkH += n) {
    // vertical block (i.e. block in vertical image imgV)
    const T *imgBlkV = dataV;
    for (int iBlkV = 0; iBlkV < step; iBlkV++, imgBlkV += blkSize) {
#endif
      // for computation of start value for iDelta
      // TODO: LUT for iDelta0?
      int iDelta0 = (iBlkV - iBlkH + step) % step;
      T *dirRow   = dirBlkH + iDelta0 * w;
      // go through all iDelta belonging to this H/V block
      for (int iDelta = iDelta0; iDelta < w;
           iDelta += step, dirRow += step_x_w) {
        // compute horizontal offset
        // TODO: LUT for off?
        int off = (iBlkH + iDelta) / step;
#if CPM2_ACCUMULATE_IN_MEMORY
#if CPM2_SPLIT_CLEAR
        // clear is split
        memset(dirRow, 0x00, n * sizeof(T));
#endif // CPM2_SPLIT_CLEAR
       // go through all rows and accumulate in memory
        for (const T *imgRowH = imgBlkH, *imgRowV = imgBlkV;
             imgRowH < dataH + szImg; imgRowH += wSPS, imgRowV += wSPS) {
          // go through all horizontal indices in imgH
          for (int iOffH = 0; iOffH < n; iOffH += nElems) {
            // compute index in l2
            int iOffV = modn[iOffH + off];
            // load inputs
            vecH     = load<SIMD_WIDTH>(imgRowH + iOffH);
            vecV     = loadu<SIMD_WIDTH>(imgRowV + iOffV);
            T *dirEl = dirRow + iOffH;
            // compute and accumulate
            store(dirEl, adds(load<SIMD_WIDTH>(dirEl),
                              MatchOperator::match(vecV, vecH)));
          }
        }
#else  // CPM2_ACCUMULATE_IN_MEMORY
       // accumulate in register: faster!
       // (especially since dir.clear() is not executed)
       // go through all horizontal indices in imgH
        T *dirEl = dirRow;
        for (int iOffH = 0; iOffH < n; iOffH += nElems, dirEl += nElems) {
          // compute index in l2
          int iOffV = modn[iOffH + off];
          // go through all rows and accumulate in vecS
          vecS = setzero<T, SIMD_WIDTH>();
          for (const T *imgRowH = imgBlkH + iOffH, *imgRowV = imgBlkV + iOffV;
               imgRowH < dataH + szImg; imgRowH += wSPS, imgRowV += wSPS) {
            // load inputs
            vecH = load<SIMD_WIDTH>(imgRowH);
            vecV = loadu<SIMD_WIDTH>(imgRowV);
            // compute and accumulate
            vecS = adds(vecS, MatchOperator::match(vecV, vecH));
          }
          // stream_store and sfence (outside inner loop) was slower
          store(dirEl, vecS);
        }
#endif // CPM2_ACCUMULATE_IN_MEMORY
      }
    }
  }
}

// NSAD term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAbsDiffTerm2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2, int w, int wSPS,
  int step, int n, int blkSize, const int *const modn,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching2<T, SIMD_WIDTH, SIMD_ALIGN, MatchAbsDiff<T, SIMD_WIDTH>>(
    l1, l2, w, wSPS, step, n, blkSize, modn, dir);
}

// ASC+ term (in second version using min/max)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void asc2Term2(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
               const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2, int w,
               int wSPS, int step, int n, int blkSize, const int *const modn,
               SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  SIMDVec<T, SIMD_WIDTH> dirV;
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  columnPairMatching2<T, SIMD_WIDTH, SIMD_ALIGN, MatchAsc2<T, SIMD_WIDTH>>(
    l1, l2, w, wSPS, step, n, blkSize, modn, dir);
  // multiply result by 2 since the max-min expression is only half as
  // large as the one in ascDirTerm1
  T *const dirdata = dir.data;
  const int sz     = dir.size;
  for (T *dirp = dirdata; dirp < dirdata + sz; dirp += simd_elems) {
    // SIMDVec<T,SIMD_WIDTH>
    dirV = load<SIMD_WIDTH>(dirp);
    store(dirp, adds(dirV, dirV));
  }
}

// ASC direct term (slow, for tests only)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void ascDirectTerm(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
                   const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2,
                   int w, int wSPS, int step, int n, int blkSize,
                   const int *const modn,
                   SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching2<T, SIMD_WIDTH, SIMD_ALIGN, MatchAscDirect<T, SIMD_WIDTH>>(
    l1, l2, w, wSPS, step, n, blkSize, modn, dir);
}

// NSADA term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAbsDiffAbsTerm2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2, int w, int wSPS,
  int step, int n, int blkSize, const int *const modn,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching2<T, SIMD_WIDTH, SIMD_ALIGN,
                      MatchAbsDiffAbs<T, SIMD_WIDTH>>(l1, l2, w, wSPS, step, n,
                                                      blkSize, modn, dir);
}

// NSAD-NSADA mix term
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
void normSumAvgAbsDiffAbsDiffAbsTerm2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l1,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &l2, int w, int wSPS,
  int step, int n, int blkSize, const int *const modn,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &dir)
{
  columnPairMatching2<T, SIMD_WIDTH, SIMD_ALIGN,
                      MatchAvgAbsDiffAbsDiffAbs<T, SIMD_WIDTH>>(
    l1, l2, w, wSPS, step, n, blkSize, modn, dir);
}

// ---------------------------------------------------------------------------
// normalized distance measure (denominator terms are added)
// ---------------------------------------------------------------------------

// normalized distance measure where denominator terms are added
template <int invertJointMeasure, typename SPSType, typename T, int SIMD_WIDTH,
          int SIMD_ALIGN>
void computeScalePlaneNormalizedAdd1(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &ssMeasure,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &cvMeasure, int wSPS,
  const int *const modulo2wpw, const int *const shuffle, double postScale,
  SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane)
{
  if (overlap(jointMeasure, scalePlane) || overlap(ssMeasure, scalePlane) ||
      overlap(cvMeasure, scalePlane))
    throw SIMDException("computeScalePlaneNormalizedAdd1",
                        "input and output images overlap");
  if (!ssMeasure.sameSize(cvMeasure))
    throw SIMDException("computeScalePlaneNormalizedAdd1",
                        "ss/cv measure images have different size");
  if (jointMeasure.w != jointMeasure.h)
    throw SIMDException("computeScalePlaneNormalizedAdd1",
                        "joint measure image is not square-sized");
  if (jointMeasure.w != ssMeasure.w)
    throw SIMDException("computeScalePlaneNormalizedAdd1",
                        "joint/ss measure differ in width");
  const int w = jointMeasure.w;
  // number of elements in a vector relating to the input SIMDImages
  const int vecElemsIn = numSIMDVecElements<T, SIMD_WIDTH>();
  // number of input vectors required to fill single output vector
  const int numInVecs = numInputSIMDVecs<SPSType, T>();
  // number of output vectors required to take up single input vector
  const int numOutVecs = numOutputSIMDVecs<SPSType, T>();
  // number of input and output elements
  const int numElems = numSIMDVecsElements<SPSType, T, SIMD_WIDTH>();
  // width must be multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("computeScalePlaneNormalizedAdd1",
                        "width must be multiple of numElems");
  SIMDVec<T, SIMD_WIDTH> cvMeasxN, measSum[numInVecs], joint[numInVecs];
  SIMDVec<SPSType, SIMD_WIDTH> spsVal[numOutVecs];

  const int *modulo2woff = modulo2wpw;
  // destination scale plane
  scalePlane.resize(wSPS, w);
  SPSType *destP = scalePlane.data;
  // pointer to dirTerm
  const T *jointPtr = jointMeasure.data;
  // i corresponds to theta'
  for (int i = 0; i < w; i++) {
    // we take a single value from cvMeasure and expand it to
    // vecElemsIn identical values (numInVecs vectors);
    // here we add 1 to avoid division by zero!
    // TODO: should we make this epsilon offset variable?
    // SIMDVec<T,SIMD_WIDTH>
    cvMeasxN             = set1<T, SIMD_WIDTH>(cvMeasure.data[i] + 1);
    const T *ssMeasPtr   = ssMeasure.data;
    const int *shufflepj = shuffle;
    // j corresponds to theta
    for (int j = 0; j < w; j += numElems) {
      // SIMDVec<T,SIMD_WIDTH> measSum[numInVecs], joint[numInVecs];
      for (int vi = 0, off = 0; vi < numInVecs; vi++, off += vecElemsIn) {
        // sum of two values in denominator
        measSum[vi] = adds(cvMeasxN, load<SIMD_WIDTH>(ssMeasPtr + off));
        // denominator
        joint[vi] = load<SIMD_WIDTH>(jointPtr + off);
        // inversion for correlation-based joint measures from [-1,1]
        if (invertJointMeasure) joint[vi] = subs(measSum[vi], joint[vi]);
      }
      // SIMDVec<SPSType,SIMD_WIDTH> spsVal[numOutVecs];
      // division and multiplication by postScales
      fdivmul(joint, measSum, postScale, spsVal);
      // store to memory
      // TODO: would store of SPS be faster with extract instead of writing
      // TODO: to memory first?)
      // TODO: would stream operation be faster?
#if 1
      SPSType Cstore[numElems] SIMD_ATTR_ALIGNED(SIMD_ALIGN);
      store(Cstore, spsVal, numOutVecs);
      // now we do the "slanted" insertion into dest
      // TODO: store of SPS is quite time-consuming due to inefficient
      // TODO: cache use, especially if single bytes are written
      const int *modulo2woffmj = modulo2woff - j;
      // we have to go through all elements separately
      for (int k = 0; k < numElems; k++, modulo2woffmj--)
        *(destP + shufflepj[k] + *modulo2woffmj) = Cstore[k];
#else
      // test how vector store would accelerate this code:
      // 3400 us -> 2300 us!
      // code doesn't produce correct homevectors, of course
      store(destP, spsVal, numOutVecs);
      destP += numElems;
#endif
      jointPtr += numElems;
      ssMeasPtr += numElems;
      shufflepj += numElems;
    }
    modulo2woff++;
  }
}

// ---------------------------------------------------------------------------
// multi dimensional sigmoid on normalized distance measures
// (denominator terms are added)
// ---------------------------------------------------------------------------

// normalized distance measure where denominator terms are added
template <int invertJointMeasure, int NPREPROCS, typename SPSType, typename T,
          int SIMD_WIDTH, int SIMD_ALIGN>
void computeScalePlaneSigmoidNormalizedAdd1(
  const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN> &jointMeasureStk,
  const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN> &ssMeasureStk,
  const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN> &cvMeasureStk, int wSPS,
  const int *const modulo2wpw, const int *const shuffle,
  const std::vector<double> &sigmoidW, const std::vector<double> &sigmoidW0,
  double postScale, SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane)
{
  if (overlap(jointMeasureStk, scalePlane) ||
      overlap(ssMeasureStk, scalePlane) || overlap(cvMeasureStk, scalePlane))
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "input and output images overlap");
  if (!ssMeasureStk.sameSize(cvMeasureStk))
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "ss/cv measure images have different size");
  if (jointMeasureStk.w != jointMeasureStk.h)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "joint measure image is not square-sized");
  if (jointMeasureStk.w != ssMeasureStk.w)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "joint/ss measure differ in width");
  if (jointMeasureStk.numPlanes != NPREPROCS ||
      ssMeasureStk.numPlanes != NPREPROCS ||
      cvMeasureStk.numPlanes != NPREPROCS)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "input number of planes differ the defined number");
  if (sigmoidW.size() != (size_t) NPREPROCS)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "sigmoidW has wrong size");
  if (sigmoidW0.size() != (size_t) NPREPROCS)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "sigmoidW0 has wrong size");
  const int w = jointMeasureStk.w;
  // number of elements in a vector relating to the input SIMDImages
  const int vecElemsIn = numSIMDVecElements<T, SIMD_WIDTH>();
  // number of input vectors required to fill single output vector
  const int numInVecs = numInputSIMDVecs<SPSType, T>();
  // number of output vectors required to take up single input vector
  const int numOutVecs = numOutputSIMDVecs<SPSType, T>();
  // number of input and output elements
  const int numElems = numSIMDVecsElements<SPSType, T, SIMD_WIDTH>();
  // width must be multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("computeScalePlaneSigmoidNormalizedAdd1",
                        "width must be multiple of numElems");
  SIMDVec<T, SIMD_WIDTH> cvMeasxN[NPREPROCS], measSum[NPREPROCS][numInVecs],
    joint[NPREPROCS][numInVecs];
  SIMDVec<SPSType, SIMD_WIDTH> spsVal[numOutVecs];

  const int *modulo2woff = modulo2wpw;
  // destination scale plane
  scalePlane.resize(wSPS, w);
  SPSType *destP = scalePlane.data;
  // pointer to dirTerm
  const T *jointPtr[NPREPROCS]; // = jointMeasure.data;
  for (int pp = 0; pp < NPREPROCS; pp++) {
    jointPtr[pp] = jointMeasureStk[pp].data;
  }
  // i corresponds to theta'
  for (int i = 0; i < w; i++) {
    // we take a single value from cvMeasure and expand it to
    // vecElemsIn identical values (numInVecs vectors);
    // here we add 1 to avoid division by zero!
    // TODO: should we make this epsilon offset variable?
    // SIMDVec<T,SIMD_WIDTH>
    const T *ssMeasPtr[NPREPROCS];
    for (int pp = 0; pp < NPREPROCS; pp++) {
      cvMeasxN[pp]  = set1<T, SIMD_WIDTH>(cvMeasureStk[pp].data[i] + 1);
      ssMeasPtr[pp] = ssMeasureStk[pp].data;
    }
    const int *shufflepj = shuffle;
    // j corresponds to theta
    for (int j = 0; j < w; j += numElems) {
      for (int pp = 0; pp < NPREPROCS; pp++) {
        for (int vi = 0, off = 0; vi < numInVecs; vi++, off += vecElemsIn) {
          // sum of two values in denominator
          measSum[pp][vi] =
            adds(cvMeasxN[pp], load<SIMD_WIDTH>(ssMeasPtr[pp] + off));
          // denominator
          joint[pp][vi] = load<SIMD_WIDTH>(jointPtr[pp] + off);
          // inversion for correlation-based joint measures from [-1,1]
          if (invertJointMeasure)
            joint[pp][vi] = subs(measSum[pp][vi], joint[pp][vi]);
        }
      }
      // division, sigmoid and scaling as float
      fdivMsigmoidmul<NPREPROCS, numInVecs, SPSType, T, SIMD_WIDTH>(
        joint, measSum,
        // conversion to C array
        &sigmoidW[0], &sigmoidW0[0], postScale, spsVal);

      // store to memory
      // TODO: would store of SPS be faster with extract instead of writing
      // TODO: to memory first?)
      // TODO: would stream operation be faster?
#if 1
      SPSType Cstore[numElems] SIMD_ATTR_ALIGNED(SIMD_ALIGN);
      store(Cstore, spsVal, numOutVecs);
      // now we do the "slanted" insertion into dest
      // TODO: store of SPS is quite time-consuming due to inefficient
      // TODO: cache use, especially if single bytes are written
      const int *modulo2woffmj = modulo2woff - j;
      // we have to go through all elements separately
      for (int k = 0; k < numElems; k++, modulo2woffmj--)
        *(destP + shufflepj[k] + *modulo2woffmj) = Cstore[k];
#else
      // test how vector store would accelerate this code:
      // 3400 us -> 2300 us!
      // code doesn't produce correct homevectors, of course
      store(destP, spsVal, numOutVecs);
      destP += numElems;
#endif
      for (int pp = 0; pp < NPREPROCS; pp++) {
        jointPtr[pp] += numElems;
        ssMeasPtr[pp] += numElems;
      }
      shufflepj += numElems;
    }
    modulo2woff++;
  }
}

// ---------------------------------------------------------------------------
// normalized distance measure (denominator terms are added), swizzled version
// ---------------------------------------------------------------------------

// 0 is faster
#define CSPNA2_LOOP_AS_IN_CPM2 0

// normalized distance measure where denominator terms are added
template <int invertJointMeasure, typename SPSType, typename T, int SIMD_WIDTH,
          int SIMD_ALIGN>
void computeScalePlaneNormalizedAdd2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &ssMeasure,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &cvMeasure, int w, int wSPS,
  int stepAlpha, int nAlpha, int blkSize, const int *const modNAlpha,
  const int *const shuffleModulo, const int *const unshuffle, double postScale,
  SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane)
{
  if (overlap(jointMeasure, scalePlane) || overlap(ssMeasure, scalePlane) ||
      overlap(cvMeasure, scalePlane))
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "input and output images overlap");
  if (!ssMeasure.sameSize(cvMeasure))
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "ss/cv measure images have different size");
  if (ssMeasure.h != 1)
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "ss/cv measure images should have height 1");
  if (ssMeasure.w != wSPS)
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "ss/cv measure images have invalid width");
  if (jointMeasure.w != jointMeasure.h)
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "joint measure image is not square-sized");
  if (jointMeasure.w != w)
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "joint measure has invalid dimensions");
  // number of elements in a vector relating to the input SIMDImages
  const int vecElemsIn = numSIMDVecElements<T, SIMD_WIDTH>();
  // number of input vectors required to fill single output vector
  const int numInVecs = numInputSIMDVecs<SPSType, T>();
  // number of output vectors required to take up single input vector
  const int numOutVecs = numOutputSIMDVecs<SPSType, T>();
  // number of input and output elements
  const int numElems = numSIMDVecsElements<SPSType, T, SIMD_WIDTH>();
  // width must be multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("computeScalePlaneNormalizedAdd2",
                        "width must be multiple of numElems");
  // SIMD vectors
  SIMDVec<T, SIMD_WIDTH> cvMeas[numInVecs], ssMeas[numInVecs],
    measSum[numInVecs], jointMeas[numInVecs];
  SIMDVec<SPSType, SIMD_WIDTH> spsVal[numOutVecs];
  // add 1 to cvMeasure to avoid div by zero
  // TODO: should we make this epsilon offset variable?
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> cvMeasureP1;
  adds(cvMeasure, T(1), cvMeasureP1);
  // shorthands
  const T *const cvData = cvMeasureP1.data, *const ssData = ssMeasure.data;
  // destination scale plane
  scalePlane.resize(wSPS, w);

  // 10. Feb 16 (rm): made this a C-if instead of a cpp-if to avoid
  // unused-pameter warnings
  if (CSPNA2_LOOP_AS_IN_CPM2) {
    // same loop structure as in columnPairMatching2 (version: vertical->horiz.)
    // step x w
    const int step_x_w       = stepAlpha * w;
    const int step_x_wSPS    = stepAlpha * wSPS;
    SPSType *const destData  = scalePlane.data;
    const T *const jointData = jointMeasure.data;
    // vertical block (i.e. block in vertical image imgV)
    const T *imgBlkCV = cvData;
    for (int iBlkCV = 0; iBlkCV < stepAlpha; iBlkCV++, imgBlkCV += blkSize) {
      // horizontal block (i.e. block in horizontal image imgH)
      const T *imgBlkSS = ssData;
      SPSType *destBlk  = destData;
      const T *jointBlk = jointData;
      for (int iBlkSS = 0; iBlkSS < stepAlpha; iBlkSS++, imgBlkSS += blkSize,
               jointBlk += nAlpha, destBlk += blkSize) {
        // for computation of start value for iDelta
        // TODO: LUT for iDelta0?
        int iDelta0       = (iBlkCV - iBlkSS + stepAlpha) % stepAlpha;
        SPSType *destRow  = destBlk + iDelta0 * wSPS;
        const T *jointRow = jointBlk + iDelta0 * w;
        // go through all iDelta belonging to this H/V block
        for (int iDelta = iDelta0; iDelta < w; iDelta += stepAlpha,
                 jointRow += step_x_w, destRow += step_x_wSPS) {
          // compute horizontal offset
          // TODO: LUT for off?
          int off = (iBlkSS + iDelta) / stepAlpha;
          // go through all horizontal indices in imgH
          SPSType *destEl  = destRow;
          const T *jointEl = jointRow;
          for (int iOffSS = 0; iOffSS < nAlpha;
               iOffSS += numElems, jointEl += numElems, destEl += numElems) {
            // compute index in l2
            int iOffCV       = modNAlpha[iOffSS + off];
            const T *imgSSEl = imgBlkSS + iOffSS, *imgCVEl = imgBlkCV + iOffCV;
            for (int vi = 0, vecOff = 0; vi < numInVecs;
                 vi++, vecOff += vecElemsIn) {
              // load measures
              ssMeas[vi] = load<SIMD_WIDTH>(imgSSEl + vecOff);
              cvMeas[vi] = loadu<SIMD_WIDTH>(imgCVEl + vecOff);
              // puts("cv loaded");
              jointMeas[vi] = load<SIMD_WIDTH>(jointEl + vecOff);
              // sum of two values in denominator
              measSum[vi] = adds(cvMeas[vi], ssMeas[vi]);
              // inversion for correlation-based joint measures from [-1,1]
              if (invertJointMeasure)
                jointMeas[vi] = subs(measSum[vi], jointMeas[vi]);
            }
            // division and multiplication by postScales
            fdivmul(jointMeas, measSum, postScale, spsVal);
            store(destEl, spsVal, numOutVecs);
          }
        }
      }
    }
  }

  else { // CSPNA2_LOOP_AS_IN_CPM2

    // simplified loop structure (directly using shuffle/unshuffle): faster!
    SPSType *destRow = scalePlane.data;
    // elements in joint.data are contiguous, simpler loop counting possible
    const T *jointEl = jointMeasure.data;
    // go through all row indices
    for (int iDelta = 0; iDelta < w; iDelta++, destRow += wSPS) {
      const int *const shuffleModulo_p_iDelta = shuffleModulo + iDelta;
      // go through all column indices (block-wise)
      SPSType *destBlk = destRow;
      for (const T *ssBlk = ssData; ssBlk < ssData + wSPS;
           ssBlk += blkSize, destBlk += blkSize) {
        SPSType *destEl = destBlk;
        for (const T *ssEl = ssBlk; ssEl < ssBlk + nAlpha;
             ssEl += numElems, destEl += numElems, jointEl += numElems) {
          // read multiple input vectors
          // we can compute iCV here since we have enough wrap-around columns
          // such that we can add off below
          int iSS       = ssEl - ssData;
          int iCV       = shuffleModulo_p_iDelta[unshuffle[iSS]];
          const T *cvEl = cvData + iCV;
          for (int vi = 0, off = 0; vi < numInVecs; vi++, off += vecElemsIn) {
            // load measures
            ssMeas[vi] = load<SIMD_WIDTH>(ssEl + off);
            cvMeas[vi] = loadu<SIMD_WIDTH>(cvEl + off);
            // puts("cv loaded");
            jointMeas[vi] = load<SIMD_WIDTH>(jointEl + off);
            // sum of two values in denominator
            measSum[vi] = adds(cvMeas[vi], ssMeas[vi]);
            // inversion for correlation-based joint measures from [-1,1]
            if (invertJointMeasure)
              jointMeas[vi] = subs(measSum[vi], jointMeas[vi]);
          }
          // division and multiplication by postScales
          fdivmul(jointMeas, measSum, postScale, spsVal);
          store(destEl, spsVal, numOutVecs);
        }
      }
    }

  } // CSPNA2_LOOP_AS_IN_CPM2
}

// ---------------------------------------------------------------------------
// computeScalePlaneNormalizedAdd1Tilt
// contributed by Annika Hoffmann
// ---------------------------------------------------------------------------

// normalized distance measure where denominator terms are added
template <int invertJointMeasure, typename SPSType, typename T, int SIMD_WIDTH,
          int SIMD_ALIGN>
void computeScalePlaneNormalizedAdd1Tilt(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &sscvMeasure, int wSPS,
  const int *const modulo2wpw, const int *const shuffle, double postScale,
  SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane)
{
  if (overlap(jointMeasure, scalePlane) || overlap(sscvMeasure, scalePlane))
    throw SIMDException("computeScalePlaneNormalizedAdd1Tilt",
                        "input and output images overlap");
  if (jointMeasure.w != jointMeasure.h)
    throw SIMDException("computeScalePlaneNormalizedAdd1Tilt",
                        "joint measure image is not square-sized");
  if (jointMeasure.w != sscvMeasure.w)
    throw SIMDException("computeScalePlaneNormalizedAdd1Tilt",
                        "joint/ss measure differ in width");

  const int w = jointMeasure.w;
  // number of elements in a vector relating to the input SIMDImages
  const int vecElemsIn = numSIMDVecElements<T, SIMD_WIDTH>();
  // number of input vectors required to fill single output vector
  const int numInVecs = numInputSIMDVecs<SPSType, T>();
  // number of output vectors required to take up single input vector
  const int numOutVecs = numOutputSIMDVecs<SPSType, T>();
  // number of input and output elements
  const int numElems = numSIMDVecsElements<SPSType, T, SIMD_WIDTH>();
  // width must be multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("computeScalePlaneNormalizedAdd1Tilt",
                        "width must be multiple of numElems");

  SIMDVec<T, SIMD_WIDTH> cvMeasxN, joint[numInVecs], sscvMeas[numInVecs];
  SIMDVec<SPSType, SIMD_WIDTH> spsVal[numOutVecs];
  SIMDVec<T, SIMD_WIDTH> ones = setunity<T, SIMD_WIDTH>();

  const int *modulo2woff = modulo2wpw;
  // destination scale plane
  scalePlane.resize(wSPS, w);
  SPSType *destP = scalePlane.data;
  // pointer to dirTerm
  const T *jointPtr    = jointMeasure.data;
  const T *sscvMeasPtr = sscvMeasure.data;
  // i corresponds to theta'
  for (int i = 0; i < w; i++) {
    // we take a single value from cvMeasure and expand it to
    // vecElemsIn identical values (numInVecs vectors);
    const int *shufflepj = shuffle;
    // j corresponds to theta
    for (int j = 0; j < w; j += numElems) {
      for (int vi = 0, off = 0; vi < numInVecs; vi++, off += vecElemsIn) {
        // sum in denominator: we add 1 to avoid div by 0
        // TODO: should we make this epsilon offset variable?
        // 18. Jul 19 (rm): fixed bug reported by Christoph Berganski:
        // adds instead of add to avoid overflow
        sscvMeas[vi] = adds(load<SIMD_WIDTH>(sscvMeasPtr + off), ones);
        // TODO: the above addition of ones is critical if all computations are
        // done in float arithmetic; but does not
        //  effect the results using typical type constellations; with float
        //  arithmetic here, division by zero has to be prevented differently!
        // sscvMeas[vi] = load<SIMD_WIDTH>(sscvMeasPtr +off);
        //  numerator
        joint[vi] = load<SIMD_WIDTH>(jointPtr + off);
        // inversion for correlation-based joint measures from [-1,1]
        if (invertJointMeasure) joint[vi] = subs(sscvMeas[vi], joint[vi]);
      }
      // division and multiplication by postScales
      fdivmul(joint, sscvMeas, postScale, spsVal);
      // store to memory
      SPSType Cstore[numElems] SIMD_ATTR_ALIGNED(SIMD_ALIGN);
      store(Cstore, spsVal, numOutVecs);
      // now we do the "slanted" insertion into dest
      // TODO: store of SPS is quite time-consuming due to inefficient
      // TODO: cache use, especially if single bytes are written
      const int *modulo2woffmj = modulo2woff - j;
      // we have to go through all elements separately
      for (int k = 0; k < numElems; k++, modulo2woffmj--) {
        *(destP + shufflepj[k] + *modulo2woffmj) = Cstore[k];
      }
      jointPtr += numElems;
      sscvMeasPtr += numElems;
      shufflepj += numElems;
    }
    modulo2woff++;
  }
}

// ---------------------------------------------------------------------------
// vertical edge filter
// ---------------------------------------------------------------------------

// vertical edge filter in columns
// Panorama -> Panorama
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class OutAddOn>
void verticalEdgeFilter(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &edge)
{
  // check overlap of arguments
  if (overlap(orig, edge))
    throw SIMDException("verticalEdgeFilter", "input and output image overlap");
  // pointers to original (input) image
  const T *po0, *po1;
  // pointers to edge (output) image
  T *pe0, *pe1;
  // SIMD registers for two pixels and their difference
  SIMDVec<T, SIMD_WIDTH> v0, v1, d;
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  const int w = orig.w, h = orig.h;
  // TODO: edge-filtering removes one row, but
  // TODO: which? isn't that a matter of interpretation? I think
  // TODO: the reason was that in the robot experiments the horizon
  // TODO: can be at the bottom row of the image, and it should not
  // TODO: move outside of the image
  edge.resize(w, h - 1,
              Panorama(orig.addOn.verticalResolution, orig.addOn.horizon - 1));
  // we operate with pointer pairs: *0 always points to block of
  // simd_elems columns in the top row, *1 traverses the columns
  po0 = orig.data;
  pe0 = edge.data;
  // we go through all columns in jumps of simd_elems columns
  for (int x = 0; x < w; x += simd_elems) {
    // we start at the top of the simd_elems-column block
    po1 = po0;
    pe1 = pe0;
    // load values of first row
    v0 = load<SIMD_WIDTH>(po1);
    // go to next row
    po1 += w;
    // go through remaining rows
    for (int y = 1; y < h; y++) {
      // load values of next row
      v1 = load<SIMD_WIDTH>(po1);
      // compute difference (signed types only!)
      d = sub(v1, v0);
      // store difference
      store(pe1, d);
      // keep values of this row for next loop
      v0 = v1;
      // go to next row
      po1 += w;
      pe1 += w;
    }
    // skip simd_elems columns (within top row)
    po0 += simd_elems;
    pe0 += simd_elems;
  }
}

// ---------------------------------------------------------------------------
// vertical edge filter tilt - consider invalid pixels
// contributed by C. Berganski
// ---------------------------------------------------------------------------

// vertical edge filter in columns
// Panorama -> Panorama
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class OutAddOn>
void verticalEdgeFilterTilt(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &edge)
{
  // check overlap of arguments
  if (overlap(orig, edge))
    throw SIMDException("verticalEdgeFilterTilt",
                        "input and output image overlap");

  // pointers to original (input) image
  const T *po0, *po1;
  // pointers to edge (output) image
  T *pe0, *pe1;
  // SIMD registers for two pixels and their difference
  SIMDVec<T, SIMD_WIDTH> v0, v1, d, invvec, invpix;

  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  const int w = orig.w, h = orig.h;
  // TODO: edge-filtering removes one row, but
  // TODO: which? isn't that a matter of interpretation? I think
  // TODO: the reason was that in the robot experiments the horizon
  // TODO: can be at the bottom row of the image, and it should not
  // TODO: move outside of the image
  edge.resize(w, h - 1,
              Panorama(orig.addOn.verticalResolution, orig.addOn.horizon - 1));

  // Prepare invalid vector (vector of type max-value)
  invvec = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::max());

  // we operate with pointer pairs: *0 always points to block of
  // simd_elems columns in the top row, *1 traverses the columns
  po0 = orig.data;
  pe0 = edge.data;
  // we go through all columns in jumps of simd_elems columns
  for (int x = 0; x < w; x += simd_elems) {
    // we start at the top of the simd_elems-column block
    po1 = po0;
    pe1 = pe0;
    // load values of first row
    v0 = load<SIMD_WIDTH>(po1);
    // go to next row
    po1 += w;
    // go through remaining rows
    for (int y = 1; y < h; y++) {
      // load values of next row
      v1 = load<SIMD_WIDTH>(po1);

      // check for invalid pixels (mask); suffices to check the
      // second row as the mask extending to the bottom of
      // the image cannot have invalid pixels in the upper row (v0)
      // but valid in the lower (v1)
      invpix = cmpeq(v1, invvec);

      // compute difference (signed types only!) for valid pixels only
      // extend mask for invalid
      d = ifelse(invpix, invvec, sub(v1, v0));

      // store difference
      store(pe1, d);
      // keep values of this row for next loop
      v0 = v1;
      // go to next row
      po1 += w;
      pe1 += w;
    }
    // skip simd_elems columns (within top row)
    po0 += simd_elems;
    pe0 += simd_elems;
  }
}

// ---------------------------------------------------------------------------
// magnify and scale
// ---------------------------------------------------------------------------

// vertical magnification around horizon with transformation from Tin to Tout
// 27. Feb 18 (rm): verticalResolution and horizon taken from input image,
// not passed as parameters
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class OutAddOn>
void magnifyAndScale(
  const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
  int interpolation, double magnifyScale, double multiplyScale,
  SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &outImg)
{
  // 27. Feb 18 (rm)
  double verticalResolution = img.addOn.verticalResolution;
  double horizon            = img.addOn.horizon;
  double sourceRowDouble;
  const Tin *src;
  Tout *dst;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  const int w = img.w, h = img.h;
  /*
  printf("magnifyImage: %d x %d, v = %.10g, h = %.10g\n",
         w, h, verticalResolution, horizon);
  */
  if ((magnifyScale < 0.0) || (magnifyScale > 1.0))
    throw SIMDException("magnifyAndScale", "magnifyScale out of range [0,1]");
  if (overlap(img, outImg))
    throw SIMDException("magnifyAndScale", "input and ouput image overlap");
  // the horizon must not lie outside the image since otherwise
  // magnification would produce undefined image rows
  // 6. Mar 13 (rm): we decrease the horizon after filtering, so the upper
  // limit was changed from 0.0 to 1.0
  if ((horizon < 0.0) || (horizon > h - 1))
    throw SIMDException("magnifyAndScale", "horizon lies outside the image");
  // w must be a multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("magnifyAndScale",
                        "width must be a multiple of numElems");

  SIMDVec<Tin, SIMD_WIDTH> inVecs1[numInVecs], inVecs2[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];
  outImg.resize(
    w, h,
    Panorama(img.addOn.verticalResolution * magnifyScale, img.addOn.horizon));

  switch (interpolation) {
  case 0:
    // ================ nearest neighbor interpolation ==================
    if (verticalResolution <= 0.0) {
      // negative values or zero indicate: use approximative solution
      // source row index as double
      // 0.5 added for rounding
      // before in the loop:
      // sourceRowDouble = 0.5 + horizon - magnifyScale * (horizon - y)
      sourceRowDouble = 0.5 + horizon - horizon * magnifyScale;
      // copy rows (no interpolation)
      for (int y = 0; y < h; y++) {
        src = img.data + w * static_cast<int>(sourceRowDouble);
        dst = outImg.data + w * y;
        for (int x = 0; x < w; x += numElems) {
          // load numInVecs input vectors
          load(src, inVecs1, numInVecs);
          // process input vectors to output vectors (in float arithmetic)
          fmul(inVecs1, multiplyScale, outVecs);
          // store numOutVec output vectors
          store(dst, outVecs, numOutVecs);
          src += numElems;
          dst += numElems;
        }
        sourceRowDouble += magnifyScale;
      }
    } else {
      // exact solution
      double gamma, gammas;
      double verticalResolutionInv = 1.0 / verticalResolution;
      for (int y = 0; y < h; y++) {
        // verticalResolution: angle[rad] per pixel
        gammas = (horizon - y) * verticalResolution;
        gamma  = atan(magnifyScale * tan(gammas));
        // 0.5 added for rounding
        sourceRowDouble = 0.5 + horizon - gamma * verticalResolutionInv;
        src             = img.data + w * static_cast<int>(sourceRowDouble);
        dst             = outImg.data + w * y;
        for (int x = 0; x < w; x += numElems) {
          // load numInVecs input vectors
          load(src, inVecs1, numInVecs);
          // process input vectors to output vectors (in float arithmetic)
          fmul(inVecs1, multiplyScale, outVecs);
          // store numOutVec output vectors
          store(dst, outVecs, numOutVecs);
          src += numElems;
          dst += numElems;
        }
      }
    }
    break;

  case 1:
    // ================ linear interpolation =================
    double verticalResolutionInv = 1.0 / verticalResolution;
    for (int ys = 0; ys < h; ys++) {
      // source coordinate
      double y;
      if (verticalResolution <= 0.0)
        // negative values indicate: use approximative solution
        y = horizon - magnifyScale * (horizon - ys);
      else {
        // exact solution
        double gammas = (horizon - ys) * verticalResolution;
        double gamma  = atan(magnifyScale * tan(gammas));
        y             = horizon - gamma * verticalResolutionInv;
      }
      int y1 = ::floor(y);
      int y2 = y1 + 1;
      // interpolation weight
      double weight = y2 - y;
      // TODO: better explanation required for the case y2 >= h
      // y2 >= h is only allowed ...
      if (y2 >= h) {
        // ... if weight == 1.0 ...
        if (weight != 1.0)
          throw SIMDException("magnifyAndScale",
                              "y2 >= h is only allowed for weight == 1.0");
        // ... in this case y2 is set to y1
        y2 = y1;
      }
      const Tin *p1, *p2;
      Tout *ps;
      // go throught the two rows and write to result row
      for (p1 = img.data + y1 * w, p2 = img.data + y2 * w,
          ps = outImg.data + ys * w;
           p1 < img.data + y1 * w + w;
           p1 += numElems, p2 += numElems, ps += numElems) {
        // load numInVecs input vectors
        load(p1, inVecs1, numInVecs);
        load(p2, inVecs2, numInVecs);
        // process input vectors
        fwaddmul(inVecs1, inVecs2, weight, multiplyScale, outVecs);
        // store numOutVec output vectors
        store(ps, outVecs, numOutVecs);
      }
    }
    break;
  }
}

// ---------------------------------------------------------------------------
// magnify and scale for classes considering tilt
// after transformation from Tin to Tout, the invalid values for Tin are
// transformed to the invalid values for Tout, not just computed with
// multiplyScale
// contributed by Annika Hoffmann
// ---------------------------------------------------------------------------

// TODO: linear interpolation (interpolation=1) is not implemented yet,
// TODO: is difficult since invalid values should not be considered
// vertical magnification around horizon with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class OutAddOn>
void magnifyAndScaleTilt(
  const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
  int interpolation, double magnifyScale, double multiplyScale,
  SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &outImg)
{
  // only interpolation == 0 is implemented
  assert(interpolation == 0);
  // 27. Feb 18 (rm)
  double verticalResolution = img.addOn.verticalResolution;
  double horizon            = img.addOn.horizon;
  double sourceRowDouble;
  const Tin *src;
  Tout *dst;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  const int w = img.w, h = img.h;

  if ((magnifyScale < 0.0) || (magnifyScale > 1.0))
    throw SIMDException("magnifyAndScaleTilt",
                        "magnifyScale out of range [0,1]");
  if (overlap(img, outImg))
    throw SIMDException("magnifyAndScaleTilt", "input and ouput image overlap");
  // the horizon must not lie outside the image since otherwise
  // magnification would produce undefined image rows
  // 6. Mar 13 (rm): we decrease the horizon after filtering, so the upper
  // limit was changed from 0.0 to 1.0
  if ((horizon < 0.0) || (horizon > h - 1))
    throw SIMDException("magnifyAndScaleTilt",
                        "horizon lies outside the image");
  // w must be a multiple of numElems
  if (w % numElems != 0)
    throw SIMDException("magnifyAndScaleTilt",
                        "width must be a multiple of numElems");
  SIMDVec<Tin, SIMD_WIDTH> inVecs1[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];

  SIMDVec<Tin, SIMD_WIDTH> maxVecIn =
    set1<Tin, SIMD_WIDTH>(SIMDTypeInfo<Tin>::max());
  SIMDVec<Tin, SIMD_WIDTH> cmpVec;
  SIMDVec<Tout, SIMD_WIDTH> maxVecOut =
    set1<Tout, SIMD_WIDTH>(SIMDTypeInfo<Tout>::max());
  SIMDVec<Tout, SIMD_WIDTH> cmpVecOut, cmpVecsOut[numOutVecs];
  SIMDVec<Tin, SIMD_WIDTH> cmpVecs[numInVecs];

  // addon parameters are changed according to magnifyAndScale above
  outImg.resize(
    w, h,
    Panorama(img.addOn.verticalResolution * magnifyScale, img.addOn.horizon));

  // ================ only nearest neighbor interpolation ==================
  if (verticalResolution <= 0.0) {
    // negative values or zero indicate: use approximative solution
    // source row index as double
    // in the loop:
    // sourceRowDouble = 0.5 + horizon - magnifyScale * (horizon - y)
    sourceRowDouble = 0.5 + horizon - horizon * magnifyScale;
    // copy rows (no interpolation)
    for (int y = 0; y < h; y++) {
      src = img.data + w * static_cast<int>(sourceRowDouble);
      dst = outImg.data + w * y;
      for (int x = 0; x < w; x += numElems) {
        // load numInVecs input vectors
        load(src, inVecs1, numInVecs);

        // comparison is only defined for single vectors:
        for (int i = 0; i < numInVecs; i++) {
          cmpVecs[i] = cmpeq(inVecs1[i], maxVecIn);
        }
        convert(cmpVecs, cmpVecsOut);

        // process input vectors to output vectors
        fmul(inVecs1, multiplyScale, outVecs);

        // input maximum values (= invalid) are transformed to output
        // maximum values (= invalid)
        for (int i = 0; i < numOutVecs; i++)
          outVecs[i] = ifelse(cmpVecsOut[i], maxVecOut, outVecs[i]);

        // store numOutVec output vectors
        store(dst, outVecs, numOutVecs);
        src += numElems;
        dst += numElems;
      }
      sourceRowDouble += magnifyScale;
    }
  } else {
    // exact solution
    double gamma, gammas;
    double verticalResolutionInv = 1.0 / verticalResolution;
    for (int y = 0; y < h; y++) {
      // verticalResolution: angle[rad] per pixel
      gammas          = (horizon - y) * verticalResolution;
      gamma           = atan(magnifyScale * tan(gammas));
      sourceRowDouble = 0.5 + horizon - gamma * verticalResolutionInv;
      src             = img.data + w * static_cast<int>(sourceRowDouble);
      dst             = outImg.data + w * y;
      for (int x = 0; x < w; x += numElems) {
        // load numInVecs input vectors
        load(src, inVecs1, numInVecs);
        // comparison is only defined for single vectors:
        for (int i = 0; i < numInVecs; i++) {
          cmpVecs[i] = cmpeq(inVecs1[i], maxVecIn);
        }
        convert(cmpVecs, cmpVecsOut);
        fmul(inVecs1, multiplyScale, outVecs);
        // process input vectors to output vectors
        for (int i = 0; i < numOutVecs; i++)
          outVecs[i] = ifelse(cmpVecsOut[i], maxVecOut, outVecs[i]);

        store(dst, outVecs, numOutVecs);
        src += numElems;
        dst += numElems;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// copy
// ---------------------------------------------------------------------------

// copy with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copy(const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img,
          SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(img, outImg) && overlap(img, outImg))
    throw SIMDException("copy",
                        "input and output image are not identical but overlap");
  const size_t size    = img.size;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  // size must be a multiple of numElems
  if (size % numElems != 0)
    throw SIMDException("copy", "size must be a multiple of numElems");
  SIMDVec<Tin, SIMD_WIDTH> inVecs[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];
  outImg.resize(img.w, img.h, img.addOn);
  const Tin *src = img.data;
  Tout *dst      = outImg.data;
  for (size_t i = 0; i < size; i += numElems) {
    // load numInVecs input vectors
    load(src, inVecs, numInVecs);
    // process input vectors to output vectors
    convert(inVecs, outVecs);
    // store numOutVec output vectors
    store(dst, outVecs, numOutVecs);
    src += numElems;
    dst += numElems;
  }
}

// ---------------------------------------------------------------------------
// copy and scale
// ---------------------------------------------------------------------------

// copy and scale with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copyAndScale(const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img,
                  double multiplyScale,
                  SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(img, outImg) && overlap(img, outImg))
    throw SIMDException("copyAndScale",
                        "input and output image are not identical but overlap");
  const size_t size    = img.size;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  // size must be a multiple of numElems
  if (size % numElems != 0)
    throw SIMDException("copyAndScale", "size must be a multiple of numElems");
  SIMDVec<Tin, SIMD_WIDTH> inVecs[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];
  outImg.resize(img.w, img.h, img.addOn);
  const Tin *src = img.data;
  Tout *dst      = outImg.data;
  for (size_t i = 0; i < size; i += numElems) {
    // load numInVecs input vectors
    load(src, inVecs, numInVecs);
    // process input vectors to output vectors (in float arithmetic)
    fmul(inVecs, multiplyScale, outVecs);
    // store numOutVec output vectors
    store(dst, outVecs, numOutVecs);
    src += numElems;
    dst += numElems;
  }
}

// each plane in stack has its own multiplyScale, so
// mulitplyScale[stk->numPlanes] is assumed
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copyAndScale(const SIMDStack<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &stk,
                  const std::vector<double> &multiplyScale,
                  SIMDStack<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outStk)
{
  if (multiplyScale.size() != (size_t) stk.numPlanes)
    throw SIMDException("copyAndScale", "multiplyScale vector has wrong size");
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(stk, outStk) && overlap(stk, outStk))
    throw SIMDException("copyAndScale",
                        "input and output stack are not identical but overlap");
  outStk.resize(stk);
  for (int i = 0; i < stk.numPlanes; i++)
    copyAndScale(stk[i], multiplyScale[i], outStk[i]);
}

// ---------------------------------------------------------------------------
// copy and scale for tilted images (contributed by Annika Hoffmann)
// ---------------------------------------------------------------------------

// copy and scale with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copyAndScaleTilt(
  const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img,
  double multiplyScale,
  SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(img, outImg) && overlap(img, outImg))
    throw SIMDException("copyAndScaleTilt",
                        "input and output image are not identical but overlap");
  const size_t size    = img.size;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  // size must be a multiple of numElems
  if (size % numElems != 0)
    throw SIMDException("copyAndScaleTilt",
                        "size must be a multiple of numElems");

  SIMDVec<Tin, SIMD_WIDTH> maxVecIn =
    set1<Tin, SIMD_WIDTH>(SIMDTypeInfo<Tin>::max());
  SIMDVec<Tin, SIMD_WIDTH> cmpVec;
  SIMDVec<Tout, SIMD_WIDTH> maxVecOut =
    set1<Tout, SIMD_WIDTH>(SIMDTypeInfo<Tout>::max());
  SIMDVec<Tout, SIMD_WIDTH> cmpVecOut, cmpVecsOut[numOutVecs];
  SIMDVec<Tin, SIMD_WIDTH> cmpVecs[numInVecs];
  SIMDVec<Tin, SIMD_WIDTH> inVecs[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];

  outImg.resize(img.w, img.h, img.addOn);
  const Tin *src = img.data;
  Tout *dst      = outImg.data;
  for (size_t i = 0; i < size; i += numElems) {
    // load numInVecs input vectors
    load(src, inVecs, numInVecs);
    // comparison is only defined for single vectors:
    //  1. Dec 21 (rm): i -> j (shadows outer j)
    for (int j = 0; j < numInVecs; j++) {
      cmpVecs[j] = cmpeq(inVecs[j], maxVecIn);
    }
    convert(cmpVecs, cmpVecsOut);
    // process input vectors to output vectors
    fmul(inVecs, multiplyScale, outVecs);
    // set the invalid pixels
    //  1. Dec 21 (rm): i -> j (shadows outer j)
    for (int j = 0; j < numOutVecs; j++)
      outVecs[j] = ifelse(cmpVecsOut[j], maxVecOut, outVecs[j]);
    // store numOutVec output vectors
    store(dst, outVecs, numOutVecs);
    src += numElems;
    dst += numElems;
  }
}

// ---------------------------------------------------------------------------
// shift and scale
// ---------------------------------------------------------------------------

// shift and scale with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copyShiftAndScale(
  const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img, double shift,
  double scale, SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(img, outImg) && overlap(img, outImg))
    throw SIMDException("copyShiftAndScale",
                        "input and output image are not identical but overlap");
  const size_t size    = img.size;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  // size must be a multiple of numElems
  if (size % numElems != 0)
    throw SIMDException("copyShiftAndScale",
                        "size must be a multiple of numElems");
  SIMDVec<Tin, SIMD_WIDTH> inVecs[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];
  outImg.resize(img.w, img.h, img.addOn);
  const Tin *src = img.data;
  Tout *dst      = outImg.data;
  for (size_t i = 0; i < size; i += numElems) {
    // load numInVecs input vectors
    load(src, inVecs, numInVecs);
    // process input vectors to output vectors (in float arithmetic)
    faddmul(inVecs, shift, scale, outVecs);
    // store numOutVec output vectors
    store(dst, outVecs, numOutVecs);
    src += numElems;
    dst += numElems;
  }
}

// ---------------------------------------------------------------------------
// scale and shift
// ---------------------------------------------------------------------------

// shift and scale with transformation from Tin to Tout
template <typename Tout, typename Tin, int SIMD_WIDTH, int SIMD_ALIGN,
          class AnyAddOn1, class AnyAddOn2>
void copyScaleAndShift(
  const SIMDImage<Tin, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &img, double scale,
  double shift, SIMDImage<Tout, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(img, outImg) && overlap(img, outImg))
    throw SIMDException("copyScaleAndShift",
                        "input and output image are not identical but overlap");
  const size_t size    = img.size;
  const int numInVecs  = numInputSIMDVecs<Tout, Tin>();
  const int numOutVecs = numOutputSIMDVecs<Tout, Tin>();
  const int numElems   = numSIMDVecsElements<Tout, Tin, SIMD_WIDTH>();
  // size must be a multiple of numElems
  if (size % numElems != 0)
    throw SIMDException("copyScaleAndShift",
                        "size must be a multiple of numElems");
  SIMDVec<Tin, SIMD_WIDTH> inVecs[numInVecs];
  SIMDVec<Tout, SIMD_WIDTH> outVecs[numOutVecs];
  outImg.resize(img.w, img.h, img.addOn);
  const Tin *src = img.data;
  Tout *dst      = outImg.data;
  for (size_t i = 0; i < size; i += numElems) {
    // load numInVecs input vectors
    load(src, inVecs, numInVecs);
    // process input vectors to output vectors (in float arithmetic)
    fmuladd(inVecs, scale, shift, outVecs);
    // store numOutVec output vectors
    store(dst, outVecs, numOutVecs);
    src += numElems;
    dst += numElems;
  }
}

// ---------------------------------------------------------------------------
// insert
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnySrcAddOn,
          class AnyDstAddOn>
void insert(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnySrcAddOn> &src,
            int x0, int y0,
            SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyDstAddOn> &dst)
{
  if ((x0 + src.w > dst.w) || (y0 + src.h > dst.h))
    throw SIMDException("insert", "size of destination image exceeded");
  for (int yd = y0, ys = 0; ys < src.h; yd++, ys++)
    memcpy(dst[yd] + x0, src[ys], src.w * sizeof(T));
}

// ---------------------------------------------------------------------------
// minimum / maximum
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
T max(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  SIMDVec<T, SIMD_WIDTH> maxv = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::min());
  const int nElems            = numSIMDVecElements<T, SIMD_WIDTH>();
  for (T *d = img.data; d < img.data + img.size; d += nElems)
    maxv = max(maxv, load<SIMD_WIDTH>(d));
  return hmax(maxv);
}

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
T maxExcept(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img,
            T exceptVal)
{
  SIMDVec<T, SIMD_WIDTH> maxv = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::min()),
                         excv = set1<T, SIMD_WIDTH>(exceptVal), v, isexcv;
  const int nElems            = numSIMDVecElements<T, SIMD_WIDTH>();
  for (T *d = img.data; d < img.data + img.size; d += nElems) {
    v      = load<SIMD_WIDTH>(d);
    isexcv = cmpeq(v, excv);
    maxv   = ifelse(isexcv, maxv, max(maxv, v));
  }
  T maxvbuf[nElems];
  storeu(maxvbuf, maxv);
  T maxs = maxvbuf[0];
  for (int i = 1; i < nElems; i++) maxs = std::max(maxs, maxvbuf[i]);
  return maxs;
}

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
T min(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  SIMDVec<T, SIMD_WIDTH> minv = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::max());
  const int nElems            = numSIMDVecElements<T, SIMD_WIDTH>();
  for (T *d = img.data; d < img.data + img.size; d += nElems)
    minv = min(minv, load<SIMD_WIDTH>(d));
  return hmin(minv);
}

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
T maxAbs(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  SIMDVec<T, SIMD_WIDTH> maxv = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::min());
  const int nElems            = numSIMDVecElements<T, SIMD_WIDTH>();
  for (T *d = img.data; d < img.data + img.size; d += nElems)
    maxv = max(maxv, abs(load<SIMD_WIDTH>(d)));
  return hmax(maxv);
}

// ---------------------------------------------------------------------------
// replace pixels oldVal with newVal
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
void replace(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img,
             T oldVal, T newVal,
             SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &res)
{
  SIMDVec<T, SIMD_WIDTH> oldv = set1<T, SIMD_WIDTH>(oldVal),
                         newv = set1<T, SIMD_WIDTH>(newVal), v;
  res                         = img;
  const int nElems            = numSIMDVecElements<T, SIMD_WIDTH>();
  for (T *d = res.data; d < res.data + res.size; d += nElems) {
    v = load<SIMD_WIDTH>(d);
    store(d, ifelse(cmpeq(oldv, v), newv, v));
  }
}

// ---------------------------------------------------------------------------
// load and save pgm files
// ---------------------------------------------------------------------------

// load from FILE*
// cannot handle comments in PGM headers (complicated)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
bool loadPGM(FILE *f, SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AddOn> &img,
             const AddOn &addOn = AddOn())
{
  int w, h, mx, ret;
  // do not write \n at the end of the format specification,
  // will also read other special chars after the \n!
  // "A sequence of white-space characters (space, tab, newline, etc.;
  // see isspace(3)).  This directive matches  any  amount  of  white
  // space, including none, in the input."
  ret = fscanf(f, "P5 %d %d %d", &w, &h, &mx);
  if (ret != 3) return false;
  if (mx > SIMDBYTE_MAX) return false;
  // instead we read a single byte to skip the single whitespace
  ret = fgetc(f);
  if (ret == EOF) return false;
  // create SIMDByte image
  // SIMD_WIDTH=SIMD_ALIGN=1 so that each image fits (TODO: necessary here?)
  SIMDImage<SIMDByte, 1, 1, AddOn> imgB;
  imgB.resize(w, h, addOn);
  // read data
  if (!imgB.load(f)) return false;
  // assign to argument image
  img = imgB;
  return true;
}

// load from filename (char*)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
bool loadPGM(const char *filename,
             SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AddOn> &img,
             const AddOn &addOn = AddOn())
{
  FILE *f = fopen(filename, "r");
  if (!f) return false;
  bool res = loadPGM(f, img, addOn);
  fclose(f);
  return res;
}

// load from filename (std::string)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
bool loadPGM(const std::string &filename,
             SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AddOn> &img,
             const AddOn &addOn = AddOn())
{
  return loadPGM(filename.c_str(), img, addOn);
}

// save to FILE*
// the user has to guarantee that the pixels of img fit into a SIMDByte
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
bool savePGM(FILE *f,
             const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  int ret;
  ret = fprintf(f, "P5\n%d %d\n%d\n", img.w, img.h, SIMDBYTE_MAX);
  if (ret < 0) return false;
  // note that this calls the copy constructor, not operator=
  // SIMD_WIDTH=SIMD_ALIGN=1 so that each image fits
  SIMDImage<SIMDByte, 1, 1> imgB = img;
  return imgB.save(f);
}

// save to filename (char*)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
bool savePGM(const char *filename,
             const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  FILE *f = fopen(filename, "w");
  if (!f) return false;
  bool res = savePGM(f, img);
  fclose(f);
  return res;
}

// save to filename (std::string)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn>
bool savePGM(const std::string &filename,
             const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn> &img)
{
  return savePGM(filename.c_str(), img);
}

// ---------------------------------------------------------------------------
// rotateHor: rotate horizontally (cyclic)
// ---------------------------------------------------------------------------

template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void rotateHor(
  const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg, int shift,
  SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  if (overlap(inImg, outImg))
    throw SIMDException("rotateHor", "input and output image overlap");
  const int w = inImg.w, h = inImg.h;
  if ((shift < 0) || (shift >= w))
    throw SIMDException("rotateHor", "shift parameter must be in [0,w)");
  inImg.addOn.isInvariantToHorizontalReorganization();
  outImg.resize(w, h, inImg.addOn);
  ImgType *ip = inImg.data, *op = outImg.data;
  const int wms = w - shift;
  for (int j = 0; j < h; j++, ip += w, op += w) {
    memcpy(op, ip + shift, wms * sizeof(ImgType));
    memcpy(op + wms, ip, shift * sizeof(ImgType));
  }
}

// ---------------------------------------------------------------------------
// swizzle (cyclic)
// ---------------------------------------------------------------------------

template <int N, typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void swizzleCyclic(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
                   int nAuxVecCols,
                   SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  if (overlap(inImg, outImg))
    throw SIMDException("swizzleCyclic", "input and output image overlap");
  const int w = inImg.w, h = inImg.h;
  const int elems = numSIMDVecElements<T, SIMD_WIDTH>();
  const int n = N, n_x_elems = n * elems, blockSize = w / n;
  const int blockSizePlus = blockSize + nAuxVecCols * elems;
  const int wPlus         = w + nAuxVecCols * n_x_elems;
  // printf("w = %d, h = %d, elems = %d, n = %d, "
  //        "n_x_elems = %d, blockSize = %d\n",
  //	    w, h, elems, n, n_x_elems, blockSize);
  if (w % n_x_elems != 0)
    throw SIMDException("swizzleCyclic", "invalid width (w % n_x_elems != 0)");
  SIMDVec<T, SIMD_WIDTH> vecs[n];
  // outImg contains n_x_elems additional columns
  inImg.addOn.isInvariantToHorizontalReorganization();
  outImg.resize(wPlus, h, inImg.addOn);
  // blockOff is the offset to the start of each block
  int blockOff[n];
  for (int k = 0, off = 0; k < n; k++, off += blockSizePlus) blockOff[k] = off;
  // go through all rows
  for (int j = 0; j < h; j++) {
    // start of rows in inImg and outImg
    T *idp = inImg[j], *outRow = outImg[j], *odp = outRow;
    // go along this row
    for (int i = 0; i < w; i += n_x_elems, idp += n_x_elems, odp += elems) {
      // load n vectors
      load(idp, vecs, n);
      // swizzle vectors
      swizzle<N>(vecs);
      // store n vectors to different blocks (with same offset
      // relative to block, contained in odp)
      for (int k = 0; k < n; k++) store(odp + blockOff[k], vecs[k]);
    }
    // complement wrap-around columns in each block
    for (int k = 0; k < n; k++) {
      T *obp = outRow + blockOff[k];
      for (int m = 0; m < nAuxVecCols; m++, obp += elems)
        load_store<SIMD_WIDTH>(obp, obp + blockSize);
    }
  }
}

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void swizzleCyclic(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
                   int n, int nAuxVecCols,
                   SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  switch (n) {
  case 1: swizzleCyclic<1>(inImg, nAuxVecCols, outImg); break;
  case 2: swizzleCyclic<2>(inImg, nAuxVecCols, outImg); break;
  case 3: swizzleCyclic<3>(inImg, nAuxVecCols, outImg); break;
  case 4: swizzleCyclic<4>(inImg, nAuxVecCols, outImg); break;
  case 5: swizzleCyclic<5>(inImg, nAuxVecCols, outImg); break;

  default:
    throw SIMDException("swizzleCyclic", "invalid n (only 1..5 supported)");
  };
}

// ---------------------------------------------------------------------------
// unSwizzleCyclic (serial, slow)
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void unSwizzleCyclic(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg, int n,
  int nAuxVecCols, SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  const int wPlus = inImg.w, h = inImg.h;
  const int elems     = numSIMDVecElements<T, SIMD_WIDTH>();
  const int n_x_elems = n * elems, blockSizePlus = wPlus / n;
  // const int blockSize = blockSizePlus - nAuxVecCols * elems;
  const int w = wPlus - nAuxVecCols * n_x_elems;
  // printf("w = %d, h = %d, elems = %d, n = %d, "
  //        "n_x_elems = %d, blockSize = %d, blockSizePlus = %d, wPlus = %d\n",
  //	     w, h, elems, n, n_x_elems, blockSize, blockSizePlus, wPlus);
  outImg.resize(w, h, inImg.addOn);
  // col. index output image
  for (int io = 0; io < w; io++) {
    // col. index input image
    int ii = (io % n) * blockSizePlus + ::floor(io / n);
    for (int j = 0; j < h; j++) outImg[j][io] = inImg[j][ii];
  }
}

// ---------------------------------------------------------------------------
// numeric operations
// ---------------------------------------------------------------------------

// saturated add of a constant
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void adds(const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
          ImgType val,
          SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(inImg, outImg) && overlap(inImg, outImg))
    throw SIMDException("adds",
                        "input and output image are not identical but overlap");
  SIMDVec<ImgType, SIMD_WIDTH> valVec = set1<ImgType, SIMD_WIDTH>(val);
  outImg                              = inImg;
  const int elems = numSIMDVecElements<ImgType, SIMD_WIDTH>();
  for (ImgType *dp = outImg.data; dp < outImg.data + outImg.size; dp += elems)
    store(dp, adds(load<SIMD_WIDTH>(dp), valVec));
}

// saturated sub of a constant
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void subs(const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
          ImgType val,
          SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(inImg, outImg) && overlap(inImg, outImg))
    throw SIMDException("subs",
                        "input and output image are not identical but overlap");
  SIMDVec<ImgType, SIMD_WIDTH> valVec = set1<ImgType, SIMD_WIDTH>(val);
  outImg                              = inImg;
  const int elems = numSIMDVecElements<ImgType, SIMD_WIDTH>();
  for (ImgType *dp = outImg.data; dp < outImg.data + outImg.size; dp += elems)
    store(dp, subs(load<SIMD_WIDTH>(dp), valVec));
}

// multiplication by a constant (only for ImgType which provides mul)
template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void mul(const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
         ImgType val,
         SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(inImg, outImg) && overlap(inImg, outImg))
    throw SIMDException("mul",
                        "input and output image are not identical but overlap");
  SIMDVec<ImgType, SIMD_WIDTH> valVec = set1<ImgType, SIMD_WIDTH>(val);
  outImg                              = inImg;
  const int elems = numSIMDVecElements<ImgType, SIMD_WIDTH>();
  for (ImgType *dp = outImg.data; dp < outImg.data + outImg.size; dp += elems)
    store(dp, mul(load<SIMD_WIDTH>(dp), valVec));
}

// pixel-wise multiplication (only for ImgType which provides mul)
// contributed by Benedikt Volkmer (modified from operator*)
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
void mul(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &a,
         const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &b,
         SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &c)
{
  if (!a.sameSize(b))
    throw SIMDException("mul", "input images have different size");
  c.resize(a.w, a.h, a.addOn);
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  for (size_t p = 0; p < c.size; p += simd_elems) {
    SIMDVec<T, SIMD_WIDTH> a_v = load<SIMD_WIDTH>(a.data + p);
    SIMDVec<T, SIMD_WIDTH> b_v = load<SIMD_WIDTH>(b.data + p);
    SIMDVec<T, SIMD_WIDTH> c_v = a_v * b_v;
    store(c.data + p, c_v);
  }
}

// ---------------------------------------------------------------------------
// horizontal cyclic binomial filter
// ---------------------------------------------------------------------------

// FUNCTION IS OBSOLETE, USE horizontalBinomialFilterCyclic INSTEAD!

// example: w = 12, elems = 4
// in row:    |  0  1  2  3  |  4  5  6  7  |  8  9 10 11 |
// left/next: | 11  0  1  2  |  3  4  5  6  |  7  8  9 10 | 11  0  1  2
// center:    |  0  1  2  3  |  4  5  6  7  |  8  9 10 11 |
// right:     |  1  2  3  4  |  5  6  7  8  |  9 10 11  0 |

template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN, class AnyAddOn1,
          class AnyAddOn2>
void binomialHorCyclic(
  const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn1> &inImg,
  SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, AnyAddOn2> &outImg)
{
  // input and output images can be identical,
  // otherwise they shouldn't overlap
  if (!same(inImg, outImg) && overlap(inImg, outImg))
    throw SIMDException("binomialHorCyclic",
                        "input and output image are not identical but overlap");
  const int w = inImg.w, h = inImg.h;
  inImg.addOn.isInvariantToHorizontalReorganization();
  outImg.resize(w, h, inImg.addOn);
  const int elems = numSIMDVecElements<ImgType, SIMD_WIDTH>();
  // ImgType rowBuf[w + elems] SIMD_ATTR_ALIGNED(SIMD_WIDTH);
  SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> rowBuf(w + elems, 1);
  const ImgType *inRowPtr = inImg.data;
  ImgType *outRowPtr      = outImg.data;
  SIMDVec<ImgType, SIMD_WIDTH> left, next, center, right;
  for (int j = 0; j < h; j++, inRowPtr += w, outRowPtr += w) {
    // make a wrap-around copy of the row:
    // - first element is last element
    rowBuf.data[0] = inRowPtr[w - 1];
    // - then come the remaining elements (w - 1 would be sufficient,
    //   but maybe w helps memcpy() to find an efficient vectorized
    //   solution?)
    memcpy(rowBuf.data + 1, inRowPtr, w * sizeof(ImgType));
    // - the first elements (-1, 0, 1, ...) go to the end
    load_store<SIMD_WIDTH>(rowBuf.data, rowBuf.data + w);
    // process
    ImgType *bufPtr = rowBuf.data, *outPtr = outRowPtr;
    // left pixels
    left = load<SIMD_WIDTH>(bufPtr);
    for (int i = 0; i < w; i += elems, outPtr += elems) {
      // load next vector of pixels
      bufPtr += elems;
      next = load<SIMD_WIDTH>(bufPtr);
      // center pixels
      center = alignre<1>(next, left);
      // right pixels
      right = alignre<2>(next, left);
      // filter
      store(outPtr, avg(avg(left, right), center));
      // next round
      left = next;
    }
  }
}

// ---------------------------------------------------------------------------
// combineImageWithMask
// contributed by Annika Hoffmann
// ---------------------------------------------------------------------------

// use a mask image to set all "invalid" pixels to the maximum value
// of the current type
// imageIn and imageOut are in IMGType format (SIMDFloat)
// imgMask is in SIMDFloat format
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AddOnIn,
          class AddOnOut>
void combineImageWithMask(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOnIn> &imgIn,
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOnIn> &imgMask,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOnOut> &imgOut)
{
  const int w = imgIn.w, h = imgIn.h;
  imgOut.resize(w, h, imgIn.addOn);

  // vector for result computation
  SIMDVec<T, SIMD_WIDTH> resVec;
  SIMDVec<T, SIMD_WIDTH> maskVec;
  SIMDVec<T, SIMD_WIDTH> inVec;
  SIMDVec<T, SIMD_WIDTH> maxVec;

  maxVec               = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::max());
  const int simd_elems = SIMD_WIDTH / sizeof(T);

  // pointer to input image: we work with two pointers, one (pi0)
  // pointing to the top row, the other (pi1) traversing the column block
  // pi0 traverses top row in jumps of simd_elems values
  const T *pi0 = imgIn.data;
  // pointer to mask image
  const T *pmask0 = imgMask.data;

  // pas traverses result (single row) in jumps of simd_elems values
  T *presult = imgOut.data;
  // 1 needs to be accepted as true (does not work for Byte-images with 255 as
  // maximal value) could the mask be something different than float image?
  SIMDVec<T, SIMD_WIDTH> cmpVec = set1<T, SIMD_WIDTH>(1);
  SIMDVec<T, SIMD_WIDTH> testVec;

  // go through top column in jumps of simd_elems
  for (int x = 0; x < w; x += simd_elems) {
    // pi1 traverses block of 8 columns, starting at the top row
    const T *pi1    = pi0;
    const T *pmask1 = pmask0;
    T *presult1     = presult;
    // go through column block
    for (int y = 0; y < h; y++) {
      inVec   = load<SIMD_WIDTH>(pi1);
      maskVec = load<SIMD_WIDTH>(pmask1);
      testVec = cmpge(maskVec, cmpVec);
      resVec  = ifelse(testVec, inVec, maxVec);
      pi1 += w;
      pmask1 += w;
      // store the result vector in the output image
      store(presult1, resVec);
      presult1 += w;
    }
    // next block of 8 columns
    pi0 += simd_elems;
    pmask0 += simd_elems;
    presult += simd_elems;
  }
}

// ---------------------------------------------------------------------------
// sumOfAbsAccuInvalid
// contributed by Annika Hoffmann
// ---------------------------------------------------------------------------

// sum of absolute values in all columns (accumulated version)
// row 1: sum of absolute values of first row,
// row h: sum of absolute values of all rows
// absSum has to be of size w x h
template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn1,
          class AddOn2>
void sumOfAbsAccuInvalid(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn1> &img,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn2> &absSum)
{
  SIMDVec<T, SIMD_WIDTH> as;
  // create vector with invalid pixels (which are always the max of
  // the used type
  SIMDVec<T, SIMD_WIDTH> invalidVec, compareVec, summandVec, zero, inVec;
  zero       = setzero<T, SIMD_WIDTH>();
  invalidVec = set1<T, SIMD_WIDTH>(SIMDTypeInfo<T>::max());

  // number of elements in vector
  const int simd_elems = SIMD_WIDTH / sizeof(T);
  const int w = img.w, h = img.h;
  // 2-dimensional image with all column sums, starting with row 1
  // up to row height
  absSum.resize(w, h, img.addOn);
  // temporary absSum, summed up from row 1 up to y
  // pointer to input image: we work with two pointers, one (pi0)
  // pointing to the top row, the other (pi1) traversing the column block
  // pi0 traverses top row in jumps of simd_elems values
  const T *pi0 = img.data;
  // pas traverses result  in jumps of simd_elems values
  T *pas = absSum.data;

  // go through top column in jumps of simd_elems
  for (int x = 0; x < w; x += simd_elems) {
    // pi1 traverses block of 8 columns, starting at the top row
    const T *pi1 = pi0;
    T *pas1      = pas + x;
    // sum of absolute values = 0
    as = setzero<T, SIMD_WIDTH>();
    // go through column block
    for (int y = 0; y < h; y++) {
      inVec = load<SIMD_WIDTH>(pi1);
      // add absolute values with saturation
      compareVec = cmpeq(inVec, invalidVec);
      summandVec = ifelse(compareVec, zero, inVec);
      as         = adds(as, abs(summandVec));
      // go to next row
      // store absolute sum for each additional pixel
      store(pas1, as);
      pi1 += w;
      pas1 += w;
    }
    // next block of 8 columns
    pi0 += simd_elems;
  }
}

// ---------------------------------------------------------------------------
// vertical edge filter [+1,0,-1]
// contributed by Benedikt Volkmer
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class InAddOn,
          class OutAddOn>
void verticalEdgeFilter2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &edge)
{
  // overlap does not interfere
  ////check overlap
  // if (overlap(orig, edge))
  //     throw SIMDException("verticalEdgeFilter2",
  //             "input and output image overlap");
  const int w = orig.w, h = orig.h;
  edge.resize(w, h, orig.addOn);
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  SIMDVec<T, SIMD_WIDTH> top_v, bottom_v, lastBottom_v, diff_v;
  for (int p = 0; p < w; p += simd_elems) {
    // special handling of first row [+1,-1]
    top_v    = load<SIMD_WIDTH>(orig.data + p);
    bottom_v = load<SIMD_WIDTH>(orig.data + w + p);
    // difference of first and second row
    diff_v = sub(top_v, bottom_v);
    // store in first row
    store(edge.data + p, diff_v);
    // top = top
    lastBottom_v = bottom_v;
    // traverse rows except last
    for (int row = 1; row < (h - 1); row++) {
      // load new bottom vectors
      bottom_v = load<SIMD_WIDTH>(orig.data + (row + 1) * w + p);
      // difference between top and bottom
      // (top*1+current*0+bottom*(-1)=top-bottom)
      diff_v = sub(top_v, bottom_v);
      // store to edge
      store(edge.data + row * w + p, diff_v);
      // initialise next loop
      top_v        = lastBottom_v;
      lastBottom_v = bottom_v;
    }
    // special handling of last row [+1,-1]
    // difference of first and second row
    diff_v = sub(top_v, bottom_v);
    // store in last row
    store(edge.data + (w * (h - 1)) + p, diff_v);
  }
}

// ---------------------------------------------------------------------------
// vertical binomial filter [1/4,1/2,1/4]
// contributed by Benedikt Volkmer
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class InAddOn,
          class OutAddOn>
void verticalBinomialFilter(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &binom)
{
  // Overlap does not interfere
  ////check overlap
  // if (overlap(orig, binom))
  //     throw SIMDException("verticalBinomialFilter",
  //             "input and output image overlap");
  const int w = orig.w, h = orig.h;
  binom.resize(w, h, orig.addOn);
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  SIMDVec<T, SIMD_WIDTH> top_v, center_v, bottom_v, avg_v;
  for (int p = 0; p < w; p += simd_elems) {
    // special handling of first row [1/2,1/2], i.e. mirroring second row
    top_v    = load<SIMD_WIDTH>(orig.data + p);
    center_v = load<SIMD_WIDTH>(orig.data + w + p);
    // average of first and second row
    avg_v = avg(top_v, center_v);
    // store in first row
    store(binom.data + p, avg_v);
    // top = top
    // center = center
    // traverse rows except last
    for (int row = 1; row < (h - 1); row++) {
      // load new bottom vectors
      bottom_v = load<SIMD_WIDTH>(orig.data + (row + 1) * w + p);
      // average of center and average between top and bottom
      avg_v = avg(avg(top_v, bottom_v), center_v);
      // store to binom
      store(binom.data + row * w + p, avg_v);
      // initialise next loop
      top_v    = center_v;
      center_v = bottom_v;
    }
    // special handling of last row [1/2,1/2]
    // average of first and second row
    avg_v = avg(top_v, bottom_v);
    // store in last row
    store(binom.data + (w * (h - 1)) + p, avg_v);
  }
}

// ---------------------------------------------------------------------------
// horizontal edge filter [1,0,-1]
// contributed by Benedikt Volkmer
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class InAddOn,
          class OutAddOn>
void horizontalEdgeFilterCyclic2(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &edge)
{
  // overlap does not interfere
  ////check overlap
  // if (overlap(orig, edge))
  //     throw SIMDException("horizontalEdgeFilterCyclic",
  //             "input and output image overlap");
  const int w = orig.w, h = orig.h;
  edge.resize(w, h, orig.addOn);
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  for (size_t row = 0; row < edge.size; row += w) {
    SIMDVec<T, SIMD_WIDTH> next_v, left_v, right_v, previousResult_v,
      tmpResult_v, result_v;
    // calculate last vector tmp result, which need first elements shifted in
    next_v  = load<SIMD_WIDTH>(orig.data + row);
    left_v  = load<SIMD_WIDTH>(orig.data + row + w - simd_elems);
    right_v = alignre<2>(next_v, left_v);
    // l l l l l l l l|n n n n n n n n
    //     r r r r r r r r
    const SIMDVec<T, SIMD_WIDTH> lastResult_v = sub(left_v, right_v);
    //   l l l l l l l l
    // - l l l l l l n n
    previousResult_v = lastResult_v;
    left_v           = next_v;
    // calculate all vectors except the last
    // see horizontalBinomialFilterCyclic for analog explanation
    for (int p = 0; p < w - simd_elems; p += simd_elems) {
      next_v      = load<SIMD_WIDTH>(orig.data + row + p + simd_elems);
      right_v     = alignre<2>(next_v, left_v);
      tmpResult_v = sub(left_v, right_v);
      store(edge.data + row + p,
            alignre<simd_elems - 1>(tmpResult_v, previousResult_v));
      previousResult_v = tmpResult_v;
      left_v           = next_v;
    }
    // store last result
    store(edge.data + row + w - simd_elems,
          alignre<simd_elems - 1>(lastResult_v, previousResult_v));
  }
}

// ---------------------------------------------------------------------------
// horizontal binomial filter [1/4,1/2,1/4]
// contributed by Benedikt Volkmer
// ---------------------------------------------------------------------------

template <typename T, int SIMD_WIDTH, int SIMD_ALIGN, class InAddOn,
          class OutAddOn>
void horizontalBinomialFilterCyclic(
  const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &edge)
{
  // overlap does not interfere
  ////check overlap
  // if (overlap(orig, edge))
  //     throw SIMDException("horizontalBinomialFilterCyclic",
  //             "input and output image overlap");
  const int w = orig.w, h = orig.h;
  edge.resize(w, h, orig.addOn);
  const int simd_elems = numSIMDVecElements<T, SIMD_WIDTH>();
  for (size_t row = 0; row < edge.size; row += w) {
    SIMDVec<T, SIMD_WIDTH> next_v, left_v, right_v, center_v, previousResult_v,
      tmpResult_v, result_v;
    // calculate last vector tmp result, which need first elements shifted in
    next_v   = load<SIMD_WIDTH>(orig.data + row);
    left_v   = load<SIMD_WIDTH>(orig.data + row + w - simd_elems);
    center_v = alignre<1>(next_v, left_v);
    right_v  = alignre<2>(next_v, left_v);
    const SIMDVec<T, SIMD_WIDTH> lastResult_v =
      avg(avg(left_v, right_v), center_v);
    previousResult_v = lastResult_v;
    left_v           = next_v;
    // calculate all vectors except the last
    // the tmpResult is the result shifted right by 1
    // shifting it into the previous tmpResult by vector length -1 will be
    // the right result
    // example with vector length 4 applied to indices [0,1,2,3,4,5,6,7]
    // p=0: left_v = [0,1,2,3] (instead of [7,0,1,2]
    //     right_v = [2,3,4,5] (instead of [1,2,3,4]
    //    center_v = [1,2,3,4] (instead of [0,1,2,3]
    //   tmpResult = [1,2,3,4]
    //      result = [5,6,7,0] <- [1,2,3] = [0,1,2,3]
    // The last also shows, why lastResult has to be calculated first
    // this is necessary because there is just alignr not alignl
    for (int p = 0; p < w - simd_elems; p += simd_elems) {
      next_v      = load<SIMD_WIDTH>(orig.data + row + p + simd_elems);
      center_v    = alignre<1>(next_v, left_v);
      right_v     = alignre<2>(next_v, left_v);
      tmpResult_v = avg(avg(left_v, right_v), center_v);
      store(edge.data + row + p,
            alignre<simd_elems - 1>(tmpResult_v, previousResult_v));
      previousResult_v = tmpResult_v;
      left_v           = next_v;
    }
    // store last result
    store(edge.data + row + w - simd_elems,
          alignre<simd_elems - 1>(lastResult_v, previousResult_v));
  }
}

// ---------------------------------------------------------------------------
// harrisFilterCyclic
// contributed by Benedikt Volkmer (modified: SIMDFloat only)
// ---------------------------------------------------------------------------

template <int NORMALIZE, bool VARIANT, int SIMD_WIDTH, int SIMD_ALIGN,
          class InAddOn, class OutAddOn>
void harrisFilterCyclic(
  const SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &corner,
  const double k, const int binomialFilterCount)
{
  SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN> X, Y, A, B, C;
  verticalEdgeFilter2(orig, X);
  horizontalEdgeFilterCyclic2(orig, Y);
  mul(X, X, A);
  mul(Y, Y, B);
  mul(X, Y, C);
  // TODO: more cache-friendly with 3 loops?
  for (int i = 0; i < binomialFilterCount; ++i) {
    horizontalBinomialFilterCyclic(A, A);
    verticalBinomialFilter(A, A);
    horizontalBinomialFilterCyclic(B, B);
    verticalBinomialFilter(B, B);
    horizontalBinomialFilterCyclic(C, C);
    verticalBinomialFilter(C, C);
  }
  corner.resize(orig.w, orig.h, orig.addOn);
  const int simd_elems = numSIMDVecElements<SIMDFloat, SIMD_WIDTH>();
  // for variant:
  SIMDVec<SIMDFloat, SIMD_WIDTH> km1_v = set1<SIMDFloat, SIMD_WIDTH>(1.0 - k);
  SIMDVec<SIMDFloat, SIMD_WIDTH> k_v   = set1<SIMDFloat, SIMD_WIDTH>(k);
  for (size_t p = 0; p < orig.size; p += simd_elems) {
    SIMDVec<SIMDFloat, SIMD_WIDTH> det_v, trace_v, trace2_v, A_v, B_v, C_v, R_v;
    A_v      = load<SIMD_WIDTH, SIMDFloat>(A.data + p);
    B_v      = load<SIMD_WIDTH, SIMDFloat>(B.data + p);
    C_v      = load<SIMD_WIDTH, SIMDFloat>(C.data + p);
    det_v    = A_v * B_v - C_v * C_v;
    trace_v  = A_v + B_v;
    trace2_v = trace_v * trace_v;
    if (VARIANT)
      // variant: (1-k)*det-k*trace^2
      R_v = km1_v * det_v - k_v * trace2_v;
    else
      // original: det-k*trace^2
      R_v = det_v - k_v * trace2_v;
    if (NORMALIZE == 1)
      // sgn(R) * sqrt(abs(R))
      R_v = sign(sqrt(abs(R_v)), R_v);
    else if (NORMALIZE == 2)
      // sgn(R) * sqrt(sqrt(abs(R)))
      R_v = sign(sqrt(sqrt(abs(R_v))), R_v);
    store(corner.data + p, R_v);
  }
}

// ---------------------------------------------------------------------------
// harrisFilterCyclic with compile-time output selection
// contributed by Benedikt Volkmer (modified: SIMDFloat only)
// ---------------------------------------------------------------------------

typedef uint8_t HarrisOutputSelector;
#define HOS_VEDGE    1
#define HOS_HEDGE    2
#define HOS_DET      4
#define HOS_TR       8
#define HOS_RESPONSE 16
#define HOS_L1       32
#define HOS_L2       64
/*constexpr*/ bool HOShas(HarrisOutputSelector selector, uint8_t output)
{
  return selector & output;
}
/*constexpr*/ int HOScount(HarrisOutputSelector selector)
{
  const uint8_t NIBBLE_LOOKUP[16] = {0, 1, 1, 2, 1, 2, 2, 3,
                                     1, 2, 2, 3, 2, 3, 3, 4};
  return NIBBLE_LOOKUP[selector & 0x0F] + NIBBLE_LOOKUP[selector >> 4];
}
/*constexpr*/ int HOSidx(HarrisOutputSelector selector, uint8_t output)
{
  int idx = 0;
  if (HOShas(selector, HOS_VEDGE)) {
    if (HOS_VEDGE == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_HEDGE)) {
    if (HOS_HEDGE == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_RESPONSE)) {
    if (HOS_RESPONSE == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_DET)) {
    if (HOS_DET == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_TR)) {
    if (HOS_TR == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_L1)) {
    if (HOS_L1 == output)
      return idx;
    else
      idx++;
  }
  if (HOShas(selector, HOS_L2)) {
    if (HOS_L2 == output)
      return idx;
    else
      idx++;
  }
  throw SIMDException("HOSidx",
                      "Index does not exists for not selected ouputs");
}

template <int NORMALIZE, bool VARIANT, HarrisOutputSelector HOS, int SIMD_WIDTH,
          int SIMD_ALIGN, class InAddOn, class OutAddOn>
void harrisFilterCyclic(
  const SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, InAddOn> &orig,
  SIMDStack<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, OutAddOn> &output,
  const double k, const int binomialFilterCount)
{
  output.resize(HOScount(HOS), orig.w, orig.h, orig.addOn);

  SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN> X, Y, A, B, C;
  verticalEdgeFilter2(orig, X);
  horizontalEdgeFilterCyclic2(orig, Y);
  if /*constexpr*/ (HOShas(HOS, HOS_VEDGE))
    X.copyTo(output[HOSidx(HOS, HOS_VEDGE)].data);
  if /*constexpr*/ (HOShas(HOS, HOS_HEDGE))
    Y.copyTo(output[HOSidx(HOS, HOS_HEDGE)].data);

  mul(X, X, A);
  mul(Y, Y, B);
  mul(X, Y, C);
  // TODO: more cache-friendly with 3 loops?
  for (int i = 0; i < binomialFilterCount; ++i) {
    horizontalBinomialFilterCyclic(A, A);
    verticalBinomialFilter(A, A);
    horizontalBinomialFilterCyclic(B, B);
    verticalBinomialFilter(B, B);
    horizontalBinomialFilterCyclic(C, C);
    verticalBinomialFilter(C, C);
  }
  const int simd_elems = numSIMDVecElements<SIMDFloat, SIMD_WIDTH>();
  // for variant:
  SIMDVec<SIMDFloat, SIMD_WIDTH> km1_v  = set1<SIMDFloat, SIMD_WIDTH>(1.0 - k),
                                 k_v    = set1<SIMDFloat, SIMD_WIDTH>(k),
                                 two_v  = set1<SIMDFloat, SIMD_WIDTH>(2.0),
                                 zero_v = setzero<SIMDFloat, SIMD_WIDTH>();
  for (size_t p = 0; p < orig.size; p += simd_elems) {
    SIMDVec<SIMDFloat, SIMD_WIDTH> det_v, trace_v, trace2_v, A_v, B_v, C_v, R_v;
    A_v     = load<SIMD_WIDTH, SIMDFloat>(A.data + p);
    B_v     = load<SIMD_WIDTH, SIMDFloat>(B.data + p);
    C_v     = load<SIMD_WIDTH, SIMDFloat>(C.data + p);
    det_v   = A_v * B_v - C_v * C_v;
    trace_v = A_v + B_v;
    if /*constexpr*/ (HOShas(HOS, HOS_DET)) {
      if (NORMALIZE > 0)
        store(output[HOSidx(HOS, HOS_DET)].data + p, sqrt(det_v));
      else
        store(output[HOSidx(HOS, HOS_DET)].data + p, det_v);
    }
    if /*constexpr*/ (HOShas(HOS, HOS_TR)) {
      if (NORMALIZE > 0)
        store(output[HOSidx(HOS, HOS_TR)].data + p, sqrt(trace_v));
      else
        store(output[HOSidx(HOS, HOS_TR)].data + p, trace_v);
    }
    if /*constexpr*/ (HOShas(HOS, HOS_RESPONSE)) {
      trace2_v = trace_v * trace_v;
      if (VARIANT)
        // variant: (1-k)*det-k*trace^2
        R_v = km1_v * det_v - k_v * trace2_v;
      else
        // original: det-k*trace^2
        R_v = det_v - k_v * trace2_v;
      if (NORMALIZE == 1)
        // sgn(R) * sqrt(abs(R))
        R_v = sign(sqrt(abs(R_v)), R_v);
      else if (NORMALIZE == 2)
        // sgn(R) * sqrt(sqrt(abs(R)))
        R_v = sign(sqrt(sqrt(abs(R_v))), R_v);
      store(output[HOSidx(HOS, HOS_RESPONSE)].data + p, R_v);
    }
    if /*constexpr*/ (HOShas(HOS, HOS_L1) || HOShas(HOS, HOS_L2)) {
      SIMDVec<SIMDFloat, SIMD_WIDTH> sqrt_v;
      trace_v = div((trace_v), two_v); // trM/2
      //       sqrt(        (trM/2)^2       -detM  )
      sqrt_v = sqrt(sub(mul(trace_v, trace_v), det_v));
      if /*constexpr*/ (HOShas(HOS, HOS_L1)) {
        R_v = trace_v - sqrt_v; // l1 = trM/2-sqrt(trM^2/4-detM)
        if (NORMALIZE > 0)
          store(output[HOSidx(HOS, HOS_L1)].data + p, sqrt(max(R_v, zero_v)));
        else
          store(output[HOSidx(HOS, HOS_L1)].data + p, R_v);
      }
      if /*constexpr*/ (HOShas(HOS, HOS_L2)) {
        R_v = trace_v + sqrt_v; // l2 = trM/2+sqrt(trM^2/4-detM)
        if (NORMALIZE > 0)
          store(output[HOSidx(HOS, HOS_L2)].data + p, sqrt(max(R_v, zero_v)));
        else
          store(output[HOSidx(HOS, HOS_L2)].data + p, R_v);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// countNonZero
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Count non zero pixels in SIMDImage
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
SIMDInt countNonZero(
  const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &image)
{
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Construct vector filled with zeros
  SIMDVec<Type, SIMD_WIDTH> zero = setzero<Type, SIMD_WIDTH>();
  // Construct vector filled with ones
  SIMDVec<Type, SIMD_WIDTH> one = set1<Type, SIMD_WIDTH>(1);
  // Zero initialize a counter vector (count ints)
  SIMDVec<SIMDInt, SIMD_WIDTH> counter = setzero<SIMDInt, SIMD_WIDTH>();
  // Uninitialized vector later used to load comparisons
  SIMDVec<Type, SIMD_WIDTH> compare;
  // Extend compare vector to SIMDInt
  SIMDVec<SIMDInt, SIMD_WIDTH>
    compare_extended[numOutputSIMDVecs<SIMDInt, Type>()];
  // Iterate image rows (row-major order)
  for (int r = 0; r < image.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < image.w; c += simd_elems) {
      // Load next simd vector and compare to zero vector. All non
      // zeros mapped to type max value
      compare = cmpeq(zero, load<SIMD_WIDTH>(&image[r][c]));
      // Map max values to 1
      compare = ifelse(compare, zero, one);
      // Extend to SIMDInt
      convert(&compare, compare_extended);
      // Iterate over extended compare vector
      for (unsigned int i = 0; i < numOutputSIMDVecs<SIMDInt, Type>(); ++i)
        // Accumulate in counter
        counter = adds(counter, compare_extended[i]);
    }
  }
  // Add all elements in the counter to get the number of non zeros in
  // the whole image
  return hadds(counter);
}

// ---------------------------------------------------------------------------
// croppedView
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Creates a view to cropped (at top/bottom) portion of a SIMDImage.
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn = DefaultAddOn>
void croppedView(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
                 const int upper, const int lower,
                 SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Upper part to crop
  int cropOffset = upper;
  // Height top crop of image
  int cropH = upper + lower;
  // Crop by making destination image a view into the src image portion
  dst.view(src, src.w * cropOffset, src.w, src.h - cropH, src.addOn);
}

/**
 * @brief Creates a view to cropped (at top/bottom) portion of a SIMDImage with
 * Panorama addon.
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
void croppedView(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &src,
                 const int upper, const int lower,
                 SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, Panorama> &dst)
{
  // Upper part to crop
  int cropOffset = upper;
  // Height top crop of image
  int cropH = upper + lower;
  // Panorama: Horizon after cropping
  double cropHorizon = src.addOn.horizon - cropOffset;
  // Crop by making destination image a view into the src image portion
  dst.view(src, src.w * cropOffset, src.w, src.h - cropH,
           Panorama(src.addOn.verticalResolution, cropHorizon));
}

// ---------------------------------------------------------------------------
// equalizeHist
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Applies histogram equalization to a masked region of the input image.
 *
 * @note Similar to OpenCV cv::equalizeHist
 *
 * @param src Source SIMDImage of arbitrary type
 * @param dst Destination SIMDImage of same type as src
 * @param mask Mask SIMDImage of same type as src
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
void equalizeHist(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
                  SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst,
                  const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &mask)
{
  // Histogram equalization is only defined for SIMDByte images
  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> src_byte, mask_byte;
  // Copy (convert) source and mask image to SIMDByte images
  copy(src, src_byte);
  copy(mask, mask_byte);
  // Create empty destination image as SIMDByte image
  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> dst_byte(src.w, src.h,
                                                              src.addOn);
  // Use SIMDByte image implementation
  equalizeHist(src_byte, dst_byte, mask_byte);
  // Copy (convert) result back to destination (Type)
  copy(dst_byte, dst);
}

/**
 * @brief Applies histogram equalization to a masked region of the input image.
 *
 * @note Similar to OpenCV cv::equalizeHist
 *
 * @param src Source SIMDImage of arbitrary type
 * @param dst Destination SIMDImage of same type as src
 * @param mask Mask SIMDImage of same type as src
 */
template <int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
void equalizeHist(
  const SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst,
  // @formatter:off
  const SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> &mask)
{
  // @formatter:on
  // Mask and source image must match in size
  assert(mask.sameSize(src));
  // Resize the destination to fit source image size
  dst.resize(src.w, src.h, src.addOn);
  // Histogram and lookup table vectors
  std::vector<SIMDInt> histogram(SIMDTypeInfo<SIMDByte>::max() + 1, 0);
  std::vector<SIMDByte> lookup_table(SIMDTypeInfo<SIMDByte>::max() + 1, 0);
  // Fill histogram of pixel values
  for (int r = 0; r < src.h; ++r) {
    for (int c = 0; c < src.w; ++c) {
      if (mask[r][c]) { histogram[src[r][c]] += 1; }
    }
  }
  // Accumulate sum and scale to normalize histogram to sum 255
  SIMDFloat sum = 0.0, scale = 255.0 / countNonZero(mask);
  // Accumulate histogram and normalize to fill lookup-table
  for (std::size_t index = 0; index < histogram.size(); ++index) {
    sum += (SIMDFloat) histogram[index];
    lookup_table[index] = saturate_cast<SIMDByte>(scale * sum);
  }
  // Transform from source to destination image via lookup table
  for (int r = 0; r < src.h; ++r) {
    for (int c = 0; c < src.w; ++c) {
      // Apply transform only to masked region
      if (mask[r][c]) {
        dst[r][c] = lookup_table[src[r][c]];
      } else {
        dst[r][c] = src[r][c];
      }
    }
  }
}

/**
 * @brief Applies histogram equalization to the input image.
 *
 * @note Similar to OpenCV cv::equalizeHist
 *
 * @param src Source SIMDImage of arbitrary type
 * @param dst Destination SIMDImage of same type as src
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
void equalizeHist(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
                  SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Histogram equalization is only defined for SIMDByte images
  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> src_byte;
  // Copy (convert) source image to SIMDByte images
  copy(src, src_byte);
  // Create empty destination image as SIMDByte image
  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> dst_byte(src.w, src.h,
                                                              src.addOn);
  // Use SIMDByte image implementation
  equalizeHist(src_byte, dst_byte);
  // Copy (convert) result back to destination (Type)
  copy(dst_byte, dst);
}

/**
 * @brief Applies histogram equalization to the input image.
 *
 * @note Similar to OpenCV cv::equalizeHist
 *
 * @param src Source SIMDImage of arbitrary type
 * @param dst Destination SIMDImage of same type as src
 */
template <int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
void equalizeHist(const SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
                  SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Resize the destination to fit source image size
  dst.resize(src.w, src.h, src.addOn);
  // Histogram and lookup table vectors
  std::vector<SIMDInt> histogram(SIMDTypeInfo<SIMDByte>::max() + 1, 0);
  std::vector<SIMDByte> lookup_table(SIMDTypeInfo<SIMDByte>::max() + 1, 0);
  // Fill histogram of pixel values
  for (int r = 0; r < src.h; ++r) {
    for (int c = 0; c < src.w; ++c) { histogram[src[r][c]] += 1; }
  }
  // Accumulate sum and scale to normalize histogram to sum 255
  SIMDFloat sum = 0.0, scale = 255.0 / (src.w * src.h);
  // Accumulate histogram and normalize to fill lookup-table
  for (std::size_t index = 0; index < histogram.size(); ++index) {
    sum += (SIMDFloat) histogram[index];
    lookup_table[index] = saturate_cast<SIMDByte>(scale * sum);
  }
  // Transform from source to destination image via lookup table
  for (int r = 0; r < src.h; ++r) {
    for (int c = 0; c < src.w; ++c) { dst[r][c] = lookup_table[src[r][c]]; }
  }
}

// ---------------------------------------------------------------------------
// interpolate
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

// Substitute for apparently slow std::lround used by NEAREST interpolation
//  Rounds double to long via type-cast
static inline long fast_round(const double x)
{
  return static_cast<long>(x + std::copysign(0.5, x));
}

/**
 * @brief Interpolates pixel value of SIMDImage
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @tparam AddOn SIMDImage addon
 * @param src The image to interpolate from
 * @param i Pixel row index (double: might be between pixels)
 * @param j Pixel column index (double: might be between pixels)
 * @param method Interpolation method to use (NEAREST or LINEAR)
 * @param border Border value to return if the pixel is outside the image
 * @return Interpolated pixel value at (i,j) of src value Type.
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn>
static SIMD_INLINE Type
interpolate(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src, double i,
            double j, const Interpolation method = Interpolation::NEAREST,
            const Type border = SIMDTypeInfo<Type>::max())
{
  // Fill out of bounds with border value
  if (i < 0 || i >= src.h || j < 0 || j >= src.w) { return border; }
  // Select interpolation method
  switch (method) {
  // Nearest neighbor interpolation
  case Interpolation::NEAREST: {
    // If the pixel is out of the src after rounding...
    if (fast_round(i) >= src.h || fast_round(j) >= src.w) {
      // Fill with border value
      return border;
    }
    // Copy pixels from source to destination
    return src[fast_round(i)][fast_round(j)];
  }
    // (Bi-)Linear interpolation
  case Interpolation::LINEAR: {
    // Unit square surrounding mapping (x,y)
    int i0 = std::floor(i), i1 = (int) (std::floor(i) + 1), j0 = std::floor(j),
        j1 = (int) (std::floor(j) + 1);
    // If the pixel is out of the src after rounding...
    if (i1 >= src.h || j1 >= src.w) {
      // Fill with border value
      return border;
    }
    // Get pixels from unit square in source src surrounding
    // mapping (x,y)
    // @formatter:off
    Type f00 = src[i0][j0], f10 = src[i1][j0], f01 = src[i0][j1],
         f11 = src[i1][j1];
    // (Bi-)Linear interpolation formula on unit square
    return f00 * (1 - i + i0) * (1 - j + j0) + f10 * (i - i0) * (1 - j + j0) +
           f01 * (1 - i + i0) * (j - j0) + f11 * (i - i0) * (j - j0);
    // @formatter:on
  }
    // Unsupported method
  default: {
    // Adhoc exception class
    class Exception : public std::exception
    {
      // Returns hardcoded message
      const char *what() const noexcept override
      {
        return "Interpolation Method Not Implemented";
      }
    };
    // Throw the adhoc exception
    throw Exception();
  }
  }
}

// Specialization for SIMDImage with Panorama addon (allows horizontal wrap
// around)
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN>
static SIMD_INLINE Type interpolate(
  const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, ns_simd::Panorama> &src,
  double i, double j, const Interpolation method = Interpolation::NEAREST,
  const Type border = SIMDTypeInfo<Type>::max())
{
  // Fill out of bounds with border value
  //  NOTE: Only on bottom/top border - left/right wraps around
  if (i < 0 || i >= src.h) { return border; }
  // Select interpolation method
  switch (method) {
  // Nearest neighbor interpolation
  case Interpolation::NEAREST: {
    // If the pixel is out of the src after rounding...
    if (fast_round(i) >= src.h) {
      // Fill with border value
      return border;
    }
    // Wrap around horizontal index, avoid % operation
    j = (j >= src.w) ? (j - src.w) : ((j < 0) ? (j + src.w) : j);
    // Copy pixels from source to destination
    return src[fast_round(i)][fast_round(j)];
  }
    // (Bi-)Linear interpolation
  case Interpolation::LINEAR: {
    // Unit square surrounding mapping (x,y)
    int i0 = std::floor(i), i1 = (int) (std::floor(i) + 1), j0 = std::floor(j),
        j1 = (int) (std::floor(j) + 1);
    // If the pixel is out of the src after rounding...
    if (i1 >= src.h) {
      // Fill with border value
      return border;
    }
    // Wrap around horizontal index, avoid % operation
    j0 = (j0 >= src.w) ? (j0 - src.w) : ((j0 < 0) ? (j0 + src.w) : j0);
    j1 = (j1 >= src.w) ? (j1 - src.w) : ((j1 < 0) ? (j1 + src.w) : j1);
    // Need to wrap around the original interpolation coordinate to
    // not have exploding differences in interpolation formula.
    j = (j >= src.w) ? (j - src.w) : ((j < 0) ? (j + src.w) : j);
    // Get pixels from unit square in source src surrounding
    // mapping (x,y)
    // @formatter:off
    Type f00 = src[i0][j0], f10 = src[i1][j0], f01 = src[i0][j1],
         f11 = src[i1][j1];
    // (Bi-)Linear interpolation formula on unit square
    return f00 * (1 - i + i0) * (1 - j + j0) + f10 * (i - i0) * (1 - j + j0) +
           f01 * (1 - i + i0) * (j - j0) + f11 * (i - i0) * (j - j0);
    // @formatter:on
  }
    // Unsupported method
  default: {
    // Adhoc exception class
    class Exception : public std::exception
    {
      // Returns hardcoded message
      const char *what() const noexcept override
      {
        return "Interpolation Method Not Implemented";
      }
    };
    // Throw the adhoc exception
    throw Exception();
  }
  }
}

// ---------------------------------------------------------------------------
// maskImage
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Blends in as mask image into a SIMDImage
 *
 * @param src SIMDImage to be masked
 * @param mask Mask with zero pixels indicating invalid (to be masked) pixels
 * @param dst Destination image (masked image is written here)
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn = DefaultAddOn>
void maskImage(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
               const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &mask,
               SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Resize destination to fit the source image
  dst.resize(src.w, src.h, src.addOn);
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Zero initialized vector to compare the mask pixels to (zero pixel is
  // masked)
  SIMDVec<Type, SIMD_WIDTH> zero = setzero<Type, SIMD_WIDTH>();
  // Fill vector with max value of the image Type which is used to
  // indicate invalid pixels.
  SIMDVec<Type, SIMD_WIDTH> max =
    set1<Type, SIMD_WIDTH>(SIMDTypeInfo<Type>::max());
  // Uninitialized vector later used to load comparisons
  SIMDVec<Type, SIMD_WIDTH> compare;
  // Iterate image rows (row-major order)
  for (int r = 0; r < src.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < src.w; c += simd_elems) {
      // Load next pixel block from mask and compare to zero vector.
      // All non zeros mapped to type max value
      compare = cmpeq(zero, load<SIMD_WIDTH>(&mask[r][c]));
      // Blend in max value for masked pixels
      store(&dst[r][c], ifelse(compare, max, load<SIMD_WIDTH>(&src[r][c])));
    }
  }
}

// ---------------------------------------------------------------------------
// pixelSum
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Sums all pixel values over an image
 */
template <class SumType = SIMDInt, class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn   = DefaultAddOn>
SumType pixelSum(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &image)
{
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Uninitialized vector for loading pixel blocks into
  SIMDVec<Type, SIMD_WIDTH> inVec;
  // Vector to hold input pixel block extended to SumType
  SIMDVec<SumType, SIMD_WIDTH>
    inVecExtended[numOutputSIMDVecs<SumType, Type>()];
  // Zero initialized vector for accumulating parts of sum
  SIMDVec<SumType, SIMD_WIDTH> sumVec = setzero<SumType, SIMD_WIDTH>();
  // Iterate image rows (row-major order)
  for (int r = 0; r < image.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < image.w; c += simd_elems) {
      // Load next pixel block
      inVec = load<SIMD_WIDTH>(&image[r][c]);
      // Extend to SumType
      convert(&inVec, inVecExtended);
      // Iterate over extended input vector
      for (unsigned int i = 0; i < numOutputSIMDVecs<SumType, Type>(); ++i)
        // Accumulate in counter
        sumVec = adds(sumVec, inVecExtended[i]);
    }
  }
  // Add all elements in the sum to get the pixel sum over the whole
  // image
  return hadds(sumVec);
}

/**
 * @brief Sums all pixel values over an image considering a mask image
 */
template <class SumType = SIMDInt, class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn   = DefaultAddOn>
SumType pixelSum(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &image,
                 const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &mask)
{
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Uninitialized vector for loading pixel blocks into
  SIMDVec<Type, SIMD_WIDTH> inVec;
  // Vector to hold input pixel block extended to SumType
  SIMDVec<SumType, SIMD_WIDTH>
    inVecExtended[numOutputSIMDVecs<SumType, Type>()];
  // Zero initialized vector for accumulating parts of sum
  SIMDVec<SumType, SIMD_WIDTH> sumVec = setzero<SumType, SIMD_WIDTH>();
  // Zero initialized vector to compare the mask pixels to (zero pixel is
  // masked)
  SIMDVec<Type, SIMD_WIDTH> zero = setzero<Type, SIMD_WIDTH>();
  // Uninitialized vector later used to load comparisons
  SIMDVec<Type, SIMD_WIDTH> compare;
  // Iterate image rows (row-major order)
  for (int r = 0; r < image.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < image.w; c += simd_elems) {
      // Load next pixel block from mask
      compare = cmpeq(zero, load<SIMD_WIDTH>(&mask[r][c]));
      // Load next pixel block from image and set all masked pixels to
      // zero
      inVec = ifelse(compare, zero, load<SIMD_WIDTH>(&image[r][c]));
      // Extend to SumType
      convert(&inVec, inVecExtended);
      // Iterate over extended input vector
      for (unsigned int i = 0; i < numOutputSIMDVecs<SumType, Type>(); ++i)
        // Accumulate in counter
        sumVec = adds(sumVec, inVecExtended[i]);
    }
  }

  // Add all elements in the sum to get the pixel sum over the whole
  // image
  return hadds(sumVec);
}

// ---------------------------------------------------------------------------
// normalizedPixelSum
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Sums all pixel values over an image and normalizes the result by
 * dividing by the number of pixels in the image.
 */
template <typename SumType = SIMDFloat, class Image>
SumType normalizedPixelSum(const Image &image)
{
  // Normalize by number of pixels in the image
  return pixelSum<SumType>(image) / (SumType) (image.w * image.h);
}

/**
 * @brief Sums all pixel values over an image and normalizes the result by
 * dividing by the number of not masked pixels in the image.
 */
template <typename SumType = SIMDFloat, class Image>
SumType normalizedPixelSum(const Image &image, const Image &mask)
{
  // Normalize by number of pixels in the image
  return pixelSum<SumType>(image, mask) / (SumType) countNonZero(mask);
}

// ---------------------------------------------------------------------------
// substituteValue
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

/**
 * @brief Substitutes all occurrences of a value within a SIMDImage
 *
 * @param src Source image to processed
 * @param valueIn Value to be substituted
 * @param valueOut Value used to substitute
 * @param dst Destination image (processed image is written here)
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn = DefaultAddOn>
void substituteValue(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src,
                     const Type valueIn, const Type valueOut,
                     SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Resize destination to fit the source image
  dst.resize(src.w, src.h, src.addOn);
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Fill vector with the input value (value to be replaced)
  SIMDVec<Type, SIMD_WIDTH> vecValueIn = set1<Type, SIMD_WIDTH>(valueIn);
  // Fill vector with the output value (value to replace with)
  SIMDVec<Type, SIMD_WIDTH> vecValueOut = set1<Type, SIMD_WIDTH>(valueOut);
  // Uninitialized vector to load pixel blocks into
  SIMDVec<Type, SIMD_WIDTH> vecIn;
  // Uninitialized vector later used to load comparisons
  SIMDVec<Type, SIMD_WIDTH> compare;
  // Iterate image rows (row-major order)
  for (int r = 0; r < src.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < src.w; c += simd_elems) {
      // Load next pixel block into SIMD vector
      vecIn = load<SIMD_WIDTH>(&src[r][c]);
      // Compare to input value vector
      compare = cmpeq(vecValueIn, vecIn);
      // If equal to valueIn, replace with valueOut
      store(&dst[r][c], ifelse(compare, vecValueOut, vecIn));
    }
  }
}

// ---------------------------------------------------------------------------
// SIMDImage remapping
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------

// Image remapping utility
namespace Remap {
/**
 * @brief Image remapping map
 */
struct Map
{
  /**
   * @brief Width of the image mapping (width of the target image)
   */
  unsigned int w;

  /**
   * @brief Height of the image mapping (height of the target image)
   */
  unsigned int h;

  /**
   * @brief Initialized (empty) map of requested width x height
   */
  Map(const unsigned int w, const unsigned int h)
    : w {w}, h {h}, map {h, std::vector<double>(w)}
  {}

  /**
   * @brief Access mapping value at coordinates (read+write)
   */
  inline double &operator()(unsigned int x, unsigned int y)
  {
    // Delegate to std::vector []; No bound checking!
    return map[x][y];
  }

  /**
   * @brief Access mapping value at coordinates (read only)
   */
  inline double operator()(unsigned int x, unsigned int y) const
  {
    // Delegate to std::vector []; No bound checking!
    return map[x][y];
  }

private:
  // Mapping represented as 2d vector
  std::vector<std::vector<double>> map;
};
} // namespace Remap

/**
 * @brief OpenCV cv::remap like SIMDImage remapping
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @tparam AddOn1 SIMDImage addon of source image
 * @tparam AddOn2 SIMDImage addon of destination image
 * @param src Source image to remap
 * @param dst Destination image to map to
 * @param mapx Row coordinate mapping
 * @param mapy Column coordinate mapping
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn1,
          class AddOn2>
void remap(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn1> &src,
           SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn2> &dst,
           const Remap::Map &mapx, const Remap::Map &mapy,
           const Interpolation interpolation = Interpolation::NEAREST,
           const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Input and output may not be identical nor overlap
  assert(!ns_simd::same(src, dst) && !ns_simd::overlap(src, dst));
  // Map sizes for rows and columns must match
  assert(mapx.w == mapy.w && mapx.h == mapy.h);
  // Check whether output image already has requested size
  if (dst.w != (int) mapx.w || dst.h != (int) mapx.h)
    // Resize output image, keep output image addon
    dst.resize(mapx.w, mapy.h, dst.addOn);
  // Process image pixel-wise in row-major order
  for (int y = 0; y < dst.w; ++y) {
    for (int x = 0; x < dst.h; ++x) {
      // Get mapping coordinates from row/column maps
      double mx = mapx(x, y), my = mapy(x, y);
      // Copy/Interpolate pixel from source to destination
      dst[x][y] = interpolate(src, mx, my, interpolation, border);
    }
  }
}

// ---------------------------------------------------------------------------
// SIMDImage columns-/row-wise shifting (kind of special case of remapping)
// contributed by Christoph Berganski
// ---------------------------------------------------------------------------
// Image shifting utility
namespace Shift {
/**
 * @brief Row shifting map (offset table)
 */
struct ShiftRow
{
  /**
   * @brief Height/Number of rows of image to shift
   */
  unsigned int h;

  /**
   * @brief Initializes (empty) offset table of requested height
   */
  explicit ShiftRow(const unsigned int h) : h {h}, shift(h) {}

  /**
   * @brief Access shift offset value at coordinates (read+write)
   */
  inline double &operator()(unsigned int x)
  {
    // Delegate to std::vector []; No bound checking!
    return shift[x];
  }

  /**
   * @brief Access shift offset value at coordinates (read only)
   */
  inline double operator()(unsigned int x) const
  {
    // Delegate to std::vector []; No bound checking!
    return shift[x];
  }

private:
  // Shift offsets table
  std::vector<double> shift;
};

/**
 * @brief Column shifting map (offset table)
 */
struct ShiftCol
{
  /**
   * @brief Width/Number of cols of image to shift
   */
  unsigned int w;

  /**
   * @brief Initializes (empty) offset table of requested width
   */
  explicit ShiftCol(const unsigned int w) : w {w}, shift(w) {}

  /**
   * @brief Access shift offset value at coordinates (read+write)
   */
  inline double &operator()(unsigned int y)
  {
    // Delegate to std::vector []; No bound checking!
    return shift[y];
  }

  /**
   * @brief Access shift offset value at coordinates (read only)
   */
  inline double operator()(unsigned int y) const
  {
    // Delegate to std::vector []; No bound checking!
    return shift[y];
  }

private:
  // Shift offset table
  std::vector<double> shift;
};
} // namespace Shift

/**
 * @brief Shifts SIMDImage columns
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @tparam AddOn1 SIMDImage addon of source image
 * @tparam AddOn2 SIMDImage addon of destination image
 * @param src Source image to shift
 * @param dst Destination image to map to
 * @param shift Column shift offset table
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn1,
          class AddOn2>
void shift(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn1> &src,
           SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn2> &dst,
           const Shift::ShiftCol &shift,
           const Interpolation interpolation = Interpolation::NEAREST,
           const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Input and output may not be identical nor overlap
  assert(!ns_simd::same(src, dst) && !ns_simd::overlap(src, dst));
  // Offset table size must match with image height
  assert((int) shift.w == src.w);
  // Output image needs to have the same size as input
  if (dst.w != src.w || dst.h != src.h)
    // Resize output image, keep output image addon
    dst.resize(src.w, src.h, dst.addOn);
  // Process image pixel-wise in row-major order
  for (int y = 0; y < dst.w; ++y) {
    for (int x = 0; x < dst.h; ++x) {
      // Get offset from table to shift along row
      double mx = x + shift(y);
      // Copy/Interpolate pixel from source to destination
      dst[x][y] = interpolate(src, mx, y, interpolation, border);
    }
  }
}

/**
 * @brief Shifts SIMDImage rows
 * @tparam Type Value type of the SIMDImage
 * @tparam SIMD_WIDTH Width of SIMD vectors
 * @tparam SIMD_ALIGN Memory alignment of SIMD vector
 * @tparam AddOn1 SIMDImage addon of source image
 * @tparam AddOn2 SIMDImage addon of destination image
 * @param src Source image to shift
 * @param dst Destination image to map to
 * @param shift Row shift offset table
 * @param interpolation Interpolation method to use for remapping non-integer
 * coordinates
 * @param border Border value to insert if a pixel is outside the image
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN, class AddOn1,
          class AddOn2>
void shift(const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn1> &src,
           SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn2> &dst,
           const Shift::ShiftRow &shift,
           const Interpolation interpolation = Interpolation::NEAREST,
           const Type border                 = SIMDTypeInfo<Type>::max())
{
  // Input and output may not be identical nor overlap
  assert(!ns_simd::same(src, dst) && !ns_simd::overlap(src, dst));
  // Offset table size must match with image height
  assert((int) shift.h == src.h);
  // Output image needs to have the same size as input
  if (dst.w != src.w || dst.h != src.h)
    // Resize output image, keep output image addon
    dst.resize(src.w, src.h, dst.addOn);
  // Process image pixel-wise in row-major order
  for (int y = 0; y < dst.w; ++y) {
    for (int x = 0; x < dst.h; ++x) {
      // Get offset from table to shift along row
      double my = y + shift(x);
      // Copy/Interpolate pixel from source to destination
      dst[x][y] = interpolate(src, x, my, interpolation, border);
    }
  }
}

/**
 * @brief Computes pixel-wise absolute difference of SIMDImages
 */
template <class Type, int SIMD_WIDTH, int SIMD_ALIGN,
          class AddOn = DefaultAddOn>
void absoluteDifference(
  const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src1,
  const SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &src2,
  SIMDImage<Type, SIMD_WIDTH, SIMD_ALIGN, AddOn> &dst)
{
  // Source images must match in size and addon
  assert(src1.w == src2.w && src1.h == src2.h && src1.addOn == src2.addOn);
  // Resize destination to fit the source image
  dst.resize(src1.w, src1.h, src1.addOn);
  // Number of elements in SIMD vector of given width and type
  int simd_elems = SIMD_WIDTH / sizeof(Type);
  // Iterate image rows (row-major order)
  for (int r = 0; r < src1.h; ++r) {
    // Iterate images columns in steps of simd_elems
    for (int c = 0; c < src1.w; c += simd_elems) {
      // Compute absolute difference of next pixel block from source
      // images and store to destination image
      store(&dst[r][c], absDiff(load<SIMD_WIDTH>(&src1[r][c]),
                                load<SIMD_WIDTH>(&src2[r][c])));
    }
  }
}
} // namespace ns_simd

#endif
