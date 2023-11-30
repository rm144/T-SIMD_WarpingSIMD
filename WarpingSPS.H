// ===========================================================================
//
// WarpingSPS.H --
// C++ template for scale-plane stack
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
#ifndef WARPING_SPS_H_
#define WARPING_SPS_H_

#include "SIMDImage.H"
#include "SIMDVecAll.H"
#include "WarpingParameter.H"

#include <cstddef>
#include <stdio.h>
#include <string.h>

namespace simd {

// -------------------------------------------------------------------------
// SymmetricScaleFactors
// -------------------------------------------------------------------------

// neighboring scale factors are related by constant factor;
// scale factors are inversion-symmetric

class SymmetricScaleFactors
{
public:
  // vector storing all scale parameters
  std::vector<double> scaleFactors;

  SymmetricScaleFactors(int nScalePlanes, double maxScaleFactor)
  {
    scaleFactors.resize(nScalePlanes);
    // center index (nScalePlanes needs to be odd)
    int nc = (nScalePlanes - 1) / 2;
    // initialize scale factors
    scaleFactors[nc] = 1.0;
    // if we have more than 1 scale plane
    if (nc != 0) {
      // factor between scaleFactors
      double fac = exp(log(maxScaleFactor) / nc);
      // set scaleFactors
      double fn = 1.0, fp = 1.0;
      for (int i = 1; i <= nc; i++) {
        fn /= fac;
        fp *= fac;
        scaleFactors[nc - i] = fn;
        scaleFactors[nc + i] = fp;
      }
    }
  }
};

// -------------------------------------------------------------------------
// SingleScaleFactor (for Min2Warping)
// -------------------------------------------------------------------------

class SingleScaleFactor
{
public:
  // vector storing all scale parameters
  std::vector<double> scaleFactors;

  SingleScaleFactor()
  {
    scaleFactors.resize(1);
    scaleFactors[0] = 1.0;
  }
};

// -------------------------------------------------------------------------
// WarpingSPS
// -------------------------------------------------------------------------

// SIMD_WIDTH is the number of vector elements processed in
// parallel in the SPS by a SIMD vector extension;
// SIMD_ALIGN is the alignment requirement of the SIMD vector
// extension (in bytes)
template <typename SPSType, size_t SIMD_WIDTH, size_t SIMD_ALIGN>
class WarpingSPS
{
public:
  // a reference to the joint parameter set
  const WarpingParameter<SPSType, SIMD_WIDTH> &param;
  // we intentionally provide free access to the scale-plane stack
  // from the outside
  SIMDStack<SPSType, SIMD_WIDTH, SIMD_ALIGN> stack;

  WarpingSPS(const WarpingParameter<SPSType, SIMD_WIDTH> &param,
             int numScaleFactors)
    : param(param), stack(numScaleFactors, param.wSPS, param.w)
  {
    // puts("WarpingSPS::WarpingSPS");
  }

  virtual ~WarpingSPS() {}

  bool saveAll(FILE *f) const { return stack().save(f); }

  // couldn't call this load, leads to conflicts with SIMD template load
  bool loadAll(FILE *f) { return stack().load(f); }

  bool sameSize(const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps) const
  {
    // 2 SPS have the same size if the 2 stacks have the same size
    return stack.sameSize(sps.stack);
  }

  void copyFrom(const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps)
  {
    if (!sameSize(sps))
      throw SIMDException("WarpingSPS::copyFrom", "size mismatch");
    // we just copy data, no resize operation
    stack.copy(sps.stack.data);
  }

  // fill in simd_sps_elems more elements in reshuffled scale-plane stack
  void complement()
  {
    SIMDVec<SPSType, SIMD_WIDTH> copy;
    SPSType *src, *dst;
    int off;
    const int w = param.w, wSPS = param.wSPS;

    for (int p = 0; p < stack.numPlanes; p++) {
      // offset to parts
      off = 0;
      // pointer into scale plane
      SPSType *data = stack[p].data;
      // go through all stepAlpha parts
      for (int s = 0; s < param.stepAlpha; s++) {
        // pick plane start
        src = data + off;
        dst = src + param.nAlpha;
        // go through entire first column in each part
        for (int i = 0; i < w; i++) {
          // printf("p = %d, s = %d, i = %d, src = %p, dst = %p (%d)\n",
          //	 p, s, i, src, dst, param.simd_sps_elems);
          // copy entire first column in each part
          copy = load<SIMD_WIDTH>(src);
          store(dst, copy);
          // for (int k = 0; k < param.simd_sps_elems; k++) dst[k] = src[k];
          src += wSPS;
          dst += wSPS;
        }
        off += param.blkSize;
      }
    }
  }

  // inversion of scale-plane stack (ss-cv => cv-ss)
  void invertFrom(const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &spsFrom)
  {
    if (!sameSize(spsFrom))
      throw SIMDException("WarpingSPS::invertFrom", "size mismatch");
    //  1. Oct 14 (rm)
    int w = param.w, wSPS = param.wSPS, simd_sps_elems = param.simd_sps_elems;
    int blkSize = param.blkSize, nAlpha = param.nAlpha;
    // 10. Oct 14 (rm): we just need shuffleModuloPW
    // int *shuffleModulo = param.shuffleModulo;
    int *shuffleModuloPW = param.shuffleModuloPW;
    int *unshuffle       = param.unshuffle;
    int *modulo2wpw      = param.modulo2wpw;
    // (in-place inversion may not be possible, not sure)
    // original scale-plane stack has coordinates
    // iThetaSS, iDeltaCV
    // inverted scale-plane stack has coordinates
    // iThetaCV, iDeltaSS
    int S = stack.numPlanes;
    SIMDVec<SPSType, SIMD_WIDTH> d;
    // 22. Dec 13 (rm): changed loop order: plane -> row -> column
    for (int s = 0; s < S; s++) {
      SPSType *planeTo   = stack[s].data;
      SPSType *planeFrom = spsFrom.stack[(S - 1) - s].data;
      for (int iDeltaCV = 0; iDeltaCV < w; iDeltaCV++) {
        SPSType *planeToRow   = planeTo + modulo2wpw[-iDeltaCV];
        SPSType *planeFromRow = planeFrom + modulo2wpw[iDeltaCV];
        /*
          int *shuffleModulo_plus_iDeltaCV
          = param.shuffleModulo + iDeltaCV;
          for (int iThetaSS = 0; iThetaSS < param.w; iThetaSS++) {
          // compute new coordinates
          // thetaCV = deltaCV + thetaSS
          // int iThetaCV = modulo[iDeltaCV + iThetaSS];
          // deltaSS = -deltaCV;
          // int iDeltaSS = modulopw[-iDeltaCV];
          // copy data (from SPS copy to SPS which is used in computation)
          // note that we invert the order of the scale planes s -> (S-1)-s
          // scalePlane1[s][shuffle[iThetaCV] + wSPS * iDeltaSS] =
          //  scalePlane0[(S-1)-s][shuffle[iThetaSS] + wSPS * iDeltaCV];
          // 22. Dec 13 (rm): faster via table?
          // scalePlane1[s][shuffleModulo[iDeltaCV + iThetaSS]
          //	       + modulo2wpw[-iDeltaCV]] =
          // scalePlane0[(S-1)-s][shuffle[iThetaSS]
          //		       + modulo2wpw[iDeltaCV]];
          planeToRow[shuffleModulo_plus_iDeltaCV[iThetaSS]] =
          planeFromRow[param.shuffle[iThetaSS]];
          }
        */
        //  1. Oct 14 (rm)
        for (int blkStart = 0; blkStart < wSPS; blkStart += blkSize) {
          for (int iDstS = blkStart; iDstS < blkStart + nAlpha;
               iDstS += simd_sps_elems) {
            // 10. Oct 14 (rm):
            // now shuffleModuloPW (negative indices) instead of w +
            // int iSrcS = shuffleModulo[w + unshuffle[iDstS] - iDeltaCV];
            int iSrcS = shuffleModuloPW[unshuffle[iDstS] - iDeltaCV];
            // SIMDVec<SPSType,SIMD_WIDTH>
            d = loadu<SIMD_WIDTH>(planeFromRow + iSrcS);
            store(planeToRow + iDstS, d);
          }
        }
      }
    }
    // complement
    complement();
  }

  // get a compass estimate from the (shuffled) SPS
  // TODO: may be overgeneralized: since we add, we can probably assume that
  // TODO: sizeof(CompassType) >= sizeof(SPSType)
  // TODO: (decided to leave this in as an example)
  template <typename CompassType>
  void compassEstimate(
    SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> &compass) const
  {
    const int wSPS           = param.wSPS;
    const int blkSize        = param.blkSize;
    const int simd_sps_elems = param.simd_sps_elems;
    const int simd_compass_elems =
      numSIMDVecElements<CompassType, SIMD_WIDTH>();
    // SIMD_WIDTH / sizeof(CompassType);
    const int nAlpha     = param.nAlpha;
    const int nPsi       = param.nPsi;
    const int numPlanes  = stack.numPlanes;
    const int numInVecs  = numInputSIMDVecs<CompassType, SPSType>();
    const int numOutVecs = numOutputSIMDVecs<CompassType, SPSType>();
    const int numElems =
      numSIMDVecsElements<CompassType, SPSType, SIMD_WIDTH>();
    if (blkSize % numElems != 0)
      throw SIMDException("WarpingSPS:CompassEstimate",
                          "blkSize has to be a multiple of numElems");
    if (nPsi % simd_compass_elems != 0)
      throw SIMDException("WarpingSPS:CompassEstimate",
                          "nPsi has to be a multiple of simd_compass_elems");
    SIMDVec<CompassType, SIMD_WIDTH> sum[simd_compass_elems];
    SIMDVec<SPSType, SIMD_WIDTH> minv[numInVecs];
    SIMDVec<CompassType, SIMD_WIDTH> minvC[numOutVecs];
    compass.resize(nPsi, 1);
    // we compute a compass estimate for nPsi steps in blocks of
    // simd_compass_elems elements
    for (int iPsi0 = 0; iPsi0 < nPsi; iPsi0 += simd_compass_elems) {
      // we compute simd_compass_elems sums and then add them horizontally
      // SIMDVec<CompassType,SIMD_WIDTH> sum[simd_compass_elems];
      for (int ic = 0, iPsi = iPsi0; ic < simd_compass_elems; ic++, iPsi++) {
        sum[ic] = setzero<CompassType, SIMD_WIDTH>();
        // row offset: the psi-index is subtracted in the warp
        // equations, that's why we have to invert it here
        int rowOff = param.modulo2wpw[-param.jPsiVec[iPsi]];
        // we sum over each SPS row, but leave out the fill-in elements:
        // we have stepAlpha blocks with blkSize elements each in the SPS!
        for (int block = 0; block < wSPS; block += blkSize)
          for (int iAlpha = 0; iAlpha < nAlpha; iAlpha += numElems) {
            // compute minimum through all scale planes
            int off = rowOff + block + iAlpha;
            for (int vi = 0, ioff = 0; vi < numInVecs;
                 vi++, ioff += simd_sps_elems) {
              // printf("block %d, iAlpha = %d, off %u\n", block, iAlpha, off);
              minv[vi] = load<SIMD_WIDTH>(stack[0].data + off + ioff);
              for (int p = 1; p < numPlanes; p++)
                minv[vi] =
                  min(minv[vi], load<SIMD_WIDTH>(stack[p].data + off + ioff));
            }
            // now we convert to CompassType
            convert(minv, minvC);
            // add to sum[ic]
            for (int vo = 0; vo < numOutVecs; vo++)
              sum[ic] = adds(sum[ic], minvC[vo]);
          }
      }
      // add (saturated) all sums horizontally and store into compass
      store(compass.data + iPsi0, hadds(sum));
    }
  }

  // get a compass estimate from the (shuffled) SPS (quantil)
  void quantilCompassEstimate(
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &compass,
    double quantilFraction) const
  {
    const int wSPS           = param.wSPS;
    const int blkSize        = param.blkSize;
    const int nAlpha         = param.nAlpha;
    const int nPsi           = param.nPsi;
    const int quantilElement = int(quantilFraction * (param.w - 1));
    const int numPlanes      = stack.numPlanes;
    std::vector<SPSType> rowVec;
    SIMDVec<SPSType, SIMD_WIDTH> minv;
    // 18. May 16 (rm): was min.elements, clang++ warned about it
    SPSType minValues[SIMDVec<SPSType, SIMD_WIDTH>::elements];
    compass.resize(nPsi, 1);
    // go through all psi values
    for (int iPsi = 0; iPsi < nPsi; iPsi++) {
      rowVec.clear();
      // row offset: the psi-index is subtracted in the warp
      // equations, that's why we have to invert it here
      int rowOff = param.modulo2wpw[-param.jPsiVec[iPsi]];
      // we leave out the fill-in elements:
      // we have stepAlpha blocks with blkSize elements each in the SPS!
      for (int block = 0; block < wSPS; block += blkSize)
        for (int iAlpha = 0; iAlpha < nAlpha; iAlpha += minv.elements) {
          // compute minimum through all scale planes
          int off = rowOff + block + iAlpha;
          minv    = load<SIMD_WIDTH>(stack[0].data + off);
          for (int p = 1; p < numPlanes; p++)
            minv = min(minv, load<SIMD_WIDTH>(stack[p].data + off));
          store(minValues, minv);
          rowVec.insert(rowVec.end(), minValues, minValues + minv.elements);
        }
      // determine and store quantil element
      nth_element(rowVec.begin(), rowVec.begin() + quantilElement,
                  rowVec.end());
      compass.data[iPsi] = rowVec[quantilElement];
    }
  }

  bool validPlane(int i) const { return (i >= 0) && (i < stack.numPlanes); }

  // get plane as it is (shuffled, theta-delta order)
  template <typename ImgType>
  void getPlaneRaw(int i, SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &dest,
                   double shift = 0.0, double scale = 1.0) const
  {
    if (!validPlane(i))
      throw SIMDException("WarpingSPS<SPSType>::getPlane",
                          "scale plane index out of range");
    copyShiftAndScale(stack[i], shift, scale, dest);
  }

  // get unshuffled plane
  template <typename ImgType>
  void getPlaneUnshuffled(int i,
                          SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &dest,
                          double shift = 0.0, double scale = 1.0) const
  {
    if (!validPlane(i))
      throw SIMDException("WarpingSPS::getPlaneUnshuffled",
                          "scale plane index out of range");
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> unshuffled;
    unshufflePlane(stack[i], unshuffled);
    copyShiftAndScale(unshuffled, shift, scale, dest);
  }

  // get unshuffled plane in theta-theta' order
  template <typename ImgType>
  void getPlaneUnshuffledTT(int i,
                            SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &dest,
                            double shift = 0.0, double scale = 1.0) const
  {
    if (!validPlane(i))
      throw SIMDException("WarpingSPS::getPlaneUnshuffledTT",
                          "scale plane index out of range");
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> unshuffledTT;
    unshufflePlaneTT(stack[i], unshuffledTT);
    copyShiftAndScale(unshuffledTT, shift, scale, dest);
  }

  template <typename ImgType>
  void getPlane(int i, int format,
                SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &dest,
                double shift = 0.0, double scale = 1.0) const
  {
    switch (format) {
    case 0: getPlaneRaw(i, dest, shift, scale); break;
    case 1: getPlaneUnshuffled(i, dest, shift, scale); break;
    case 2: getPlaneUnshuffledTT(i, dest, shift, scale); break;
    default: throw SIMDException("WarpingSPS::getPlane", "invalid format");
    }
  }

  // compute minimum over all scale planes and return as Image
  template <typename ImgType>
  void getPlaneMinimumRaw(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &min,
                          double shift = 0.0, double scale = 1.0) const
  {
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> planeMinArray;
    planeMinimum(planeMinArray);
    copyShiftAndScale(planeMinArray, shift, scale, min);
  }

  // compute minimum over all scale planes and return as Image (unshuffled)
  template <typename ImgType>
  void getPlaneMinimumUnshuffled(
    SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &min, double shift = 0.0,
    double scale = 1.0) const
  {
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> planeMinArray;
    planeMinimum(planeMinArray);
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> unshuffled;
    unshufflePlane(planeMinArray, unshuffled);
    copyShiftAndScale(unshuffled, shift, scale, min);
  }

  // compute minimum over all scale planes and return as Image (unshuffled,
  // theta-theta' order)
  template <typename ImgType>
  void getPlaneMinimumUnshuffledTT(
    SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &min, double shift = 0.0,
    double scale = 1.0) const
  {
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> planeMinArray;
    planeMinimum(planeMinArray);
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> unshuffledTT;
    unshufflePlaneTT(planeMinArray, unshuffledTT);
    copyShiftAndScale(unshuffledTT, shift, scale, min);
  }

  template <typename ImgType>
  void getPlaneMinimum(int format,
                       SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &min,
                       double shift = 0.0, double scale = 1.0) const
  {
    switch (format) {
    case 0: getPlaneMinimumRaw(min, shift, scale); break;
    case 1: getPlaneMinimumUnshuffled(min, shift, scale); break;
    case 2: getPlaneMinimumUnshuffledTT(min, shift, scale); break;
    default:
      throw SIMDException("WarpingSPS::getPlaneMinimum", "invalid format");
    }
  }

  template <typename CompassType, typename ImgType>
  void getCompassImage(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &res,
                       double shift = 0.0, double scale = 1.0) const
  {
    SIMDImage<CompassType, SIMD_WIDTH, SIMD_ALIGN> compass;
    compassEstimate<CompassType>(compass);
    copyShiftAndScale(compass, shift, scale, res);
  }

  template <typename ImgType>
  void getQuantilCompassImage(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &res,
                              double quantilFraction, double shift = 0.0,
                              double scale = 1.0) const
  {
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> compass;
    quantilCompassEstimate(compass, quantilFraction);
    copyShiftAndScale(compass, shift, scale, res);
  }

  void planeMinimum(WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &spsMin) const
  {
    if ((param.w != spsMin.param.w) || (param.wSPS != spsMin.param.wSPS) ||
        (spsMin.stack.numPlanes != 1))
      throw SIMDException("WarpingSPS::planeMinimum", "invalid size of spsMin");
    planeMinimum(spsMin.stack[0]);
  }

  //------------------------------------------------------------------------

protected:
  void planeMinimum(SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &minimum) const
  {
    SIMDVec<SPSType, SIMD_WIDTH> minv;
    SPSType *minp;
    int simd_sps_elems = param.simd_sps_elems;

    minimum.resize(param.wSPS, param.w);
    minp = minimum.data;
    for (int e = 0; e < param.wSPS * param.w; e += simd_sps_elems) {
      minv = load<SIMD_WIDTH>(stack[0].data + e);
      for (int p = 1; p < stack.numPlanes; p++)
        minv = min(minv, load<SIMD_WIDTH>(stack[p].data + e));
      store(minp, minv);
      minp += simd_sps_elems;
    }
  }

  // unshuffle single plane
  void unshufflePlane(
    const SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &plane,
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &planeUnshuffled) const
  {
    planeUnshuffled.resize(param.w, param.w);
    for (int iTheta = 0; iTheta < param.w; iTheta++)
      for (int iDelta = 0; iDelta < param.w; iDelta++)
        planeUnshuffled[iDelta][iTheta] = plane[iDelta][param.shuffle[iTheta]];
  }

  // unshuffle single plane and put into theta-theta' order
  void unshufflePlaneTT(
    const SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &plane,
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &planeUnshuffledTT) const
  {
    planeUnshuffledTT.resize(param.w, param.w);
    for (int iTheta = 0; iTheta < param.w; iTheta++)
      for (int iDelta = 0; iDelta < param.w; iDelta++) {
        int iThetaPrime = param.modulo[iTheta + iDelta];
        planeUnshuffledTT[iThetaPrime][iTheta] =
          plane[iDelta][param.shuffle[iTheta]];
      }
  }
};

} // namespace simd

#endif
