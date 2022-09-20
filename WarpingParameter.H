// ===========================================================================
// 
// WarpingParameter.H --
// warping parameters and some lookup tables bundled in a template class
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

#ifndef _WARPING_PARAMETER_H_
#define _WARPING_PARAMETER_H_

#include <cmath>
#include <string.h>
#include "SIMDVecAll.H"
#include "SIMDException.H"

// 20. Sep 22 (Jonas Keller): on Windows M_PI is not defined, so we
// define it here
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ns_simd {

  // we need SPSType (type of scale-plane stack elements) and the vector
  // width in bytes (SIMD_WIDTH) since these influence the dimensions of
  // the SPS which are stored in WarpingParameter
  template <typename SPSType, int SIMD_WIDTH>
  class WarpingParameter
  {
  public:
    // the number of SPS (scale-plane stack) elements processed
    // in parallel in a SIMD vector (determines wSPS and blkSize)
    // TODO: WarpingParameter: check if use of these members is still correct
    const int simd_sps_elems;
    // width of snapshot and current view, also scale planes are width x width
    const int w;
    // angular step size
    const double dTheta;
    // number of steps in search direction alpha (movement direction angle)
    const int nAlpha;
    // number of steps in search direction psi (compass angle)
    const int nPsi;
    // step width for alpha
    const int stepAlpha;
    // step width for psi
    const int stepPsi;
    // step angle for alpha
    const double dAlpha;
    // step angle for psi
    const double dPsi;
    // width of shuffled scale-plane stack
    const int wSPS;
    // size of a single contiguous block in shuffled SPS
    const int blkSize;
  
    // vector storing all discrete alphas
    double *alphaVec;
    // vector storing all discrete psis
    double *psiVec;
    // vector storing pixel indices of all alphas
    int *jAlphaVec;
    // vector storing pixel indices of all psis
    int *jPsiVec;

    // modulo tables for faster computation
    // (we intentionally provide free access to these tables) 
    // modulo2w: modulo w and multiplied with wSPS
    int *modulo2w, *modulo2wpw;
    // modulo: just modulo w
    int *modulo, *modulopw;
    // modulo nAlpha
    int *modNAlpha;
    // shuffle table
    int *shuffle;
    // shuffle table combined with modulo table
    int *shuffleModulo, *shuffleModuloPW;
    // inverse shuffle table
    int *unshuffle;

    WarpingParameter(int width, int nAlpha, int nPsi) :
      // simd_sps_elems(SIMD_WIDTH / sizeof(SPSType)),
      simd_sps_elems(NUM_SIMDVEC_ELEMENTS(SPSType,SIMD_WIDTH)),
      w(width), dTheta((2.0 * M_PI) / w), nAlpha(nAlpha), nPsi(nPsi),
      stepAlpha(w / nAlpha), stepPsi(w / nPsi),
      dAlpha((2.0 * M_PI) / nAlpha), dPsi((2.0 * M_PI) / nPsi),
      // wSPS: width of shuffled scale-plane stack: we need simd_sps_elems
      // more so that we can read simd_sps_elems values in parallel
      // whereever we are in the sps rows
      wSPS(w + stepAlpha * simd_sps_elems), blkSize(nAlpha + simd_sps_elems),
      alphaVec(0), psiVec(0), jAlphaVec(0), jPsiVec(0),
      modulo2w(0), modulo2wpw(0), modulo(0), modulopw(0), modNAlpha(0),
      shuffle(0), shuffleModulo(0), shuffleModuloPW(0), unshuffle(0)
    {
      // TODO: WarpingParameter: are the checks (w/nAlpha/nPsi...) ok?
      // TODO: move the checks somewhere else?

      // SPS is sampled in nAlpha / nPsi steps
      if (w % nAlpha != 0)
	throw SIMDException("WarpingParameter::WarpingParameter",
			    "width w needs to be a multiple of nAlpha");
      if (w % nPsi != 0)
	throw SIMDException("WarpingParameter::WarpingParameter",
			    "width w needs to be a multiple of nPsi");
      // each block in the shuffled SPS contains nAlpha elements which
      // have to be read in parallel
      if (nAlpha % simd_sps_elems != 0)
	throw SIMDException("WarpingCompound::WarpingCompound",
			    "nAlpha needs to be a multiple of simd_sps_elems");

      // compass computation is done via horizontal instructions, so we need
      // to process simd_sps_elems in parallel, but this we only test if the
      // compass is used

      createModuloTables();
      createParameterArrays();
      // create shuffle table (after modulo-tables!)
      createShuffle();
    }

    virtual ~WarpingParameter()
    {
      destroyModuloTables();
      destroyParameterArrays();
      destroyShuffle();
    }

    // ------------------------------------------------------------------------

  protected:

    void createModuloTables()
    {
      // for our special order in dest (theta, delta) we need this modulo table
      // modulo2w = new int [2*w];
      // we need to multiply by wSPS to skip rows in the scale plane
      // for (int m = 0; m < w; m++) modulo2w[m] = m * wSPS;
      // shift to the middle of the array so that we can handle negative indices
      // modulo2wpw = modulo2w + w;
      // duplicate
      // memcpy(modulo2wpw, modulo2w, w * sizeof(int));
      // new version 20. Nov 15 (rm), 23. Nov 15 (rm): more space
      modulo2w = new int [6 * w];
      for (int m = 0; m < 6 * w; m++) modulo2w[m] = (m % w) * wSPS;
      modulo2wpw = modulo2w + 3 * w;
      // just plain modulo w table:
      // modulo = new int [2*w];
      // for (int i = 0; i < 2 * w; i++) modulo[i] = i % w;
      // modulopw = modulo + w;
      // new version 20. Nov 15 (rm), 23. Nov 15 (rm): more space
      modulo = new int [6 * w];
      for (int i = 0; i < 6 * w; i++) modulo[i] = i % w;
      modulopw = modulo + 3 * w;
      // modulo nAlpha (triple copy), new version 20. Nov 15 (rm)
      modNAlpha = new int [3*nAlpha];
      for (int i = 0; i < 3 * nAlpha; i++) modNAlpha[i] = i % nAlpha;
      // memcpy(modNAlpha + nAlpha, modNAlpha, nAlpha * sizeof(int));
      // memcpy(modNAlpha + 2 * nAlpha, modNAlpha, nAlpha * sizeof(int));
    }

    void destroyModuloTables()
    {
      if (modulo2w) {
	delete [] modulo2w;
	modulo2w = 0;
	modulo2wpw = 0;
      }
      if (modulo) {
	delete [] modulo;
	modulo = 0;
	modulopw = 0;
      }
      if (modNAlpha) {
	delete [] modNAlpha;
	modNAlpha = 0;
      }
    }

    void createShuffle()
    {
      // resorting of horizontal index for SSE2 implementation of search phase
      shuffle = new int [w];
      // forward shuffle
      for (int j = 0; j < w; j++)
	// we need simd_sps_elems words more in each section (blkSize)
	shuffle[j] = (j % stepAlpha) * blkSize + j / stepAlpha;
      // 30. Jul 11 (rm): combined shuffle-modulo table
      // 10. Oct 14 (rm): shuffleModuloPW
      // 20. Nov 15 (rm), 23. Nov 15 (rm): more space
      shuffleModulo = new int [6 * w];
      for (int j = 0; j < 6 * w; j++)
	shuffleModulo[j] = shuffle[modulo[j]];
      shuffleModuloPW = shuffleModulo + 3 * w;
      //  1. Oct 14 (rm): unshuffle
      unshuffle = new int [wSPS];
      for (int js = 0; js < wSPS; js++)
	// TODO: unshuffle: is % w necessary?
	unshuffle[js] = ((js % blkSize) * stepAlpha + js / blkSize) % w;
    }

    void destroyShuffle()
    {
      if (shuffle) {
	delete [] shuffle;
	shuffle = 0;
      }
      if (shuffleModulo) {
	delete [] shuffleModulo;
	shuffleModulo = 0;
	shuffleModuloPW = 0;
      }
      if (unshuffle) {
	delete [] unshuffle;
	unshuffle = 0;
      }
    }

    void createParameterArrays()
    {
      // resize tables
      jAlphaVec = new int [nAlpha];
      alphaVec = new double [nAlpha];
      jPsiVec = new int [nPsi];
      psiVec = new double [nPsi];
      // initialize tables (modified from WarpTableSet)
      for (int iAlpha = 0; iAlpha < nAlpha; iAlpha++) {
	jAlphaVec[iAlpha] = iAlpha * stepAlpha;
	alphaVec[iAlpha] = jAlphaVec[iAlpha] * dTheta;
      }
      for (int iPsi = 0; iPsi < nPsi; iPsi++) {
	jPsiVec[iPsi] = iPsi * stepPsi;
	psiVec[iPsi] = jPsiVec[iPsi] * dTheta;
      }
    }

    void destroyParameterArrays()
    {
      if (jAlphaVec) {
	delete [] jAlphaVec;
	jAlphaVec = 0;
      }
      if (alphaVec) {
	delete [] alphaVec;
	alphaVec = 0;
      }
      if (jPsiVec) {
	delete [] jPsiVec;
	jPsiVec = 0;
      }
      if (psiVec) {
	delete [] psiVec;
	psiVec = 0;
      }
    }

  };

}

#endif
