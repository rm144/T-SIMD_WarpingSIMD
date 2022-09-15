// ===========================================================================
// 
// WarpingSearchCollection.H --
// collection of search algorithms for warping
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

#ifndef _WARPING_SEARCH_COLLECTION_H_
#define _WARPING_SEARCH_COLLECTION_H_

#include "WarpingSearch.H"

namespace ns_simd {

  // ###########################################################################
  // full matching
  // ###########################################################################

  // ===========================================================================
  // PARX
  // ===========================================================================

  // PARX denotes search order: psi, alpha, rho, x
  template <typename SPSType, 
	    typename MatchType,
	    int SIMD_WIDTH, 
	    int SIMD_ALIGN>
  class WarpingFull_PARX : 
    public WarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
  {
    virtual void search
    (const WarpingTemplate &warpTemplate,
     const WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps, 
     MinWarpingMatch<MatchType,SIMD_WIDTH,SIMD_ALIGN> &match)
    {
      // we assume that sizeof(SPSType) <= sizeof(MatchType)
      // (since we add up SPSType elements in MatchType sum)
      // (just to provide the user with a clear compiler message)
      typeSizeLEQ<SPSType,MatchType>();
      
      // 30. Sep 14 (rm): local copies (somewhat faster, surprisingly)
      //  1. Jun 15 (rm): now const
      const int w = warpTemplate.w, nRho = warpTemplate.nRho;
      const int nPsi = sps.param.nPsi, nAlpha = sps.param.nAlpha;
      int *const shuffleModulo = sps.param.shuffleModulo;
      int *const modulo2wpw = sps.param.modulo2wpw;
      int *const jPsiVec = sps.param.jPsiVec;
      int *const jAlphaVec = sps.param.jAlphaVec;
      MatchType *const _match = match.match.data;
      SPSType *const _stack = sps.stack.data;
      const int simd_sps_elems = sps.param.simd_sps_elems;
      const MatchType maxMatch = SIMDTypeInfo<MatchType>::max();
      // faster if not defined inside loops?
      SIMDVec<SPSType,SIMD_WIDTH> v;
      SIMDVecs<NumSIMDVecs<MatchType,SPSType>::out, MatchType, SIMD_WIDTH>
	vM, minSum, sum;
      
      // we write simd_sps_elems results to match array in one go
      if (nAlpha % simd_sps_elems != 0)
	throw SIMDException("WarpingFull_PARX::search",
			    "nAlpha has to be a multiple of simd_sps_elems");
      // here the match array is accumulated
      // TODO: entire clear in the beginning or split (when needed) - which 
      // TODO: is faster?
      match.clear();
      MatchType *matchp = _match;
      // go through all psi parameters
      for (int iPsi = 0; iPsi < nPsi; iPsi++) {
	// precompute part of (jY - jPsi) % w (times wSPS)
	int *modulo2wpwMinusPsi = modulo2wpw - jPsiVec[iPsi];
	// go through all alpha parameters
	for (int iAlpha = 0; iAlpha < nAlpha; iAlpha += simd_sps_elems) {
	  // precompute part of shuffled (jX + jAlpha) % w
	  int *shuffleModuloPlusAlpha = shuffleModulo + jAlphaVec[iAlpha];
	  // minimal sum over all rhos
	  set1(minSum, maxMatch);
	  // go through all rho parameters
	  for (int iRho = 0; iRho < nRho; iRho++) {
	    // initialize sum to zero
	    setzero(sum);
	    // go through all x values
	    for (int jX = 0; jX < w; jX++) {
	      TemplateSPOffset sp = warpTemplate.tmplt[iRho][jX].sp;
	      TemplateJY jY = warpTemplate.tmplt[iRho][jX].jy;
	      // compute address in scale-plane stack
	      SPSType *vp = _stack	     // sps
		+ sp			     // sp in sps
		+ shuffleModuloPlusAlpha[jX] // column in sp
		+ modulo2wpwMinusPsi[jY];    // element in column
	      // load (unaligned!)
	      v = loadu<SIMD_WIDTH>(vp);
	      // extend and sum up
	      sum = adds(sum, extend<MatchType>(v));
	    }
	    // update minimum
	    minSum = min(minSum, sum);
	  }
	  store(matchp, minSum);
	  matchp += simd_sps_elems;
	}
      }
    }
  };

  // ###########################################################################
  // warping curve (structure taken from full-PARX, but no SIMD code)
  // ###########################################################################

  template <typename SPSType, 
	    typename MatchType,
	    int SIMD_WIDTH, 
	    int SIMD_ALIGN>
  void warpingCurve(const WarpingTemplate &warpTemplate,
		    const WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
		    int iAlpha, int iPsi,
		    std::vector<int> &jDeltaVec,
		    std::vector<int> &planeIndexVec,
		    int &iRhoMin,
		    MatchType &matchSumMin)
  {
    // 30. Sep 14 (rm): local copies (somewhat faster, surprisingly)
    //  1. Jun 15 (rm): now const
    const int w = warpTemplate.w, wSPS = warpTemplate.wSPS;
    const int nRho = warpTemplate.nRho;
    int *const modulo = sps.param.modulo;
    int *const shuffleModulo = sps.param.shuffleModulo;
    int *const modulopw = sps.param.modulopw;
    SPSType *const _stack = sps.stack.data;
    const MatchType maxMatch = SIMDTypeInfo<MatchType>::max();
    // resize result vectors
    jDeltaVec.resize(w);
    planeIndexVec.resize(w);
    // jAlpha, jPsi:
    const int jAlpha = sps.param.jAlphaVec[iAlpha];
    const int jPsi = sps.param.jPsiVec[iPsi];
    // minimal values
    matchSumMin = maxMatch;
    iRhoMin = -1;
    // go through all rho values
    for (int iRho = 0; iRho < nRho; iRho++) {
      // initialize sum to zero
      MatchType sum = 0;
      // go through all x values
      for (int jX = 0; jX < w; jX++) {
	TemplateSPOffset sp = warpTemplate.tmplt[iRho][jX].sp;
	TemplateJY jY = warpTemplate.tmplt[iRho][jX].jy;
	// compute address in scale-plane stack
	SPSType *vp = _stack		// sps
	  + sp				// sp in sps
	  + shuffleModulo[jX + jAlpha]	// column in sp
	  + wSPS * modulopw[(int) jY - jPsi];	// element in column
	// extend and sum up
	sum += MatchType(*vp);
      }
      // update minimum
      if (sum < matchSumMin) {
	matchSumMin = sum;
	iRhoMin = iRho;
      }
    }
    // this should never happen
    assert(iRhoMin != -1);
    // copy row from warp-template
    for (int jX = 0; jX < w; jX++) {
      int jTheta = modulo[jX + jAlpha];
      jDeltaVec[jTheta] = 
	modulopw[(int) warpTemplate.tmplt[iRhoMin][jX].jy - jPsi];
      planeIndexVec[jTheta] = warpTemplate.tmplt[iRhoMin][jX].sp / (wSPS * w);
    }
  }

  // ###########################################################################
  // partial matching
  // ###########################################################################

  // ===========================================================================
  // PARX
  // ===========================================================================

  // PARX denotes search order: psi, alpha, rho, x
  template <typename SPSType, 
	    typename MatchType,
	    int SIMD_WIDTH, 
	    int SIMD_ALIGN>
  class WarpingPartial_PARX : 
    public WarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
    virtual void search
    (const WarpingTemplate &warpTemplate,
     const WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
     const MinWarpingSearchRange<SPSType,SIMD_WIDTH,SIMD_ALIGN> &range,
     MinWarpingMatch<MatchType,SIMD_WIDTH,SIMD_ALIGN> &match,
     int jAlphaOff, int jPsiOff)
    {
      // we assume that sizeof(SPSType) <= sizeof(MatchType)
      // (since we add up SPSType elements in MatchType sum)
      // (just to provide the user with a clear compiler message)
      typeSizeLEQ<SPSType,MatchType>();

      //  1. Oct 14 (rm): local copies (sometimes faster, not here though)
      //  3. Jun 15 (rm): now const
      const int w = warpTemplate.w, nRho = warpTemplate.nRho;
      const int nAlpha = sps.param.nAlpha;
      int *const shuffleModuloPW = sps.param.shuffleModuloPW;
      int *const modulo2wpw = sps.param.modulo2wpw;
      int *const jAlphaVec = sps.param.jAlphaVec;
      int *const jPsiVec = sps.param.jPsiVec;
      MatchType *const _match = match.match.data;
      SPSType *const _stack = sps.stack.data;
      const int simd_sps_elems = sps.param.simd_sps_elems;
      const MatchType maxMatch = SIMDTypeInfo<MatchType>::max();
      // faster if not defined inside loops
      SIMDVec<SPSType,SIMD_WIDTH> v;
      SIMDVecs<NumSIMDVecs<MatchType,SPSType>::out, MatchType, SIMD_WIDTH>
	zero, vM, minSum, sum;

      setzero(zero);
      // we write simd_sps_elems results to match array in one go
      if (nAlpha % simd_sps_elems != 0)
	throw SIMDException("WarpingPartial_PARX::search",
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
	int iPsi = range.psiIndices[kPsi];
	MatchType *matchRow = _match + iPsi * nAlpha;
	// go through all alpha parameters
	// note that alphaIndices contain indices with a spacing of
	// simd_sps_elems!
	for (unsigned int kAlpha = 0; kAlpha < range.alphaIndices.size(); 
	     kAlpha++) {
	  int iAlpha = range.alphaIndices[kAlpha];
	  // set consecutive values to 0
	  store(matchRow + iAlpha, zero);
	}
      }
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
	// go through all alpha parameters
	// note that alphaIndices contain indices with a spacing of 8!
	for (unsigned int kAlpha = 0; kAlpha < range.alphaIndices.size(); 
	     kAlpha++) {
	  // get index
	  int iAlpha = range.alphaIndices[kAlpha];
	  // precompute part of shuffled (jX + jAlpha) % w
	  // 20. Nov 15 (rm): shuffleModuloPW instead of shuffleModulo,
	  // offset added
	  int *shuffleModuloPlusAlpha 
	    = shuffleModuloPW + jAlphaVec[iAlpha] + jAlphaOff;
	  // minimal sum over all rhos
	  set1(minSum, maxMatch);
	  // go through all rho parameters
	  for (int iRho = 0; iRho < nRho; iRho++) {
	    // initialize sum to zero
	    setzero(sum);
	    // go through all x values
	    for (int jX = 0; jX < w; jX++) {
	      TemplateSPOffset sp = warpTemplate.tmplt[iRho][jX].sp;
	      TemplateJY jY = warpTemplate.tmplt[iRho][jX].jy;
	      // compute address in scale-plane stack
	      SPSType *vp = _stack	     // sps
		+ sp			     // sp in sps
		+ shuffleModuloPlusAlpha[jX] // column in sp
		+ modulo2wpwMinusPsi[jY];    // element in column
	      // load (unaligned!)
	      v = loadu<SIMD_WIDTH>(vp);
	      // extend and sum up
	      sum = adds(sum, extend<MatchType>(v));
	    }
	    // update minimum
	    minSum = min(minSum, sum);
	  }
	  MatchType *matchp = matchRow + iAlpha;
	  store(matchp, minSum);
	}
      }
    }
  };

} // namespace
    
#endif
