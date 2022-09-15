// ===========================================================================
// 
// MinWarpingSearch.H --
// classes for full and partial MinWarping search algorithms
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

#ifndef _MINWARPING_SEARCH_H_
#define _MINWARPING_SEARCH_H_

#include <cmath>
#include "MinWarpingSearchRange.H"
#include "MinWarpingTemplate.H"
#include "WarpingSPS.H"
#include "MinWarpingMatch.H"

namespace ns_simd {

  // ===========================================================================
  // MinWarping full search algorithms
  // ===========================================================================

  template <typename SPSType, 
	    typename MatchType,
	    int SIMD_WIDTH, 
	    int SIMD_ALIGN>
  class MinWarpingFull
  {
  public:
    virtual ~MinWarpingFull() {}
    virtual void 
    search(const MinWarpingTemplate &minTemplate,
	   const WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps, 
	   MinWarpingMatch<MatchType,SIMD_WIDTH,SIMD_ALIGN> &match) = 0;
  };


  // ===========================================================================
  // MinWarping partial search algorithms
  // ===========================================================================

  template <typename SPSType,
	    typename MatchType,
	    int SIMD_WIDTH, 
	    int SIMD_ALIGN>
  class MinWarpingPartial
  {
  public:
    virtual ~MinWarpingPartial() {}
    virtual void 
    search(const MinWarpingTemplate &minTemplate,
	   const WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
	   const MinWarpingSearchRange<SPSType,SIMD_WIDTH,SIMD_ALIGN> &range,
	   MinWarpingMatch<MatchType,SIMD_WIDTH,SIMD_ALIGN> &match,
	   int jAlphaOff = 0, int jPsiOff = 0) = 0;
  };

}

#endif
