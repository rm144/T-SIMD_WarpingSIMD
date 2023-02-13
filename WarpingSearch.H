// ===========================================================================
//
// WarpingSearch.H --
// class for full and partial Warping search algorithms
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

#ifndef _WARPING_SEARCH_H_
#define _WARPING_SEARCH_H_

#include "MinWarpingMatch.H"       // uses MinWarpingMatch
#include "MinWarpingSearchRange.H" // uses MinWarpingSearchRange
#include "WarpingSPS.H"
#include "WarpingTemplate.H"

namespace ns_simd {

// ===========================================================================
// Warping full search algorithms
// ===========================================================================

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingFull
{
public:
  virtual ~WarpingFull() {}
  virtual void search(
    const WarpingTemplate &warpingTemplate,
    const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
    MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN> &match) = 0;
};

// ===========================================================================
// Warping partial search algorithms
// ===========================================================================

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingPartial
{
public:
  virtual ~WarpingPartial() {}
  virtual void search(
    const WarpingTemplate &warpingTemplate,
    const WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
    const MinWarpingSearchRange<SPSType, SIMD_WIDTH, SIMD_ALIGN> &range,
    MinWarpingMatch<MatchType, SIMD_WIDTH, SIMD_ALIGN> &match,
    int jAlphaOff = 0, int jPsiOff = 0) = 0;
};

} // namespace ns_simd

#endif
