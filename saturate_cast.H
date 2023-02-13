// ===========================================================================
//
// saturate_cast.H --
// Christoph Berganski
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

#ifndef SATURATE_CAST_H
#define SATURATE_CAST_H

// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"

// std::min and std::max
#include <algorithm>

// WarpingSIMD / SIMDImage namespace
namespace simd {
/**
 * @brief Type cast with saturation. This means: Values outside of the target
 * type range are clipped to the target type min or max value.
 */
template <typename Tout, typename Tin>
Tout saturate_cast(const Tin v)
{
  // Target Type min value in input type
  Tin min = (Tin) TypeInfo<Tout>::min();
  // Target Type max value in input type
  Tin max = (Tin) TypeInfo<Tout>::max();
  // Clip value to [min,max] before casting to Tout
  return (Tout) std::min(std::max(v, min), max);
}
} // namespace simd

#endif /*SATURATE_CAST_H*/
