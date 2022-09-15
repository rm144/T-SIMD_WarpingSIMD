#ifndef SATURATE_CAST_H
#define SATURATE_CAST_H

// SIMD typedefs and TypeInfo
#include "SIMDTypes.H"

// WarpingSIMD / SIMDImage namespace
namespace ns_simd {
/**
 * @brief Type cast with saturation. This means: Values outside of the target
 * type range are clipped to the target type min or max value.
 */
template<typename Tout, typename Tin>
    Tout saturate_cast(const Tin v) {
        // Target Type min value in input type
        Tin min = (Tin) SIMDTypeInfo<Tout>::min();
        // Target Type max value in input type
        Tin max = (Tin) SIMDTypeInfo<Tout>::max();
        // Clip value to [min,max] before casting to Tout
        return (Tout) std::min(std::max(v, min), max);
    }
}

#endif /*SATURATE_CAST_H*/
