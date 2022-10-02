// ===========================================================================
// 
// TransformMatrix.H --
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

#ifndef TRANSFORM_H
#define TRANSFORM_H

// Fixed size generic array
#include <array>
// Math functions: sin(), cos()
#include <cmath>

// Prior to C++17 there are no inline variables
#if __cplusplus >= 201703L
// If inline variables are available, it is probably better to actually use
// inline for global constants: The variable (constant) has external linkage and
// shares one address across all translation units.
#define GLOBAL_CONSTANT    inline const
#else
// If not available, just fall back to static having one instance per
// translation unit.
#define GLOBAL_CONSTANT    static const
#endif

/**
 * @brief Homogenous transformations (4d matrices)
 */
namespace transform {
/**
 * @brief Homogenous transformation represented as 4x4 matrix
 */
using Matrix = std::array<std::array<double, 4>, 4>;

/**
 * @brief Homogenous coordinate vector (4d vector)
 */
using Vector = std::array<double, 4>;

/**
 * @brief The identity transformation
 */
GLOBAL_CONSTANT Matrix Id = {{
// @formatter:off
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
}};
// @formatter:on

/**
 * @brief Basic X-Axis rotation transformation
 */
static inline Matrix Rx(const double x) {
    // @formatter:off
    return {{
        {1,     0,       0,  0},
        {0, cos(x), -sin(x), 0},
        {0, sin(x),  cos(x), 0},
        {0,     0,       0,  1}
    }};
    // @formatter:on
}

/**
 * @brief Basic Y-Axis rotation transformation
 */
static inline Matrix Ry(const double y) {
    // @formatter:off
    return {{
        { cos(y), 0, sin(y), 0},
        {     0,  1,     0,  0},
        {-sin(y), 0, cos(y), 0},
        {     0,  0,     0,  1}
    }};
    // @formatter:on
}

/**
 * @brief Basic Z-Axis rotation transformation
 */
static inline Matrix Rz(const double z) {
    // @formatter:off
    return {{
        {cos(z), -sin(z), 0, 0},
        {sin(z),  cos(z), 0, 0},
        {0,           0,  1, 0},
        {0,           0,  0, 1}
    }};
    // @formatter:on
}

/**
 * @brief Basic X-Axis translation transformation
 */
static inline Matrix Tx(const double x) {
    // @formatter:off
    return {{
        {1, 0, 0, x},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    // @formatter:on
}

/**
 * @brief Basic Y-Axis translation transformation
 */
static inline Matrix Ty(const double y) {
    // @formatter:off
    return {{
        {1, 0, 0, 0},
        {0, 1, 0, y},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    // @formatter:on
}

/**
 * @brief Basic Z-Axis translation transformation
 */
static inline Matrix Tz(const double z) {
    // @formatter:off
    return {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, z},
        {0, 0, 0, 1}
    }};
    // @formatter:on
}

/**
 * @brief Transformation Matrix-Matrix multiplication (composition)
 */
static inline Matrix operator*(const Matrix &A, const Matrix &B) {
    // Initialize empty 4x4 matrix
    // @formatter:off
    Matrix C = {{
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0}
    }};
    // @formatter:on
    // Compute product (accumulate in C)
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            for (std::size_t k = 0; k < 4; ++k)
                C[i][j] += A[i][k] * B[k][j];
    // Return accumulator matrix
    return C;
}

/**
 * @brief Transformation Matrix-Vector multiplication
 */
static inline Vector operator*(const Matrix &A, const Vector &x) {
    // Initialize empty 4x1/1x4 vector
    Vector y = {{0.0, 0.0, 0.0, 0.0}};
    // Compute product (accumulate in y)
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            y[i] += A[i][j] * x[j];
    // Return accumulator vector
    return y;
}
}

#endif /*TRANSFORM_H*/
