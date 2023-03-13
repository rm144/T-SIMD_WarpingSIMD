# T-SIMD

\tableofcontents

# Introduction

**T-SIMD** is a low-level C++ template SIMD library which wraps built-in vector data types and built-in vector intrincics in template classes and template functions or overloaded functions, respectively.

Templates parameters are the element data type of the vectors and the vector width in bytes (e.g. 16 for SSE* and NEON, 32 for AVX/AVX2). This makes it possible to flexibly change the data type and the vector instruction set for entire portions of the code. Moreover, many implementation details at the intrinsics level are hidden by **T-SIMD**. SSE*, AVX/AVX2, AVX-512, and ARM NEON vector instruction sets are currently supported. Please note the [differences on different architectures](@ref differences-on-different-architectures).

# Compiler and C++ Standard Support/Requirements

**T-SIMD** requires at least C++11.

**T-SIMD** has been tested with gcc and clang on Linux.

# Quick Start

First, download the latest release from <https://www.ti.uni-bielefeld.de/html/people/moeller/tsimd_warpingsimd.html> (section "Software Download").
Either download the single header library version `tsimd_sh.H` or download the full archive, unpack it and copy the contents of the folder starting with "CODE" into your project.

To use the **T-SIMD** library, simply include the single-header file `tsimd_sh.H` or the `tsimd.H` header from the archive.

# Example

The following example shows how to use **T-SIMD** to compute the sum of two vectors
and print the result to the standard output.

```cpp
#include "tsimd.H" // or "tsimd_sh.H" if you use the single-header version
#include <stdio.h>

int main() {
  simd::Vec<simd::Float> a = simd::setzero<simd::Float>();
  simd::Vec<simd::Float> b = simd::setzero<simd::Float>();
  simd::Vec<simd::Float> sum = simd::add(a, b);
  simd::print("%g ", c);
  puts("");
  return 0;
}
```

# Preprocessor Definitions

The following preprocessor definitions can be defined to change the behavior of **T-SIMD** (**WARNING**: These definitions **must** be defined before including `tsimd.H` or any other **T-SIMD** header file):

* `SIMD_ALIGN_CHK` : If this macro is defined, **T-SIMD** will check whether the alignment of the data pointers passed to the functions is correct.
* `MAX_SIMD_WIDTH` : This macro can be used to limit the maximum vector width that **T-SIMD** will use. Must be an integer of at least 16.
* `SIMDVEC_SANDBOX` : If this macro is defined, **T-SIMD** functions will not execute their intended operations, but instead will print what function was called to the standard output. This is useful for debugging purposes.
# Differences on Different Architectures {#differences-on-different-architectures}

**T-SIMD** tries to behave as similar as possible on different architectures. However, there are some differences that are due to the different instruction sets and the different implementations of the intrinsics. The following list gives an overview of the differences:

* Denormalized floating-point numbers are flushed to zero on ARMv7. This is not the case on all other supported architectures.
* The `simd::div` and `simd::sqrt` functions are approximations on ARM in contrast to the exact results on all other supported architectures.
* Converting from floating-point to integer values results in slightly different results on ARM.

There may be other slight differences that are not listed here.

# Author and Contributors

**T-SIMD** was originally written by [Ralf Möller](http://www.ti.uni-bielefeld.de/html/people/moeller/) (moeller@ti.uni-bielefeld.de).

Other contributors are: Jonas Keller, Markus Vieth, Adam Marschall, Jan-Lukas Wolf, Lukas Schiermeier, and Moritz Breipohl.

Their specific contributions are noted in the source code.

# License

**T-SIMD** is licensed under the license contained in the files `LICENSE` and `LICENSE.doc` .

# Contact and Bug Reports

Please report any bugs or errors found to [Ralf Möller](http://www.ti.uni-bielefeld.de/html/people/moeller/) (moeller@ti.uni-bielefeld.de).
