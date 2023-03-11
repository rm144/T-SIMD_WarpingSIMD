// ===========================================================================
//
// SIMDException.H --
// exception class for SIMD code (with location and error string)
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

#ifndef SIMD_EXCEPTION_H_
#define SIMD_EXCEPTION_H_

#include <exception>
#include <string>

namespace simd {

class SIMDException : public std::exception
{
public:
  const std::string exc, loc, err, msg;

  SIMDException(const std::string &loc, const std::string &err)
    : exception(), exc("SIMDException"), loc(loc), err(err),
      msg(exc + " at " + loc + ": " + err)
  {}

  virtual ~SIMDException() throw() {}

  virtual const char *what() const throw() { return msg.c_str(); }
};

} // namespace simd

#endif
