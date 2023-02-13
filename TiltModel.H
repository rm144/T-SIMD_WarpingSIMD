// ===========================================================================
//
// TiltModel.H --
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

#ifndef TILT_MODEL_H
#define TILT_MODEL_H

// Math functions for tilt model conversions
#include <cmath>
// Generic STL tuple and tuple utility
#include <tuple>
// Used to overload stream insertion operator
#include <ostream>

// Strong type angular units
#include "Angle.H"

// Homogenous transformations
#include "TransformMatrix.H"

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {}
// Export a short, flat name to this namespace
namespace tilt = ns_simd::tilt;

// Panoramic image tilt types and functions namespace
namespace ns_simd::tilt {
/**
 * @brief Roll-Pitch representation tilt model (X-Y-Axis tilt)
 * @tparam Unit Angular unit of tilt angles
 */
template <class Unit>
struct RollPitch
{
  // Tilt angles about the X- and Y-axis in [Unit]
  Angle<Unit> delta_X, delta_Y;

public:
  /**
   * @brief Initializes zero tilt
   */
  RollPitch() : delta_X {0}, delta_Y {0} {}

  /**
   * @brief Initializes tilt from individual angles as base type double
   * @param delta_X Tilt about the X-Axis in [Unit]
   * @param delta_Y Tilt about the Y-Axis in [Unit]
   */
  RollPitch(const double delta_X, const double delta_Y)
    : delta_X {delta_X}, delta_Y {delta_Y}
  {}

  /**
   * @brief Initializes tilt from individual angles
   * @param delta_X Tilt about the X-Axis in [Unit]
   * @param delta_Y Tilt about the Y-Axis in [Unit]
   */
  RollPitch(const Angle<Unit> &delta_X, const Angle<Unit> &delta_Y)
    : delta_X {delta_X}, delta_Y {delta_Y}
  {}

  /**
   * @brief Initializes tilt from std::tuple of angles as base type double
   * @param angles Tilt about X- and Y-Axis in [Unit] packed in
   * std::tuple {delta_X, delta_Y}
   */
  explicit RollPitch(const std::tuple<double, double> &angles)
    : delta_X {std::get<0>(angles)}, delta_Y {std::get<1>(angles)}
  {}

  /**
   * @brief Initializes tilt from std::tuple of angles
   * @param angles Tilt about X- and Y-Axis in [Unit] packed in
   * std::tuple {delta_X, delta_Y}
   */
  explicit RollPitch(const std::tuple<Angle<Unit>, Angle<Unit>> &angles)
    : delta_X {std::get<0>(angles)}, delta_Y {std::get<1>(angles)}
  {}

  /**
   * @brief Converts to Roll-Pitch representation of other angular unit
   * @tparam Target Angular unit to convert to
   */
  template <class Target>
  operator RollPitch<Target>() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // NOTE: Relies on implicit angular unit conversion of Angle
    return RollPitch<Target>(delta_X, delta_Y);
  }

  /**
   * @brief Converts to other tilt model (incl. Unit conversion)
   * @tparam Target Tilt model to convert to
   * @tparam TargetUnit Angular unit to use with Target tilt model
   */
  template <template <class> class Target, class TargetUnit>
  operator Target<TargetUnit>() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // Wrap with simple unit conversion of Target model
    return Target<TargetUnit> {// Conversion must be defined for Radian
                               RollPitch<Radian>((*this))};
  }

  /**
   * @brief Converts tilt model to homogenous transformation matrix
   */
  operator transform::Matrix() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // Make matrix-matrix product operator available
    using transform::operator*;
    // Multiply X- and Y-Axis transformations
    return transform::Rx((double) delta_X) * transform::Ry((double) delta_Y);
  }

  /**
   * @brief Stream insertion operator overload to print tilt model
   * @param lhs Left-hand side of the operator: Stream to insert the tilt
   * model angles into
   * @param rhs Right-hand side: The tilt model to insert into the stream
   * @return Reference to the output stream (left-hand side) after the
   * insertion operation.
   */
  friend std::ostream &operator<<(std::ostream &lhs, const RollPitch &rhs)
  {
    return lhs << rhs.delta_X << " " << rhs.delta_Y;
  }
};

/**
 * @brief Axis-Angle representation tilt model (axis [actually angle to axis]
 * + tilt about axis tilt)
 * @tparam Unit Angular unit of tilt angles
 */
template <class Unit>
struct AxisAngle
{
  // Angle to tilt axis and tilt angle about that axis in [Unit]
  Angle<Unit> theta_R, phi;

public:
  /**
   * @brief Initializes zero tilt
   */
  AxisAngle() : theta_R {0}, phi {0} {}

  /**
   * @brief Initializes tilt from individual angles as base type double
   * @param theta_R Angle to tilt axis [Unit]
   * @param phi Angle of tilt about the tile axis [Unit]
   */
  AxisAngle(const double theta_R, const double phi)
    : theta_R {theta_R}, phi {phi}
  {}

  /**
   * @brief Initializes tilt from individual angles
   * @param theta_R Angle to tilt axis [Unit]
   * @param phi Angle of tilt about the tile axis [Unit]
   */
  AxisAngle(const Angle<Unit> &theta_R, const Angle<Unit> &phi)
    : theta_R {theta_R}, phi {phi}
  {}

  /**
   * @brief Initializes tilt from std::tuple of angles as base type double
   * @param angles Angle to tilt axis and tilt angle in [Unit] packed in
   * std::tuple {theta_R, phi}
   */
  explicit AxisAngle(const std::tuple<double, double> &angles)
    : theta_R {std::get<0>(angles)}, phi {std::get<1>(angles)}
  {}

  /**
   * @brief Initializes tilt from std::tuple of angles
   * @param angles Angle to tilt axis and tilt angle in [Unit] packed in
   * std::tuple {theta_R, phi}
   */
  explicit AxisAngle(const std::tuple<Angle<Unit>, Angle<Unit>> &angles)
    : theta_R {std::get<0>(angles)}, phi {std::get<1>(angles)}
  {}

  /**
   * @brief Converts to Axis-Angle representation of other angular unit
   * @tparam Target Angular unit to convert to
   */
  template <class Target>
  operator AxisAngle<Target>() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // NOTE: Relies on implicit angular unit conversion of Angle
    return AxisAngle<Target>(theta_R, phi);
  }

  /**
   * @brief Converts to other tilt model (incl. Unit conversion)
   * @tparam Target Tilt model to convert to
   * @tparam TargetUnit Angular unit to use with Target tilt model
   */
  template <template <class> class Target, class TargetUnit>
  operator Target<TargetUnit>() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // Wrap with simple unit conversion of Target model
    return Target<TargetUnit> {// Conversion must be defined for Radian
                               AxisAngle<Radian>((*this))};
  }

  /**
   * @brief Converts tilt model to homogenous transformation matrix
   */
  operator transform::Matrix() const
  { // NOLINT: This is an
    // intentional implicit conversion.
    // Make matrix-matrix product operator available
    using transform::operator*;
    // Multiply Z- and X-Axis transformations
    return transform::Rz((double) theta_R) * transform::Rx((double) phi);
  }

  /**
   * @brief Stream insertion operator overload to print tilt model
   * @param lhs Left-hand side of the operator: Stream to insert the tilt
   * model angles into
   * @param rhs Right-hand side: The tilt model to insert into the stream
   * @return Reference to the output stream (left-hand side) after the
   * insertion operation.
   */
  friend std::ostream &operator<<(std::ostream &lhs, const AxisAngle &rhs)
  {
    return lhs << rhs.theta_R << " " << rhs.phi;
  }
};

// Specialization of Roll-Pitch to Axis-Angle model conversion (given for
// Radian, others derived automatically)
template <>
template <>
RollPitch<Radian>::operator AxisAngle<Radian>() const
{
  return AxisAngle<Radian> {
    // Rotation Axis theta_R
    Angle<Radian>(atan2((double) delta_Y, (double) delta_X)),
    // Rotation angle phi
    Angle<Radian>(acos(cos((double) delta_X) * cos((double) delta_Y)))};
}

// Specialization of Axis-Angle to Roll-Pitch model conversion (given for
// Radian, others derived automatically)
template <>
template <>
AxisAngle<Radian>::operator RollPitch<Radian>() const
{
  return RollPitch<Radian> {
    // Roll angle delta_X
    Angle<Radian>((double) phi * cos((double) theta_R)),
    // Pitch angle delta_Y
    Angle<Radian>((double) phi * sin((double) theta_R))};
}
} // namespace ns_simd::tilt

#endif /*TILT_MODEL_H*/
