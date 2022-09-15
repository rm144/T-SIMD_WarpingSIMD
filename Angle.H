#ifndef ANGLE_H
#define ANGLE_H

// Math functions and constants (e.g. M_PI)
#include <cmath>
// Used to overload stream insertion operator
#include <ostream>

/**
 * @brief Radian angular Unit tag to use with strong type Angle<Unit>
 */
struct Radian {};

/**
 * @brief Degree angular Unit tag to use with strong type Angle<Unit>
 */
struct Degree {};

/**
 * @brief Arc Second angular Unit tag to use with strong type Angle<Unit>
 */
struct ArcSec {};

/**
 * @brief Arc Minute angular Unit tag to use with strong type Angle<Unit>
 */
struct ArcMin {};

/**
 * @brief Angular unit conversion base template
 * @tparam Unit Unit of the input angle
 * @tparam Target Desired unit of the output angle
 * @tparam Type Base type of the angle
 * @param x Angle value in [Unit]
 * @return Converted angle value of same base Type in [Target]
 */
template<class Unit, class Target, class Type = double>
    inline Type convert(const Type x, Unit = Unit{}, Target = Target{}) {
        // Assume all Units to be convertible to/from Degree. Convert all units
        // vai conversion chain Unit->Degree->Target. Thus it suffices to
        // specify conversion to/from Degree to enable arbitrary conversion for
        // new/custom angular units.
        return convert(convert(x, Unit{}, Degree{}), Degree{}, Target{});
    }

// Special case to avoid conversion of the same Unit
template<class Type = double, class Unit>
    inline Type convert(const Type x, Unit = Unit{}, Unit = Unit{}) {
        // No conversion necessary
        return x;
    }

/**
 * @brief Macro to instantiate generic angular unit conversions using a
 * converting factor.
 */
#define CONVERT(Unit, Target, Factor)                                      \
    template<class Type> inline Type convert(const Type x, Unit, Target) { \
        return static_cast<Type>((Factor)) * x;                            \
    }

// Converts Radian to Degree
CONVERT(Radian, Degree, (180.0 / M_PI))

// Converts Degree to Radian
CONVERT(Degree, Radian, (M_PI / 180.0))

// Converts ArcSec to Degree
CONVERT(ArcSec, Degree, (1.0 / 3600.0))

// Converts Degree to ArgSec
CONVERT(Degree, ArcSec, (1.0 * 3600.0))

// Converts ArcMin to Degree
CONVERT(ArcMin, Degree, (1.0 / 60.0))

// Converts Degree to ArgMin
CONVERT(Degree, ArcMin, (1.0 * 60.0))

/**
 * @brief Strong type augmenting the wrapped type by angular unit and unit
 * conversions.
 * @tparam Unit Angular unit type selector
 * @tparam Type The wrapped base type
 */
template<class Unit, class Type = double>
    class Angle {
        // Need access to all other Angle types to implement conversions
        template<class, class>
            friend
            class Angle;

        // The wrapped angle value of the base Type
        Type x;

    public:
        /**
         * @brief Default constructs a strong type angle
         */
        Angle() = default;

        /**
         * @brief Constructs a strong type angle given a value of the base
         * Type.
         * @param x Angle value in [Unit]
         */
        explicit Angle(const Type x)
            : x{x} {}

        /**
         * @brief Converts back to the base type
         * @return The angle value as the base Type
         */
        explicit operator Type() const {
            return x;
        }

        /**
         * @brief Accesses reference to the wrapped base Type
         * @return A base Type reference to the wrapped angle value
         */
        explicit operator Type &() {
            return x;
        }

        /**
         * @brief Casts to value of arbitrary other type
         * @tparam Other Type to cast to
         * @return Angle value type-casted to Other
         */
        template<class Other>
            explicit operator Other() const {
                return static_cast<Other>(x);
            }

        /**
         * @brief Converts to strong typed angle of another Unit
         */
        template<class Target>
            operator Angle<Target, Type>() const { // NOLINT: This is an
                // intentional implicit conversion.
                return Angle<Target, Type>(convert<Unit, Target>(x));
            }

        /**
         * @brief Converts to strong typed angle of another base Type
         */
        template<class Target>
            operator Angle<Unit, Target>() const { // NOLINT: This is an
                // intentional implicit conversion.
                return Angle<Unit, Target>(static_cast<Target>(x));
            }

        /**
         * @brief Converts to strong typed angle of another base Type and Unit
         * @note Unit conversion is done first, then the converted angle is type
         * casted.
         */
        template<class TargetUnit, class TargetType>
            operator Angle<TargetUnit, TargetType>() const { // NOLINT: This is
                // an intentional implicit conversion.
                return Angle<TargetUnit, TargetType>(
                    static_cast<TargetType>(convert<Unit, TargetUnit>(x))
                );
            }

        /**
         * @brief Adds angle of the same Unit and Type
         * @return Reference to the left-hand side (*this)
         */
        Angle &operator+=(const Angle &rhs) {
            // Add to the angle value and return reference to this
            return x += rhs.x, (*this);
        }

        /**
         * @brief Subtracts angle of the same Unit and Type
         * @return Reference to the left-hand side (*this)
         */
        Angle &operator-=(const Angle &rhs) {
            // Subtract to the angle value and return reference to this
            return x -= rhs.x, (*this);
        }

        /**
         * @brief Scales angle by scalar of base Type
         */
        Angle &operator*=(const Type &rhs) {
            // Multiply with the angle value and return reference to this
            return x *= rhs, (*this);
        }

        /**
         * @brief Divides angle by scalar of base Type
         */
        Angle &operator/=(const Type &rhs) {
            // Divide the angle value and return reference to this
            return x /= rhs, (*this);
        }

        /**
         * @brief Adds two angle of the same Unit and Type
         * @return The sum of the angles having the same Unit and Type
         */
        friend Angle operator+(Angle lhs, const Angle &rhs) {
            // Reuse compound assignment, return left-hand side after addition
            return lhs += rhs, lhs;
        }

        /**
         * @brief Subtracts two angle of the same Unit and Type
         * @return The difference of the angles having the same Unit and Type
         */
        friend Angle operator-(Angle lhs, const Angle &rhs) {
            // Reuse compound assignment, return left-hand side after
            // subtraction
            return lhs -= rhs, lhs;
        }

        /**
         * @brief Left multiplies angle with scalar of base Type
         */
        friend Angle operator*(const Type &lhs, Angle rhs) {
            // Reuse compound assignment, return right-hand side after
            // multiplication
            return rhs *= lhs, rhs;
        }

        /**
         * @brief Divides angle by scalar of base Type from the right
         */
        friend Angle operator/(Angle lhs, const Type &rhs) {
            // Reuse compound assignment, return left-hand side after
            // division
            return lhs /= rhs, lhs;
        }

        /**
         * @brief Unary plus operator on Angle<Unit>
         */
        Angle operator+() const {
            // Apply unary plus to wrapped value
            return Angle(+x);
        }

        /**
         * @brief Unary minus operator on Angle<Unit>
         */
        Angle operator-() const {
            // Apply unary minus to wrapped value
            return Angle(-x);
        }

        /**
         * @brief Compares two angles <
         */
        template<class Other, class OtherType>
            friend bool
            operator<(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x < ((Angle) rhs).x;
            }

        /**
         * @brief Compares two angles >
         */
        template<class Other, class OtherType>
            friend bool
            operator>(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x > ((Angle) rhs).x;
            }

        /**
         * @brief Compares two angles <=
         */
        template<class Other, class OtherType>
            friend bool
            operator<=(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x <= ((Angle) rhs).x;
            }

        /**
         * @brief Compares two angles >=
         */
        template<class Other, class OtherType>
            friend bool
            operator>=(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x >= ((Angle) rhs).x;
            }

        /**
         * @brief Compares two angles for equality
         */
        template<class Other, class OtherType>
            friend bool
            operator==(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x == ((Angle) rhs).x;
            }

        /**
         * @brief Compares two angles for inequality
         */
        template<class Other, class OtherType>
            friend bool
            operator!=(const Angle &lhs, const Angle<Other, OtherType> &rhs) {
                // Compare the wrapped values after converting right-hand side
                // to this (left-hand side) type
                return lhs.x != ((Angle) rhs).x;
            }

        /**
         * @brief Stream insertion operator overload to print angles
         * @param lhs Left-hand side of the operator: Stream to insert the angle
         * into
         * @param rhs Right-hand side: The angle to insert into the stream
         * @return Reference to the output stream (left-hand side) after the
         * insertion operation.
         */
        friend std::ostream &operator<<(std::ostream &lhs, const Angle &rhs) {
            return lhs << rhs.x;
        }

        /**
         * @brief Stream extraction operator overload to read angles
         * @param lhs Left-hand side of the operator: Stream to extract the
         * angle from
         * @param rhs Right-hand side: The angle to extract into
         * @return Reference to the input stream (left-hand side) after the
         * extraction operation.
         */
        friend std::istream &operator>>(std::istream &lhs, Angle &rhs) {
            return lhs >> rhs.x;
        }
    };

/**
 * @brief Absolute value of an angle
 */
template<class Unit, class Type>
    Angle<Unit, Type> fabs(const Angle<Unit, Type> angle) {
        return Angle<Unit, Type>(fabs((Type)angle));
    }

/**
 * @brief Constructs a strong typed angle in radian from long double literal
 * @param x Angle value in radian
 * @return Strong typed radian angle with base type long double of x
 */
Angle<Radian, long double> operator "" _rad(long double x) {
    return Angle<Radian, long double>(x);
}

/**
 * @brief Constructs a strong typed angle in degree from long double literal
 * @param x Angle value in degree
 * @return Strong typed degree angle with base type long double of x
 */
Angle<Degree, long double> operator "" _deg(long double x) {
    return Angle<Degree, long double>(x);
}

#endif /*ANGLE_H*/