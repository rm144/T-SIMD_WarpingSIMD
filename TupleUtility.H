// ===========================================================================
// 
// TupleUtility.H --
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

#ifndef TUPLE_UTILITY_H
#define TUPLE_UTILITY_H

// Tuple implementation of the STL
#include <tuple>
// Standard utility function (e.g. std::index_sequence)
#include <utility>

/**
 * @brief Custom operators operating on std::tuple
 */
namespace tuple_operators {
// Adds two tuples. Auxiliary template unpacking the index sequence
template<class Tuple, std::size_t... Is>
    Tuple add(const Tuple &t1, const Tuple &t2, std::index_sequence<Is...>) {
        return {(std::get<Is>(t1) + std::get<Is>(t2))...};
    }

// Adds two tuples.
template<class... Types>
    auto add(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return add(t1, t2, std::index_sequence_for<Types...>{});
    }

// Tuple addition operator overload.
template<class... Types>
    auto
    operator+(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return add(t1, t2);
    }

// Adds const Tuple &to left-hand side. Auxiliary template unpacking the index
// sequence
template<class Tuple, std::size_t... Is>
    Tuple &add_to(Tuple &t1, const Tuple &t2, std::index_sequence<Is...>) {
        return ((std::get<Is>(t1) += std::get<Is>(t2)), ..., t1);
    }

// Adds const Tuple &to left-hand side.
template<class... Types>
    auto &add_to(std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return add_to(t1, t2, std::index_sequence_for<Types...>{});
    }

// Tuple addition operator overload.
template<class... Types>
    auto operator+=(std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return add_to(t1, t2);
    }

// Subtracts two tuples. Auxiliary template unpacking the index sequence
template<class Tuple, std::size_t... Is>
    Tuple sub(const Tuple &t1, const Tuple &t2, std::index_sequence<Is...>) {
        return {(std::get<Is>(t1) - std::get<Is>(t2))...};
    }

// Subtracts two tuples.
template<class... Types>
    auto sub(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return sub(t1, t2, std::index_sequence_for<Types...>{});
    }

// Tuple subtraction operator overload.
template<class... Types>
    auto
    operator-(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return sub(t1, t2);
    }

// Multiplies tuple elements by scalar value. Auxiliary template unpacking
// the index sequence.
template<class Type, class Tuple, std::size_t... Is>
    Tuple mul(const Type s, const Tuple &t, std::index_sequence<Is...>) {
        return {(s * std::get<Is>(t))...};
    }

// Multiplies tuple elements by scalar value.
template<class Type, class... Types>
    auto mul(const Type s, const std::tuple<Types...> &t) {
        return mul(s, t, std::index_sequence_for<Types...>{});
    }

// Tuple scalar multiplication operator overload.
template<class Type, class... Types>
    auto operator*(const Type s, const std::tuple<Types...> &t) {
        return mul(s, t);
    }

// Divides tuple elements by scalar value. Auxiliary template unpacking
// the index sequence.
template<class Type, class Tuple, std::size_t... Is>
    Tuple div(const Tuple &t, const Type s, std::index_sequence<Is...>) {
        return {(std::get<Is>(t) / s)...};
    }

// Divides tuple elements by scalar value.
template<class Type, class... Types>
    auto div(const std::tuple<Types...> &t, const Type s) {
        return div(t, s, std::index_sequence_for<Types...>{});
    }

// Tuple scalar multiplication operator overload.
template<class Type, class... Types>
    auto operator/(const std::tuple<Types...> &t, const Type s) {
        return div(t, s);
    }

// Finds the maximum element in a tuple
template<class Tuple, std::size_t I>
    auto max_element(const Tuple &t, std::index_sequence<I>) {
        // Maximum of a one element tuple
        return std::get<I>(t);
    }

// Finds the maximum element in a tuple
template<class Tuple, std::size_t I, std::size_t... Is>
    auto max_element(const Tuple &t, std::index_sequence<I, Is...>) {
        // Get the maximum element in tail of the tuple
        auto max = max_element(t, std::index_sequence<Is...>{});
        // The maximum in the tuple is either the first or the max of the tail
        return std::get<I>(t) >= max ? std::get<I>(t) : max;
    }

// Finds the maximum element in a tuple
template<class Type>
    auto max_element(const std::tuple<Type> &t) {
        return std::get<Type>(t);
    }

// Finds the maximum element in a tuple
template<class... Types>
    auto max_element(const std::tuple<Types...> &t) {
        return max_element(t, std::index_sequence_for<Types...>{});
    }

// Elementwise absolute value of a tuple
template<class... Types>
    std::tuple<Types...> fabs(const std::tuple<Types...> &t) {
        return std::apply(
            [](Types... values) { return std::tuple{fabs(values)...}; }, t
        );
    }

// Elementwise maximum of two tuples. Auxiliary template unpacking the index
// sequence
template<class Tuple, std::size_t... Is>
    Tuple max(const Tuple &t1, const Tuple &t2, std::index_sequence<Is...>) {
        return {((std::get<Is>(t1) > std::get<Is>(t2)) ? std::get<Is>(t1)
                                                       : std::get<Is>(t2))...};
    }

// Elementwise maximum of two tuples.
template<class... Types>
    auto max(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return max(t1, t2, std::index_sequence_for<Types...>{});
    }

// Elementwise minimum of two tuples. Auxiliary template unpacking the index
// sequence
template<class Tuple, std::size_t... Is>
    Tuple min(const Tuple &t1, const Tuple &t2, std::index_sequence<Is...>) {
        return {((std::get<Is>(t1) < std::get<Is>(t2)) ? std::get<Is>(t1)
                                                       : std::get<Is>(t2))...};
    }

// Elementwise minimum of two tuples.
template<class... Types>
    auto min(const std::tuple<Types...> &t1, const std::tuple<Types...> &t2) {
        return min(t1, t2, std::index_sequence_for<Types...>{});
    }

// Product of all elements in a tuple. Auxiliary template unpacking the index
// sequence
template<class Tuple, std::size_t... Is>
    auto product(const Tuple &t, std::index_sequence<Is...>) {
        return (std::get<Is>(t) * ...);
    }

// Product of all elements in a tuple.
template<class... Types>
    auto product(const std::tuple<Types...> &t) {
        return product(t, std::index_sequence_for<Types...>{});
    }

// Prints tuple to std::ostream. Auxiliary template unpacking the index sequence
template<class Tuple, std::size_t... Is>
    std::ostream &
    out(std::ostream &stream, const Tuple &t, std::index_sequence<Is...>) {
        return (((stream << std::get<Is>(t) << ' '), ...), stream);
    }

// Prints tuple to std::ostream.
template<class... Types>
    std::ostream &out(std::ostream &stream, const std::tuple<Types...> &t) {
        return out(stream, t, std::index_sequence_for<Types...>{});
    }

// Tuple stream insertion operator overload.
template<class... Types>
    std::ostream &
    operator<<(std::ostream &stream, const std::tuple<Types...> &t) {
        return out(stream, t);
    }
}

#endif /*TUPLE_UTILITY_H*/
