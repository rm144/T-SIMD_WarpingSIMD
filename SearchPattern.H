// ===========================================================================
//
// SearchPattern.H --
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

#ifndef SEARCH_PATTERN_H
#define SEARCH_PATTERN_H

// Tuple implementation of the STL
#include <tuple>
// Standard utility function (e.g. std::index_sequence)
#include <utility>
// String streams (input/output from/to strings)
#include <sstream>

// Custom tuple operators
#include "TupleUtility.H"

// Namespace containing search strategies
namespace search {
/**
 * @brief N-dimensional discrete pattern search strategy
 */
template <class... Types>
class Pattern
{
  // Number of dimensions/axes
  static constexpr auto N = std::tuple_size_v<std::tuple<Types...>>;

  /**
   * @brief Pattern search state
   *
   * @note This is a hidden subclass of the pattern search strategy and
   * should solely be constructed via the Pattern::init<Callable> method.
   *
   * @tparam Type Type of the output value of the objective function
   */
  template <class Type>
  class State
  {
    // The wrapping Pattern class needs access
    friend Pattern;

    // The current and the previous center point and pattern width
    // along each axis
    std::tuple<Types...> argmin, argmin_prev, w;
    // The value at the current center point (argmin)
    Type min;

  public:
    /**
     * @brief Initializes new pattern search state
     *
     * @param argmin Initial minimum (center) point
     * @param w Initial widths of the pattern along each axis
     * @param min Initial minimum at the center point
     */
    State(const std::tuple<Types...> &argmin, const std::tuple<Types...> &w,
          const Type &min)
      : argmin {argmin}, argmin_prev {argmin}, w {w}, min {min}
    {}

    /**
     * @brief Returns minimum value and the corresponding argument.
     *
     * The minimum is at the pattern center and has the value of the
     * objective function at the center.
     *
     * @return std::pair of {min,argmin}
     */
    std::pair<Type, std::tuple<Types...>> min_argmin() const
    {
      return std::pair {min, argmin};
    }

    /**
     * @brief Constructs all pattern points
     * @return All pattern points as iterable std::array of
     * std::tuple
     */
    auto make_points() const
    {
      // Call auxiliary template with index sequence of number of
      // axes
      return make_points(std::index_sequence_for<Types...> {});
    }

    /**
     * @brief Produces a string representation of the search state
     * @return std::string representing the search state
     */
    [[nodiscard]] std::string str() const
    {
      // Stream insertion operator overload operating on
      // std::tuple
      using namespace tuple_operators;
      // Build string representation by inserting into output
      // string stream
      std::ostringstream repr;
      // Insert coordinates, width and minimum value into
      // string representation
      repr << argmin << " " << min << " " << w << " ";
      // Get the string out of the string stream
      return repr.str();
    }

    // Private state methods
  private:
    /**
     * @brief Makes a pattern point in direction along axis
     * @tparam Axis Pattern axis to move along
     * @param direction Direction to move (-1: down / +1: up)
     * @return Point of the pattern as std::tuple<Types...>
     */
    template <std::size_t Axis>
    auto make_point(const int direction) const
    {
      // Start at center point
      std::tuple<Types...> x = argmin;
      // Add width of axis in chosen direction to center
      std::get<Axis>(x) += direction * std::get<Axis>(w);
      // Return the point tuple
      return x;
    }

    /**
     * @brief Makes all (2N+1) pattern points along Axes in up and
     * down direction
     * @tparam Axes Pattern axes to make points on
     * @return std::array of all constructed points
     */
    template <std::size_t... Axes>
    auto make_points(std::index_sequence<Axes...>) const
    {
      // Construct array of 2*N+1 points. Each point is a
      // tuple of Types...
      return std::array<std::tuple<Types...>, 2 * N + 1> {
        // Add the current center point (the +1 point)
        argmin,
        // Make all points on given Axes in down (-1)
        // and up (+1) direction (the 2*N points)
        //  NOTE: Unpacks the Axes... pack twice into the
        //  array constructor - once per direction
        make_point<Axes>(-1)..., make_point<Axes>(+1)...};
    }
  };

  // The pattern center point, widths, thresholds and scale factors
  std::tuple<Types...> c {Types {0}...}, w {Types {0}...}, t {Types {1}...},
    s {Types {0}...};

public:
  /**
   * @brief Default initializes a pattern search
   * @note This pattern does immediately converge (as the initial width
   * is zero)
   */
  Pattern() = default;

  /**
   * @brief Constructs a pattern given a center point, widths,
   * thresholds and scale factors per axis.
   * @param c Pattern center point as tuple of Types...
   * @param w Pattern width per axis as tuple of Types...
   * @param t Pattern convergence threshold per axis as tuple of Types...
   * @param s Pattern shrinking scale factor per axis as tuple of Types...
   */
  Pattern(const std::tuple<Types...> &c, const std::tuple<Types...> &w,
          const std::tuple<Types...> &t, const std::tuple<Types...> &s)
    : c {c}, w {w}, t {t}, s {s}
  {}

  /**
   * @brief Constructs a pattern given a center point and shared
   * width, threshold and scale factor for all axes.
   * @note Type needs to be convertible to all Types...
   * @tparam Type Type for shared width, threshold and scale factor
   * @param c Pattern center point as tuple of Types...
   * @param w Pattern width shared by all axes
   * @param t Pattern convergence threshold shared by all axes
   * @param s Pattern shrinking scale factor shared by all axes
   */
  Pattern(const std::tuple<Types...> &c, const double w, const double t,
          const double s)
    : c {c}, w {Types {w}...}, t {Types {t}...}, s {Types {s}...}
  {}

  /**
   * @brief Initializes a new pattern search.
   * @tparam Callable Type of callable objective function
   * @param f The objective function (Callable)
   * @return New pattern search state
   */
  template <class Callable>
  auto init(Callable &&f) const
  {
    return State<decltype(f(Types {}...))> {
      c, w, std::apply(std::forward<Callable>(f), c)};
  }

  /**
   * @brief Evolves the pattern state by one step according to the
   * objective function.
   * @tparam Callable Type of callable objective function
   * @param f The objective function (Callable)
   * @param state The current pattern search state
   * @return New pattern search state
   */
  template <class Callable, class Type>
  auto step(Callable &&f, State<Type> state) const
  {
    // Flag indicating whether a shrink step needs to be done
    bool need_shrink = true;
    // Construct and iterate all pattern points centered at state
    for (auto &&x : state.make_points()) {
      // Skip already computed points
      if (x == state.argmin || x == state.argmin_prev) { continue; }
      // Evaluate function at new point
      auto y =
        std::apply(std::forward<decltype(f)>(f), std::forward<decltype(x)>(x));
      // Compare to current minimum
      if (y < state.min) {
        // Update tracked minimum
        state.min         = y;
        state.argmin_prev = state.argmin;
        state.argmin      = x;
        // The pattern moved - not need to shrink this step
        need_shrink = false;
      }
    }
    // Check whether shrink step is necessary
    if (need_shrink) {
      // Shrink the state (width) with scale factors s but not
      // below thresholds t
      shrink(state, std::index_sequence_for<Types...> {});
    }
    // Return modified state
    return state;
  }

  /**
   * @brief Tests whether the search using the state is done (converged).
   * @param state The current search state to test.
   * @return true is the search converged.
   */
  template <class Type>
  bool done(const State<Type> &state) const
  {
    return all_done(state, std::index_sequence_for<Types...> {});
  }

private:
  /**
   * @brief Tests whether all widths are below the corresponding
   * threshold.
   * @return true if all width below the threshold: w < t
   */
  template <std::size_t... Axes, class Type>
  bool all_done(const State<Type> &state, std::index_sequence<Axes...>) const
  {
    // Compare each width to corresponding threshold. Fold
    // comparisons linked with && over all Axes
    return ((std::get<Axes>(state.w) < std::get<Axes>(t)) && ...);
  }

  /**
   * @brief Shrinks specific axis of the pattern
   */
  template <std::size_t Axis, class Type>
  void shrink_axis(State<Type> &state) const
  {
    // Test whether the axis is still above the threshold
    if (std::get<Axis>(state.w) >= std::get<Axis>(t)) {
      // Axis is above threshold: Shrink by scale factor
      std::get<Axis>(state.w) *= (double) std::get<Axis>(s);
    }
  }

  /**
   * @brief Shrinks the pattern in place using the scale factors s
   *
   * @note Stops shrinking per dimension if the width reaches the
   * threshold: Only shrink if w >= t
   */
  template <std::size_t... Axes, class Type>
  void shrink(State<Type> &state, std::index_sequence<Axes...>) const
  {
    // Shrink along each axis via pack expansion over the comma
    // operator.
    (shrink_axis<Axes>(state), ...);
  }

  /**
   * @brief Shrinks the pattern in place using the scale factors s
   *
   * @note Stops shrinking per dimension if the width reaches the
   * threshold: Only shrink if w >= t
   */
  template <class Type>
  void shrink(State<Type> &state) const
  {
    // Shrink all axes of the index sequence over the Types...
    shrink(state, std::index_sequence_for<Types...> {});
  }
};
} // namespace search
#endif /*SEARCH_PATTERN_H*/
