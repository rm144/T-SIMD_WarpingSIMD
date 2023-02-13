// ===========================================================================
//
// SearchNelderMead.H --
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

#ifndef SEARCH_NELDER_MEAD_H
#define SEARCH_NELDER_MEAD_H

// std::abs
#include <cstdlib>
// Standard math functions
#include <cmath>
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
 * @brief N-dimensional Nelder-Mead search strategy
 * @tparam Types Types per dimension/axis of the search space
 */
template <class... Types>
class NelderMead
{
  /**
   * @brief Nelder-Mead search (simplex) state
   * @tparam Type Type of the output value of the objective function
   */
  template <class Type>
  class State
  {
    // Wrapping search strategy needs access
    friend NelderMead;

    // Number of dimensions/axes
    static constexpr auto N = sizeof...(Types);

    // Vertices of the simplex
    std::array<std::tuple<Types...>, N + 1> xs;
    // Objective function values at the vertices
    std::array<Type, N + 1> ys;
    // Counts update steps
    std::size_t steps = 0;

  public:
    /**
     * @brief Returns minimum value and the corresponding argument.
     * @return std::pair of {min,argmin}
     */
    std::pair<Type, std::tuple<Types...>> min_argmin() const
    {
      // Start searching for the minimum value at first vertex
      //  NOTE: This is necessary as the state might be in an
      //  unsorted state when requesting the minimum. Sorting
      //  here would be a more costly and non-const operation.
      std::size_t min_index = 0;
      // Iterate over all N + 1 vertices
      for (std::size_t i = 0; i <= N; ++i) {
        // Update the minimum index if a smaller value is found
        if (ys[i] < ys[min_index]) { min_index = i; }
      }
      // Return (min,argmin) at found index
      return std::pair {ys[min_index], xs[min_index]};
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
      // Iterate over all N + 1 vertices
      for (std::size_t i = 0; i <= N; ++i) {
        // Insert vertex coordinates + function value into
        // string representation
        repr << xs[i] << " " << ys[i] << " ";
      }
      // Get the string out of the string stream
      return repr.str();
    }

  private:
    /**
     * @brief Initializes a Nelder-Mead search state with simplex
     * vertices and corresponding function values.
     */
    State(std::array<std::tuple<Types...>, N + 1> xs,
          std::array<Type, N + 1> ys)
      : xs {xs}, ys {ys}
    {
      sort();
    }

    /**
     * @brief Sorts the state according to the current values of the
     * objective function.
     */
    void sort()
    {
      // Bubble Sort: At most N passes over the vertices are
      // necessary for sorting
      for (std::size_t n = 0; n <= N; ++n) {
        // Swapped at least one pair this pass?
        bool swapped = false;
        // Iterate all vertices
        for (std::size_t i = 1; i <= (N - n); ++i) {
          // Compare the pair of adjacent vertices
          if (ys[i - 1] > ys[i]) {
            // Out of order need to swap
            std::swap(ys[i - 1], ys[i]);
            std::swap(xs[i - 1], xs[i]);
            // Swapped, thus do not stop this pass
            swapped = true;
          }
        }
        // If no swap this pass sorting is done
        if (!swapped) { break; }
      }
    }

    /**
     * @brief Computes the centroid of the first N (not N+1)
     * vertices of the simplex.
     *
     * @return The centroid vertex x0 (mean of x1-xN)
     */
    std::tuple<Types...> centroid() const
    {
      // Arithmetic operator overloads operating on std::tuple
      using namespace tuple_operators;
      // Start with zero initialized tuple accumulator
      std::tuple<Types...> sum = {Types {0}...};
      // Iterate over the first N (not N+1) points
      for (std::size_t i = 0; i < N; ++i) {
        // Sum vertices
        sum += xs[i];
      }
      // Compute mean by dividing by number of vertices
      return sum / N;
    }

    /**
     * @brief Computes the sample standard deviation of the current
     * function values at the vertices.
     * @note This depends o the value Type to provide most of the
     * arithmetic operations or at least be implicitly convertible
     * to some arithmetic type.
     * @return Standard deviation of the function values
     */
    auto std() const
    {
      // Start accumulating the values
      decltype(Type {} + Type {}) sum = 0.0;
      // Iterate the function values
      for (std::size_t i = 0; i <= N; ++i) { sum += ys[i]; }
      // Divide by number of simplex points to get the mean
      decltype(Type {} + Type {}) mean = sum / (N + 1);
      // Start accumulating standard deviation
      sum = 0;
      // Iterate the function values
      for (std::size_t i = 0; i <= N; ++i) {
        sum += (ys[i] - mean) * (ys[i] - mean);
      }
      // Divide by number of simplex points to get the squared
      // standard deviation - square root to get the standard
      // deviation.
      return std::sqrt(sum / (N + 1));
    }

    /**
     * @brief Computes the axis-aligned bounding box (actually the
     * side lengths) of the simplex.
     */
    auto bounding_box() const
    {
      // Arithmetic operator overloads operating on std::tuple
      using namespace tuple_operators;
      // Start with minimum and maximum initialized to first
      // vertex
      auto min = xs[0], max = xs[0];
      // Iterate over all N + 1 vertices
      for (std::size_t i = 1; i <= N; ++i) {
        // Find elementwise (coordinate wise) minimum and
        // maximum
        min = tuple_operators::min(xs[i], min);
        max = tuple_operators::max(xs[i], max);
      }
      // Side lengths of the bounding box are the differences
      // between coordinate wise upper and lower bound
      return max - min;
    }
  };

  // Number of dimensions/axes
  static constexpr auto N = sizeof...(Types);

  // (Initial) Vertices of the simplex
  std::array<std::tuple<Types...>, sizeof...(Types) + 1> xs {};
  // Reflection, expansion, contraction, shrink coefficient and
  // convergence threshold
  double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5, t = 1.0;
  // Maximum number of steps before forcing termination
  std::size_t max_steps = 0;

public:
  /**
   * @brief Default initializes a Nelder-Mead search
   * @note This simplex does immediately converge (as all simplex points
   * are the same)
   */
  NelderMead() = default;

  /**
   * @brief Initializes the Nelder-mead search strategy simplex.
   * @param xs Initial vertices of the search simplex
   * @param alpha Reflection coefficient (sensible default is 1.0)
   * @param gamma Expansion coefficient (sensible default is 2.0)
   * @param ro Contraction coefficient (sensible default is 0.5)
   * @param sigma Shrink coefficient (sensible default is 0.5)
   * @param t Convergence threshold
   * @param max_steps Fail-safe condition: stop after exceeding max_steps
   */
  NelderMead(const std::array<std::tuple<Types...>, sizeof...(Types) + 1> &xs,
             const double alpha, const double gamma, const double rho,
             const double sigma, const double t, const std::size_t max_steps)
    : xs {xs}, alpha {alpha}, gamma {gamma}, rho {rho}, sigma {sigma}, t {t},
      max_steps {max_steps}
  {}

  /**
   * @brief Initializes a new Nelder-Mead search.
   * @tparam Callable Type of callable objective function
   * @param f The objective function (Callable)
   * @return New search state (current simplex)
   */
  template <class Callable>
  auto init(Callable &&f) const
  {
    // Prepare array of function values at initial vertices
    std::array<decltype(f(Types {}...)), N + 1> ys;
    // Iterate the initial simplex vertices
    for (std::size_t i = 0; i <= N; ++i) {
      // Evaluate the function at the vertex
      ys[i] = std::apply(std::forward<Callable>(f), xs[i]);
    }
    // Construct initial search state
    return State<decltype(f(Types {}...))> {xs, ys};
  }

  /**
   * @brief Evolves the search state (simplex) state by one step according
   * to the objective function.
   *
   * @tparam Callable Type of callable objective function
   * @param f The objective function (Callable)
   * @param state The current search state
   * @return New search state (current simplex)
   */
  template <class Callable, class Type>
  auto step(Callable &&f, State<Type> state) const
  {
    // Arithmetic operator overloads operating on std::tuple
    using namespace tuple_operators;
    // 0. Count the step
    state.steps += 1;
    // 1. Order according to the values at the vertices
    state.sort();
    // 2. Calculate the centroid of the simplex
    std::tuple<Types...> x0 = state.centroid();
    // 3. Compute the reflected point
    std::tuple<Types...> xr = x0 + alpha * (x0 - state.xs[N]);
    // Evaluate the function at the reflected point
    auto yr = std::apply(std::forward<decltype(f)>(f), xr);
    // Compare reflected point to second worst and best point
    if (state.ys[0] <= yr && yr < state.ys[N - 1]) {
      // Replace the worst with the reflected point
      state.xs[N] = xr;
      state.ys[N] = yr;
      // Stop here and return modified search state
      return state;
    }
    // 4. Compare the reflected point to the best point again
    if (yr <= state.ys[0]) {
      // Compute the expanded point
      std::tuple<Types...> xe = x0 + gamma * (xr - x0);
      // Evaluate the function at the expanded point
      auto ye = std::apply(std::forward<decltype(f)>(f), xe);
      // Compare the expanded to the reflected point
      if (ye < yr) {
        // Replace the worst point with the expanded point
        state.xs[N] = xe;
        state.ys[N] = ye;
        // Stop here and return modified search state
        return state;
      } else {
        // Replace the worst with the reflected point
        state.xs[N] = xr;
        state.ys[N] = yr;
        // Stop here and return modified search state
        return state;
      }
    }
    // 5. Compute the contracted point
    std::tuple<Types...> xc = x0 + rho * (state.xs[N] - x0);
    // Evaluate the function at the contracted point
    auto yc = std::apply(std::forward<decltype(f)>(f), xc);
    // Compare the contracted to the worst point
    if (yc < state.ys[N]) {
      // Replace the worst with the contracted point
      state.xs[N] = xc;
      state.ys[N] = yc;
      // Stop here and return modified search state
      return state;
    }
    // 6. Shrink the simplex by replacing all but the best point
    for (std::size_t i = 1; i < N; ++i) {
      // Shrink by adding and scaling difference to best
      auto xi = sigma * state.xs[i] + (1.0 - sigma) * state.xs[0];
      // Evaluate the function at the vertex
      auto yi = std::apply(std::forward<Callable>(f), xi);
      // Update state with new vertex and function value
      state.xs[i] = xi;
      state.ys[i] = yi;
    }
    // Return the modified state
    return state;
  }

  /**
   * @brief Tests whether the search using the state is done (converged).
   * @param state The current search state to test.
   * @return true if the search converged.
   */
  template <class Type>
  bool done(const State<Type> &state) const
  {
    // Test whether the maximum amount of steps reached or the
    // bounding box is smaller than the threshold.
    return state.steps > max_steps ||
           all_done(state.bounding_box(), std::index_sequence_for<Types...> {});
  }

private:
  /**
   * @brief Tests for convergence by comparing the side lengths of the
   * axis-aligned bounding box to the convergence threshold.
   * @tparam Axes Index sequence of axes to test
   * @param bounding_box Tuple of side lengths of the bounding box
   * @return true if all side lengths below the threshold: w < t
   */
  template <std::size_t... Axes>
  bool all_done(const std::tuple<Types...> bounding_box,
                std::index_sequence<Axes...>) const
  {
    // Compare each side length to threshold. Fold comparisons
    // linked with && over all Axes
    return ((std::get<Axes>(bounding_box) < Types {t}) && ...);
  }
};
} // namespace search
#endif /*SEARCH_NELDER_MEAD_H*/
