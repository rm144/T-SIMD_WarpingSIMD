#ifndef SEARCH_H
#define SEARCH_H

// std::tuple - Tuple of arbitrary values
#include <tuple>
// std::pair - Pair of two values
#include <utility>
// Infinity
#include <limits>

/**
 * @brief Namespace containing search strategies
 */
namespace search {
/**
 * @brief Minimizes the objective function using a generic, search strategy.
 * @tparam Objective Type of the objective function to minimize
 * @tparam Method Type of the optimization method to use
 * @tparam Types Argument types of the objective function (Objective and Method
 * must agree on these)
 * @param f Objective function object to minimize
 * @param searcher Search method instance to use (contains the parameters)
 * @return std::pair of the found min and argmin
 */
template<class Objective, template<class...> class Method, class... Types>
    auto minimize(Objective &&f, Method<Types...> searcher) {
        // Initialize the searcher state
        auto state = searcher.init(f);
        // Iterate until converged to some solution
        while (!searcher.done(state)) {
            // Do one step of the search method
            state = searcher.step(f, state);
        }
        // Return a pair of the found minimum (output) and the corresponding
        // argument (input) of f
        return state.min_argmin();
    }

// Maximizes the objective function f(xs...) by minimizing -f(xs...)
template<class Objective, template<class...> class Method, class... Types>
    auto maximize(Objective &&f, Method<Types...> searcher) {
        // Construct the negated objective
        auto F = [f](Types... xs) {
            return -f(xs...);
        };
        // Minimize the negated objective
        auto[min, argmin] = minimize(F, searcher);
        // Return negated minimum and the argument
        return std::pair{-min, argmin};
    }

/**
 * @brief Tag to select unconstrained search with similar signature as the
 * constrained search
 */
class Unconstrained {};

/**
 * @brief Minimizes the objective function using a generic, search strategy.
 * @note Explicit tag for selecting unconstrained search
 * @tparam Objective Type of the objective function to minimize
 * @tparam Method Type of the optimization method to use
 * @tparam Types Argument types of the objective function (Objective and Method
 * must agree on these)
 * @param f Objective function object to minimize
 * @param searcher Search method instance to use (contains the parameters)
 * @return std::pair of the found min and argmin
 */
template<class Objective, template<class...> class Method, class... Types>
    auto minimize(Objective &&f, Method<Types...> searcher, Unconstrained) {
        // Just forward arguments to unconstrained minimizer
        return minimize(f, searcher);
    }

// Explicitly unconstrained maximizes the objective function
template<class Objective, template<class...> class Method, class... Types>
    auto maximize(Objective &&f, Method<Types...> searcher, Unconstrained) {
        // Just forward arguments to unconstrained maximizer
        return maximize(f, searcher);
    }

// Minimizes the objective f subject to the (boolean) constraint g
template<class Objective, class Method, class Constraint>
    auto minimize(Objective &&f, Method searcher, Constraint g) {
        // Create new function converting the boolean constraint into not
        // calling f + returning infinity
        auto F = [f, g](auto... xs) {
            // Test the constraint before evaluating f
            if (g(xs...)) {
                // Return f evaluated at xs...
                return f(xs...);
            }
            // Return infinity of appropriate type
            return std::numeric_limits<decltype(f(xs...))>::infinity();
        };
        // Optimize the wrapper objective
        return minimize(F, searcher);
    }

// Maximizes the objective f subject to the boolean constraint g
template<class Objective, class Method, class Constraint>
    auto maximize(Objective &&f, Method searcher, Constraint g) {
        // Construct the negated objective
        auto F = [f](auto &&... xs) {
            return -f(std::forward<decltype(xs)>(xs)...);
        };
        // Minimize the negated objective
        auto[min, argmin] = minimize(F, searcher, g);
        // Return negated minimum and the argument
        return std::pair{-min, argmin};
    }
}

#endif /*SEARCH_H*/
