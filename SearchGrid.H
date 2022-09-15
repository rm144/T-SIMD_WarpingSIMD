#ifndef SEARCH_GRID_H
#define SEARCH_GRID_H

// std::tuple - Tuple of arbitrary values
#include <tuple>
// std::pair - Pair of two values
#include <utility>

// Namespace containing search strategies
namespace search {
/**
 * @brief 1-dimensional discrete range
 */
template<class Type>
    struct Range {
        /**
         * @brief Lower bound of the range
         */
        Type min;

        /**
         * @brief Upper bound of the range
         */
        Type max;

        /**
         * @brief Step size for iterating the range
         */
        Type step;

        /**
         * @brief Iterates the range [min;max] in step size increments
         *  calling f
         */
        template<class Callable>
            constexpr inline void for_each(Callable &&f) const {
                // Iterate the range using classic for loop
                for (Type x = min; x <= max; x += step) {
                    // Call the function at current x
                    f(x);
                }
            }
    };

// Template argument deduction guide
template<class Type>
    Range(Type, Type, Type) -> Range<Type>;

/**
 * @brief N-dimensional discrete grid search domain
 */
template<class... Types>
    class Grid;

// Template argument deduction guide for constructing a Grid from Ranges
template<class... Types>
    Grid(const Range<Types> &...) -> Grid<Types...>;

// General case of the N-dimensional Grid: Recursive Grid chain
template<class Type, class... Types>
    class Grid<Type, Types...> {
        // Wrapped 1d range of the same type (this level)
        Range<Type> range;

        // (N-1)d sub-grids (recursive Grid chain)
        Grid<Types...> ranges;

    public:
        /**
         * @brief Constructs Grid<Type,Types...> from 1d Range<Type> and
         *  N-1 1d ranges Range<Types>
         */
        explicit Grid(const Range<Type> &r, const Range<Types> &... rs)
            : range{r}, ranges{rs...} {}

        /**
         * @brief Iterates the grid [min;max]^N in step size increments
         *  calling f
         */
        template<class Callable>
            constexpr inline void for_each(Callable &&f) const {
                // Iterate the wrapped (outermost) range
                range.template for_each([&](Type x) {
                    // Recursive iteration of the sub-grid
                    ranges.for_each([&](Types... xs) {
                        // Call the function at x, xs... (x0,x1,...,xN)
                        f(x, xs...);
                    });
                });
            }
    };

// Special case of the N-dimensional Grid: 1-dimensional grid
template<class Type>
    class Grid<Type> {
        // Wrapped 1d range of the same type
        Range<Type> range;

    public:
        /**
         * @brief Constructs 1d Grid<Type> from 1d Range<Type>
         */
        explicit Grid(const Range<Type> &r) : range{r} {}

        /**
         * @brief Iterates the grid [min;max] in step size increments
         *  calling f
         */
        template<class Callable>
            constexpr inline void for_each(Callable &&f) const {
                // Iterate the wrapped range
                range.template for_each(f);
            }
    };

// Exhaustive search minimizer function overload
template<class Objective, class... Types>
    auto minimize(Objective &&f, Grid<Types...> searcher) {
        // Setup (uninitialized) minimum to track
        //  NOTE: Initialization depends on default constructor
        decltype(f(Types{}...)) min;
        // Initialize argument (index) of current minimum
        //  NOTE: Initialization depends on default constructor
        std::tuple<Types...> argmin;
        // Tracked (min,argmin) is not initialized
        bool uninitialized = true;
        // Iterate all candidate solutions of the search domain
        searcher.for_each([&](auto &&... xs) {
            // Evaluate the objective function on the candidate solution
            auto y = f(xs...);
            // Compare to current minimum
            if (uninitialized || y < min) {
                // Update tracked minimum and its argument
                min = y; argmin = {xs...};
                // Tracked (min,argmin) is now properly initialized
                uninitialized = false;
            }
        });
        // Return a pair of the found minimum (output) and the corresponding
        // argument (input) of f
        return std::pair<decltype(min), decltype(argmin)>{min, argmin};
    }
}

#endif /*SEARCH_GRID_H*/
