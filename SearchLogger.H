// ===========================================================================
// 
// SearchLogger.H --
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

#ifndef SEARCH_LOGGER_H
#define SEARCH_LOGGER_H

// std::tuple - Tuple of arbitrary values
#include <tuple>
// std::pair - Pair of two values
#include <utility>
// Generic list
#include <list>
// String streams (input/output from/to strings)
#include <sstream>

// Custom tuple operators
#include "TupleUtility.H"

/**
 * @brief Namespace containing search strategies
 */
namespace search {
/**
 * @brief Generic search strategy wrapper adding logging to the search
 */
template<class Searcher>
    struct Logger : public Searcher {
        // Hold reference to list writing the state log into
        std::list<std::string> &log;

        /**
         * @brief Constructs a logging variant of the Searcher
         * @param searcher The already constructed searcher to add the logging
         * @param log List to write the string representation of the search
         * state to.
         */
        Logger(const Searcher &searcher, std::list<std::string> &log)
            : Searcher{searcher}, log{log} {};

        /**
         * @brief Initializes a new search logging the initial state.
         * @tparam Callable Type of callable objective function
         * @param f The objective function (Callable)
         * @return New search state
         */
        template<class Callable>
            auto init(Callable &&f) const {
                // Forward arguments to wrapped init method implementation
                auto state = Searcher::init(f);
                // Add the initial state to the log and return it unmodified
                return log.push_back(state.str()), state;
            }

        /**
         * @brief Updates the searcher by one step logging the state
         * @tparam Callable Type of callable objective function
         * @tparam State Generic search state (must match the Searchers state)
         * @param f The objective function (Callable)
         * @param state The current pattern search state
         * @return New of the same type as the input search state
         */
        template<class Callable, class State>
            auto step(Callable &&f, State state) const {
                // Forward arguments to wrapped step method implementation
                state = Searcher::step(f, state);
                // Insert state string representation into logging list
                return log.push_back(state.str()), state;
            }
    };

// Grid search logging requires special treatment as there is no step method
// which can be hooked
template<class... Types>
    struct Logger<search::Grid<Types...>> : public search::Grid<Types...> {
        // Hold reference to list writing the state log into
        std::list<std::string> &log;

        /**
         * @brief Constructs a logging variant of the grid searcher
         * @param grid The already constructed searcher to add the logging
         * @param log List to write the string representation of the search
         * state to.
         */
        Logger(const search::Grid<Types...> &grid, std::list<std::string> &log)
            : search::Grid<Types...>{grid}, log{log} {};
    };

// Tag for marking the maybe_log function to add NO logging
struct NoLog {};

/**
 * @brief Maybe adds logging of the searcher step input states into the given
 * list
 * @tparam Searcher Searcher type to add logging to
 * @param searcher The searcher to log
 * @param log The list to write the search log to
 * @return Searcher augmented by logging to the given list
 */
template<class Searcher>
    auto maybe_log(Searcher &&searcher, std::list<std::string> &log) {
        return Logger{searcher, log};
    }

// Overload not adding logging
template<class Searcher>
    auto maybe_log(Searcher &&searcher, NoLog = NoLog{}) {
        return std::forward<Searcher>(searcher);
    }

// Exhaustive search minimizer function overload with logging
template<class Objective, class... Types>
    auto minimize(Objective &&f, Logger<Grid<Types...>> searcher) {
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
            // Build string representation by inserting into output string
            // stream
            std::ostringstream repr;
            // Build string representation of the objective function values
            tuple_operators::operator<<(repr, std::tuple{xs...}) << " " << y << " ";
            // Add objective function input+output to log
            searcher.log.push_back(repr.str());
        });
        // Return a pair of the found minimum (output) and the corresponding
        // argument (input) of f
        return std::pair<decltype(min), decltype(argmin)>{min, argmin};
    }
}

#endif /*SEARCH_LOGGER_H*/
