// ===========================================================================
// 
// Env.H --
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

#ifndef EXP_ENV_H
#define EXP_ENV_H

// C standard library header for C++ (contains std::getenv)
#include <cstdlib>
// C++ string stream input/output
#include <sstream>
// C++ generic tuple
#include <tuple>
// C++ string
#include <string>
// Generic list of objects
#include <list>
// Generic vector of elements
#include <vector>

// Environment variable based configuration options
namespace env_config {
/**
 * @brief List tracking all (name,value,default) option triples
 * @note env_config::get inserts into this list
 */
inline std::list<std::tuple<std::string, std::string, std::string>> options;

/**
 * @brief Prints all collected environment options to output stream
 */
void dump(std::ostream &stream, const bool dump_default = true) {
    // Iterate list of collected (name, value, default) triples
    for (const auto &[name, value, default_value]: options) {
        // Test whether this value is still default and can be dropped
        if(!dump_default && value == default_value) {
            continue;
        }
        // Insert (print) triple to stream
        stream << "# " << name << " = " << value
               << " (" << default_value << ") " << std::endl;
    }
}

/**
 * @brief Input stream extraction of generic std::vector of Type
 * @tparam Type Element type of the vector
 * @param lhs Reference to the input stream to extract values from
 * @param rhs Reference to the vector to insert values into
 * @return Reference to the modified input stream
 */
template<class Type>
    std::istream &operator>>(std::istream &&lhs, std::vector<Type> &rhs) {
        // Temporary value to read into
        Type value;
        // Read values while stream is ok (!fail())
        while ((lhs >> value).ignore()) {
            // Insert value into the vector
            rhs.push_back(value);
        }
        // Return reference to modified stream
        return lhs;
    }

/**
 * @brief Output stream insertion of generic std::vector of Type
 * @tparam Type Element type of the vector
 * @param lhs Reference to the output stream to insert values into
 * @param rhs Reference to the vector to get the values from
 * @return Reference to the modified output stream
 */
template<class Type>
    std::ostream &operator<<(std::ostream &lhs, const std::vector<Type> &rhs) {
        // Iterate over all elements of the vector
        for (const Type &value: rhs) {
            // Insert element into output stream
            lhs << value << ' ';
        }
        // Return reference to modified stream
        return lhs;
    }

// Generic value to string conversion
template<class Type>
    std::string str(Type &&value) {
        return std::to_string(value);
    }

// Value to string conversion: Special case - value is a string
std::string str(std::string value) {
    return value;
}

// Generic value to string conversion of std::vector
template<class Type>
    std::string str(const std::vector<Type> value) {
        std::ostringstream stream;
        stream << value;
        return stream.str();
    }

// Generic value to string conversion of Angle<Unit,Type>
template<class Unit, class Type>
    std::string str(const Angle<Unit, Type> value) {
        std::ostringstream stream;
        stream << value;
        return stream.str();
    }

/**
 * @brief Gets environment variable config of Type if set or default value
 * @tparam Type Type of the environment config to read
 * @param name Name of the environment variabel to read
 * @param default_value Default value if the variabel is not set
 * @return Value of the config option
 */
template<class Type>
    Type get(const std::string &name, const Type default_value) {
        // Get pointer to environment variable string (if it exists)
        const char *env_var = std::getenv(name.c_str());
        // Value of the environment variable initialized to default in case it
        // is not set
        Type value;
        // If the environment variable is set (null if not set)
        if (env_var) {
            // Extract formatted value of Type from environment variable
            std::istringstream(env_var) >> value;
        } else {
            // Set default option value
            value = default_value;
        }
        // Track the option triple in global list
        options.push_back({name, str(value), str(default_value)});
        // Return the extracted value
        return value;
    }
}

// 20. Sep 22 (Jonas Keller): prefix VEC, BOOL, INT, DOUBLE and DOUBLE_VEC
// with ENV_ to avoid name conflicts with other libraries
// (e.g. windows.h defines BOOL too)

/**
 * @brief Creates a named Type environment variable configuration option
 */
#define ENV(Type, Name, Default) \
    const auto Name = env_config::get<Type>(#Name, (Default))

/**
 * @brief Creates a named Type environment variable configuration option of a
 * type extraction method env_config::get is derived from a base type.
 */
#define DERIVE_ENV(Type, From, Name, Default) \
    const Type Name = (Type)env_config::get<From>(#Name, (From)(Default));

/**
 * @brief Creates a named Type environment variable vector configuration option
 */
#define ENV_VEC(Type, Name, ...)     \
    const std::vector<Type> Name \
        = env_config::get<std::vector<Type>>(#Name, {__VA_ARGS__})

/**
 * @brief Creates a named bool environment variable configuration option
 */
#define ENV_BOOL(Name, Default) ENV(bool, Name, Default)

/**
 * @brief Creates a named int environment variable configuration option
 */
#define ENV_INT(Name, Default) ENV(int, Name, Default)

/**
 * @brief Creates a named double environment variable configuration option
 */
#define ENV_DOUBLE(Name, Default) ENV(double, Name, Default)

/**
 * @brief Creates a named vector of double environment variable configuration
 * option
 */
#define ENV_DOUBLE_VEC(Name, ...) ENV_VEC(double, Name, __VA_ARGS__)

#endif /*EXP_ENV_H*/
