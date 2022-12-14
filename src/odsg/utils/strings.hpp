#ifndef SRC_UTILS_STRINGS_HPP_INCLUDED
#define SRC_UTILS_STRINGS_HPP_INCLUDED

#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <typeinfo>     // typeid; it's not portable, but that isn't a problem here
#include <ios>          // std::boolalpha

namespace odsg {
namespace strings {


/*
 * Convert any output-streamable thing to a string.
 * From:
 *   https://isocpp.org/wiki/faq/misc-technical-issues#convert-string-to-any
 */
template<typename T>
inline std::string
to_str(const T& value) {

    std::ostringstream oss;
    if (!(oss << value))
        throw std::runtime_error(std::string("Bad convertion to string from ") + typeid(value).name());
    return oss.str();
}

template<>
inline std::string      // This specialization converts boolean values to "true" or "false"
to_str<bool>(const bool& value) {

    std::ostringstream oss;
    if (!(oss << std::boolalpha << value))
        throw std::runtime_error(std::string("Bad convertion to string from ") + typeid(value).name());
    return oss.str();
}


/*
 * Convert a string to any input-streamable thing.
 * From:
 *   https://isocpp.org/wiki/faq/misc-technical-issues#convert-string-to-any
 */
template<typename T>
inline T
str_to(const std::string& str) {

    T value;
    std::istringstream iss(str);
    if (!(iss >> value)) {
        throw std::runtime_error(std::string("Bad convertion to ") + typeid(value).name());
    }
    return value;
}

template<>
inline std::string
str_to(const std::string& str) {
    return str;
}


/*
 * Split a string formed by a collection of values separated (only) by some delimiter, and returning then the values.
 * Consecutive delimiters are treated as only one. For example:
 *   split<int>(",1, 3, 6,,, 11", ',') ---> {1, 3, 6, 11}
 *
 * It can throw an exception in many cases, e.g. having only whitespaces between comma-delimiters.
 */
template<typename T>
inline std::vector<T>
split(const std::string& str, char delimiter=',') {
    std::istringstream iss(str);
    std::string token;
    std::vector<T> values;

    while (std::getline(iss, token, delimiter)) {
        if (!token.empty())
            values.push_back(str_to<T>(token));
    }
    return values;
}


}       // namespace strings
}       // namespace odsg
#endif  // SRC_UTILS_STRINGS_HPP_INCLUDED
