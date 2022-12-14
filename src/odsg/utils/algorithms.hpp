#ifndef SRC_UTILS_ALGORITHMS_HPP_INCLUDED
#define SRC_UTILS_ALGORITHMS_HPP_INCLUDED

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>     // std::inserter
#include <functional>   // std::greater

namespace odsg {


/*
 * A collection of not-always-very-generic utilities expanding the capacities of the <algorithm> library, sort of.
 */
namespace algorithms {


/*
 * Typing shortcut for std::find. The parameters order follows the same logic.
 */
template<typename ContainerT, typename T>
inline bool
is_found(const ContainerT& container,
         const T& value) {

    return std::find(container.begin(), container.end(), value) != container.end();
}


/*
 * Generalization of is_found for more of a value to search, in unordened containers.
 */
template<typename ContainerT, typename ContainedT>
inline bool
are_found(const ContainerT& container,
          const ContainedT& values) {

    for (typename ContainedT::const_iterator it = values.begin(); it != values.end(); ++it) {
        if (!is_found(container, *it))
            return false;
    }
    return true;
}


/*
 * An alternative to std::find but using binary search over a sorted range.
 * If only checking membership is required, use std::binary_search instead.
 */
template<typename IteratorT, typename T>
inline IteratorT
binary_find(IteratorT begin,
            IteratorT end,
            const T& value) {

    IteratorT it = std::lower_bound(begin, end, value);
    return (it != end && !(value < *it)) ? it : end;
}


/*
 * Typing shortcut for std::includes. The parameters order follows the same logic.
 */
template<typename T>
inline bool
set_includes(const std::set<T>& superset,
             const std::set<T>& subset) {

    return std::includes(superset.begin(), superset.end(),
                         subset.begin(), subset.end());
}


/*
 * Typing shortcut for std::set_intersection. The parameters order follows the same logic.
 */
template<typename T>
inline void
set_intersection(const std::set<T>& lset,
                 const std::set<T>& rset,
                 std::set<T>& result) {

    std::set_intersection(lset.begin(), lset.end(),
                          rset.begin(), rset.end(),
                          std::inserter(result, result.end()));
}


/*
 * Return the number of elements in common between two sets.
 *
 * The implementation is strongly based in the 'possible implementation' for std::set_intersection from:
 *   http://en.cppreference.com/w/cpp/algorithm/set_intersection
 */
template<typename T>
inline typename std::set<T>::size_type
set_intersection_count(const std::set<T>& lset,
                       const std::set<T>& rset) {

    typename std::set<T>::const_iterator lit = lset.begin();
    typename std::set<T>::const_iterator rit = rset.begin();
    typename std::set<T>::size_type count = 0;

    while (lit != lset.end() && rit != rset.end()) {
        if (*lit < *rit) {
            ++lit;
        } else {
            if (!(*rit < *lit)) {
                ++count;
                ++lit;
            }
            ++rit;
        }
    }
    return count;
}


/*
 * Check if an unordered container (as std::vector) has all unique elements.
 */
template <typename ContainerT>
inline bool
has_unique(ContainerT container) {      // Passed by copy

    std::sort(container.begin(), container.end());
    return std::adjacent_find(container.begin(), container.end()) == container.end();
}


/*
 * Check if a container is already sorted in crecient order.
 * It invokes operator> for the elements.
 */
template <typename ContainerT>
inline bool
is_sorted(const ContainerT& container) {

    return std::adjacent_find(container.begin(), container.end(),
                              std::greater<typename ContainerT::value_type>()) == container.end();
}


/*
 * Sort of typing shorcut for std::vector<T>::insert, for a vector at the end of another.
 * Same preconditions applies, as not to try inserting a vector with itself.
 */
template<typename T>
inline void
vector_insert_back(std::vector<T>& vec, const std::vector<T>& inserting) {
    vec.reserve(vec.size() + inserting.size());
    vec.insert(vec.end(), inserting.begin(), inserting.end());
}


}       // namespace algorithms
}       // namespace odsg
#endif  // SRC_UTILS_ALGORITHMS_HPP_INCLUDED
