#ifndef BIO_METRICS_HPP_INCLUDED
#define BIO_METRICS_HPP_INCLUDED

#include <algorithm>
#include <functional>   // std::unary_function

#include <odsg/utils/algorithms.hpp>

#include "typedefs.hpp"

namespace bio_odsg {


/*
 * Returns the overlap score (OS) between two complexes.
 */
inline double
overlapScore(const Complex& complex, const Complex& otherComplex) {
    unsigned int overlap = odsg::algorithms::set_intersection_count(complex, otherComplex);

    return double(overlap * overlap) / (complex.size() * otherComplex.size());
}


/*
 * Return if a complex is matched by at least one element in the given collection.
 *
 * Two complexes are considered as matched if the overlap score between these two is equal or greater than a given
 * threshold.
 */
template<typename ContainerT>
inline bool
matched(const Complex& matchingComplex,
        const ContainerT& complexes,
        float minOverlapScore=0.2) {

    for (typename ContainerT::const_iterator it = complexes.begin(); it != complexes.end(); ++it) {
        if (overlapScore(matchingComplex, *it) >= minOverlapScore) {
            return true;
        }
    }
    return false;
}


/*
 * Unary predicate corresponding to the 'matched' function
 */
template<typename ContainerT>
class MatchedPredicate : public std::unary_function<Complex, bool> {
public:
    MatchedPredicate(const ContainerT& cc, double os): complexes(cc), overlapScore(os) {}

    bool operator()(const Complex& matchingComplex) const {
        return matched(matchingComplex, complexes, overlapScore);
    }

private:
    const ContainerT complexes;
    const double overlapScore;
};


/*
 * Count the number of complexes in the first collection matched by at least one element in the second collection.
 */
template<typename ContainerT>
inline unsigned int
countMatchedComplexes(const ContainerT& matchingComplexes,
                      const ContainerT& complexes,
                      float minOverlapScore=0.2) {

    return std::count_if(matchingComplexes.begin(),
                         matchingComplexes.end(),
                         MatchedPredicate<ContainerT>(complexes, minOverlapScore));
}


}       // namespace bio_odsg
#endif  // BIO_METRICS_HPP_INCLUDED
