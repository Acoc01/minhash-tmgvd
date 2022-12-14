#ifndef BIO_PRINTING_HPP_INCLUDED
#define BIO_PRINTING_HPP_INCLUDED

#include <ostream>
#include <functional>

#include "typedefs.hpp"

namespace bio_odsg {


namespace detail {      // Internal stuff to this header file

    struct TrivialTruePredicate: public std::unary_function<Complex, bool> {    // TODO: Get rid of this
        bool operator()(const Complex&) const { return true; }
    };
}


/*
 * Print a collection of complexes, one by line, using the protein names.
 * The given mapping MUST covers all the proteins present in the complexes.
 *
 * This first overloaded implementation requires a unary predicate to select the complexes to print;
 * The second implementation prints all of them.
 */
template<typename ContainerT, typename UnaryPredicate>
inline void
printComplexes(const ContainerT& complexes,
               const ProteinsMap& mapping,
               std::ostream& os,
               UnaryPredicate worthPrinting) {

    for (typename ContainerT::const_iterator it = complexes.begin(); it != complexes.end(); ++it) {
        const Complex& complex = *it;

        if (!worthPrinting(complex))
            continue;

        for (Complex::const_iterator protit = complex.begin(); protit != complex.end(); ++protit) {
            if (protit != complex.begin())
                os << " ";
            os << getProteinName(mapping, *protit);
        }
        os << '\n';
    }
}


template<typename ContainerT>
inline void
printComplexes(const ContainerT& complexes,
               const ProteinsMap& mapping,
               std::ostream& os) {

    printComplexes(complexes, mapping, os, detail::TrivialTruePredicate());
}


/*
 * Print a proteins mapping.
 * It's printed sorting by protein name.
 */
inline void
printMapping(const ProteinsMap& mapping, std::ostream& os) {
    for (ProteinsMap::const_iterator it = mapping.begin(); it != mapping.end(); ++it) {
        os << it->first << ' ' << it->second << '\n';
    }
}


}       // namespace bio_odsg
#endif  // BIO_PRINTING_HPP_INCLUDED
