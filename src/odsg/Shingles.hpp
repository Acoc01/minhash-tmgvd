#ifndef SRC_SHINGLES_HPP_INCLUDED
#define SRC_SHINGLES_HPP_INCLUDED

#include "Graph.hpp"

namespace odsg {


/*
 * Shingles are expected to be required only from the GraphPartitionerBySignature class.
 * The current implementation is most basic possible, supporting shingles, as defined in the literature, of size = 1
 * and a unique signature.
 */
class Shingles {
public:
    typedef unsigned int Signature;

    Shingles();

    Signature sign(const Graph::AdjacencyList&) const;

private:
    static const unsigned long bigPrime = 0x7FFFFFFF;   // Same as 2^31 - 1, that is prime

    unsigned int A;     // A, as in a X + b mod bigPrime
    unsigned int B;     // B, as in a X + b mod bigPrime
};


}       // namespace odsg
#endif  // SRC_SHINGLES_HPP_INCLUDED
