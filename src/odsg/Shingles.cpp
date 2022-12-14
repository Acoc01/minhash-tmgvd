#include "Shingles.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <cstddef>      // NULL, std::size_t
#include <string>
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time

#include <boost/functional/hash.hpp>

#include "utils/algorithms.hpp"
#include "utils/strings.hpp"

namespace odsg {


Shingles::Shingles() {
    std::srand(std::time(NULL));
    A = (std::rand() % bigPrime) + 1;
    B = (std::rand() % bigPrime) + 1;
}


Shingles::Signature
Shingles::sign(const Graph::AdjacencyList& outlinks) const {
    assert(!outlinks.empty());
    assert(algorithms::has_unique(outlinks));

    boost::hash<std::string> stringHash;
    Signature minShingleHash = bigPrime;

    for (Graph::AdjacencyList::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
        std::size_t shingleID = stringHash(strings::to_str(*vxit));

        // TODO: I need to reexamine these calculations; these aren't completely OK
        Signature shingleHash = (((unsigned long) A * (unsigned long) shingleID) + B) % bigPrime;
        if (minShingleHash > shingleHash) {
            minShingleHash = shingleHash;
        }
    }
    return minShingleHash;
}


}   // namespace odsg
