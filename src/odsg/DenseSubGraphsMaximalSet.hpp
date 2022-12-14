#ifndef SRC_DENSE_SUB_GRAPHS_MAXIMAL_SET_HPP_INCLUDED
#define SRC_DENSE_SUB_GRAPHS_MAXIMAL_SET_HPP_INCLUDED

#include <cstddef>      // std::size_t
#include <string>
#include <vector>

#include "DenseSubGraph.hpp"    // All of our container-like classes include the definition of the contained element

namespace odsg {


/*
 * A container-like class to store dense subgraphs 'maximal between themselves', i.e. no dense subgraph in the
 * collection is included (or 'superseded') in any other. Note that it implies no duplicates, what gives to the
 * class its name.
 */
class DenseSubGraphsMaximalSet {
public:

    // Types
    typedef std::vector<DenseSubGraph>::const_iterator const_iterator;

    // Constructors
    explicit DenseSubGraphsMaximalSet(bool asClique=false): dsgs(), onlyCentersMaximality(asClique) {}

    // Mutators
    bool insert(const DenseSubGraph&);

    template<typename ContainerT>
    bool insert(const ContainerT&);     // Very unlikely to be needed, except by unit tests

    // Iterators
    const_iterator begin() const { return dsgs.begin(); }
    const_iterator end() const { return dsgs.end(); }

    // Inspectors
    bool empty() const { return dsgs.empty(); }
    std::size_t size() const { return dsgs.size(); }

    void dump(const std::string& fileName, bool append, bool includeDescriptions) const;  // It can throw an exception

private:
    std::vector<DenseSubGraph> dsgs;
    bool onlyCentersMaximality;     // If true, it lets to users of this class to treat the dense subgraphs as cliques


    // Helper for insert()
    bool supersede(const DenseSubGraph&, const DenseSubGraph&) const;

    // General helpers
    typedef std::vector<DenseSubGraph>::iterator iterator;

    iterator begin() { return dsgs.begin(); }
    iterator end() { return dsgs.end(); }
};


template<typename ContainerT>
inline bool
DenseSubGraphsMaximalSet::insert(const ContainerT& somedsgs) {
    bool allInserted = true;
    for (typename ContainerT::const_iterator it = somedsgs.begin(); it != somedsgs.end(); ++it) {
        if (!insert(*it))
            allInserted = false;
    }
    return allInserted;     // true if all insertions were effective, false otherwise
}


}       // namespace odsg
#endif  // SRC_DENSE_SUB_GRAPHS_MAXIMAL_SET_HPP_INCLUDED
