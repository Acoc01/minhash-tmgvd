#include "GraphCluster.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <set>

#include "utils/algorithms.hpp"
#include "Vertex.hpp"

namespace odsg {
//==============================================================================
GraphCluster::GraphCluster( const GraphCluster& cluster) {
    innerCluster = cluster.innerCluster;
    ptr_graph    =    cluster.ptr_graph;
}

GraphCluster
GraphCluster::operator = ( const GraphCluster& cluster ) {
    innerCluster = cluster.innerCluster;
    ptr_graph    =    cluster.ptr_graph;
    return *this;
}

//==============================================================================

void
GraphCluster::insert(Graph::const_iterator it) {
    assert(!algorithms::is_found(innerCluster, it));

    innerCluster.push_back(it);
}


void
GraphCluster::merge(const GraphCluster& cluster) {
    // TODO: assertion: *this and cluster doesn't have elements in common

    algorithms::vector_insert_back(innerCluster, cluster.innerCluster);
}


unsigned long
GraphCluster::arcsCount() const {       // Implementation from Graph::arcsCount()
    unsigned long arcs = 0;

    for (const_iterator it = begin(); it != end(); ++it) {
        arcs += it->second.size();      // Add the total of outlinks of each vertex
    }

    return arcs;
}


std::size_t
GraphCluster::nodesCount() const {      // Implementation from Graph::nodesCount()
    std::set<Vertex> nodes;

    for (const_iterator it = begin(); it != end(); ++it) {
        // The set nodes doesn't store duplicates, so we can to insert all without need to care
        nodes.insert(it->first);
        nodes.insert(it->second.begin(), it->second.end());
    }

    assert(nodes.size() >= listsCount());
    return nodes.size();
}


}   // namespace odsg
