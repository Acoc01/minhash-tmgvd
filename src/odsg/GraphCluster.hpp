#ifndef SRC_GRAPH_CLUSTER_HPP_INCLUDED
#define SRC_GRAPH_CLUSTER_HPP_INCLUDED

#include <cstddef>      // std::size_t
#include <vector>

#include <boost/iterator/indirect_iterator.hpp>

#include "Graph.hpp"

namespace odsg {


/*
 * A cluster defines a non-mutable subset of the adjacency lists of a Graph object. It's useful to define an
 * implicit subgraph with a low overhead.
 * Cluster objects are expected to be created only from GraphPartitioner derived classes.
 *
 * For a mineable graph, the implicit graph defined by a cluster is mineable too; this property is recognized and
 * exploited in many places.
 *
 * The use of boost::indirect_iterator lets us to iterate over GraphCluster objects in the same way as iterating
 * over Graph objects:
 *   http://www.boost.org/doc/libs/1_58_0/libs/iterator/doc/indirect_iterator.html
 */
class GraphCluster {
public:

    // Types
    typedef boost::indirect_iterator<std::vector<Graph::const_iterator>::const_iterator>::value_type value_type;
    typedef boost::indirect_iterator<std::vector<Graph::const_iterator>::const_iterator> const_iterator;

    // Constructor
    GraphCluster(): innerCluster(), ptr_graph(NULL){}

//==============================================================================    
    GraphCluster(const Graph* graph): innerCluster(), ptr_graph(graph) {}
    
    GraphCluster( const GraphCluster& );
//==============================================================================

    // Iterators
    const_iterator begin() const { return boost::make_indirect_iterator(innerCluster.begin()); }
    const_iterator end() const { return boost::make_indirect_iterator(innerCluster.end()); }

    // Mutators
    void insert(Graph::const_iterator);
    void merge(const GraphCluster&);

    void clear() { innerCluster.clear(); }

    // Inspectors
    bool empty() const { return innerCluster.empty(); }

    std::size_t listsCount() const { return innerCluster.size(); }
    std::size_t nodesCount() const;
    unsigned long arcsCount() const;

    const Graph* get_ptrGraph() const { return ptr_graph; }
    
    GraphCluster operator = ( const GraphCluster& cluster );

private:

    std::vector<Graph::const_iterator> innerCluster;

//==============================================================================   
    const Graph* ptr_graph;
//==============================================================================    
};


}       // namespace odsg
#endif  // SRC_GRAPH_CLUSTER_HPP_INCLUDED
