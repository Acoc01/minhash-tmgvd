#include "GraphPartitioner.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro

#include "utils/algorithms.hpp"
#include "Graph.hpp"
#include "GraphCluster.hpp"

namespace odsg {


GraphPartitioner::GraphPartitioner(const Graph* g): graph(g) {
    assert(g);
}


GraphCluster
GraphPartitioner::next(unsigned long minClusterArcsCount) {
    assert(graph->isMineable());

    // When the cluster built by getNext() is not enough big, ask another and merge both of them. Keep repeating it
    // until to fulfill the desired size or getNext() returns an empty cluster, the indication that the partitioning
    // has reached the end.
    // Note that it implicitly manages subsequent calls to this method after ending too.
    
//==============================================================================    
    GraphCluster mergedCluster(graph), cluster;
//==============================================================================

    do {
        cluster = getNext();
        mergedCluster.merge(cluster);
    } while (!cluster.empty() && mergedCluster.arcsCount() < minClusterArcsCount);

    return mergedCluster;
}


//// GraphPartitionerByInitialOutlink /////////////////////////////////////////////////////////////////////////////////

GraphPartitionerByInitialOutlink::GraphPartitionerByInitialOutlink(const Graph* g)
: GraphPartitioner(g), initialOutlinks() {

    assert(graph->isMineable());
    // The case of graph being empty is managed too, implicitly

    for (Graph::const_iterator it = graph->begin(); it != graph->end(); ++it) {
        const Graph::AdjacencyList& outlinks = it->second;
        assert(!outlinks.empty());

        initialOutlinks[outlinks[0]].insert(it);
    }
    assert(initialOutlinks.size() <= graph->listsCount());

    nextInitialOutlink = initialOutlinks.begin();
}


GraphCluster
GraphPartitionerByInitialOutlink::getNext() {
    if (nextInitialOutlink == initialOutlinks.end())
        return GraphCluster();      // An empty cluster indicates to the caller to have reached the end

    std::map<Vertex, GraphCluster>::const_iterator currentInitialOutlink = nextInitialOutlink;
    ++nextInitialOutlink;

    assert(!currentInitialOutlink->second.empty());
    return currentInitialOutlink->second;
}


//// GraphPartitionerBySignature //////////////////////////////////////////////////////////////////////////////////////

GraphPartitionerBySignature::GraphPartitionerBySignature(const Graph* g)
: GraphPartitioner(g), signatures() {

    assert(graph->isMineable());
    // The case of graph being empty is managed too, implicitly

    Shingles shingle;
    for (Graph::const_iterator it = g->begin(); it != g->end(); ++it) {
        const Graph::AdjacencyList& outlinks = it->second;
        assert(!outlinks.empty());

        Shingles::Signature signature = shingle.sign(outlinks);
        signatures[signature].insert(it);
    }
    assert(signatures.size() <= graph->listsCount());

    nextSignature = signatures.begin();
}


GraphCluster
GraphPartitionerBySignature::getNext() {
    if (nextSignature == signatures.end())
        return GraphCluster();      // An empty cluster indicates to the caller to have reached the end

    std::map<Shingles::Signature, GraphCluster>::const_iterator currentSignature = nextSignature;
    ++nextSignature;

    assert(!currentSignature->second.empty());
    return currentSignature->second;
}


}   // namespace odsg
