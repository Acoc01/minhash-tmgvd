#ifndef SRC_GRAPH_PARTITIONER_HPP_INCLUDED
#define SRC_GRAPH_PARTITIONER_HPP_INCLUDED

#include <map>

#include "Vertex.hpp"
#include "Shingles.hpp"

namespace odsg {


class Graph;
class GraphCluster;

/*
 * GraphPartitioner objects define a partition over the adjacency lists of a graph: a collection of clusters with no
 * overlap between them and no leaving single adjacency lists. The lifetime of the original graph must extend past
 * any use of the respective GraphPartitioner object, as a pointer to that graph is kept.
 *
 * GraphPartitioner is an abstract class, and each derived class define a different strategy to partition graphs,
 * according to different ideas leading to build such collection of dags from where 'better' dense subgraphs can be
 * mined.
 */
class GraphPartitioner {
public:
    GraphPartitioner(const Graph*);
    virtual ~GraphPartitioner() {};

    GraphCluster next(unsigned long minClusterArcsCount);

protected:
    const Graph* const graph;

private:
    virtual GraphCluster getNext() = 0;
};


//// GraphPartitionerByInitialOutlink /////////////////////////////////////////////////////////////////////////////////

class GraphPartitionerByInitialOutlink : public GraphPartitioner {
public:
    GraphPartitionerByInitialOutlink(const Graph*);

private:
    std::map<Vertex, GraphCluster> initialOutlinks;
    std::map<Vertex, GraphCluster>::const_iterator nextInitialOutlink;

    /*virtual*/ GraphCluster getNext();
};


//// GraphPartitionerBySignature //////////////////////////////////////////////////////////////////////////////////////

class GraphPartitionerBySignature : public GraphPartitioner {
public:
    GraphPartitionerBySignature(const Graph*);

private:
    std::map<Shingles::Signature, GraphCluster> signatures;
    std::map<Shingles::Signature, GraphCluster>::const_iterator nextSignature;

    /*virtual*/ GraphCluster getNext();
};


}       // namespace odsg
#endif  // SRC_GRAPH_PARTITIONER_HPP_INCLUDED
