#ifndef SRC_GRAPH_HPP_INCLUDED
#define SRC_GRAPH_HPP_INCLUDED

#include <cstddef>      // std::size_t
#include <string>
#include <vector>
#include <map>
#include <iosfwd>
#include <functional>   // std::binary_function
#include <algorithm>    // std::sort, ...

#include      "Vertex.hpp"
#include "WGraphTypes.hpp"
#include    "WedgeMap.hpp"

namespace odsg {


class GraphCluster;

/*
 * Graph objects are managed as a map of vertexes and their respective adjacency lists (aka outlinks).
 *
 * For all the constructors is assumed that the input data is 'consistent' (or 'without duplication'), e.g.
 * no multiple adjacency lists for a same vertex, or no duplicated vertexes inside the same adjacency list.
 * It's responsibility of the caller ensuring it, otherwise it will trigger undefined behaviour.
 */
class Graph {
public:
    //// Types ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector<Vertex> AdjacencyList;
    typedef std::map<Vertex, AdjacencyList>::const_iterator const_iterator;


    //// Constructors /////////////////////////////////////////////////////////////////////////////////////////////

    /*
     * Build an empty graph; it's mineable by definition.
     */
    Graph();
    virtual ~Graph() {}

    /*
     * Built a graph from a input text file.
     * Check the contents of the graphs/ folder in the root of the project for examples of the file formatting.
     *
     * It's possible to specify if the adjacency lists in the input data came already sorted by id. It's
     * responsibility of the caller ensuring it; it's exposed as a performance boost, as it lets to skip the slow
     * sorting phase in Graph::rebuildForMining(). In case of doubt, don't set it, otherwise it will trigger
     * undefined behaviour.
     */
    explicit Graph(const std::string& fileName, bool comeSortedByVertex=false);     // It can throw an exception

    /*
     * Built a graph from a cluster.
     * It's assumed that the cluster is from a mineable graph.
     */
    explicit Graph(const GraphCluster&);

    /*
     * Built a graph from a raw map from vertexes to adjacency lists. Provide compatibility with other
     * implementations.
     */
    template<typename ContainerT>
    explicit Graph(const std::map<Vertex, ContainerT>&);


    //// Iterators ////////////////////////////////////////////////////////////////////////////////////////////////

    const_iterator begin() const { return innerGraph.begin(); }
    const_iterator end() const { return innerGraph.end(); }


    //// Mutators /////////////////////////////////////////////////////////////////////////////////////////////////

    /*
     * Do the graph 'mineable', i.e. rebuild it to be suitable to construct a dag from it.
     * Any graph can be rebuilt to be mineable.
     *
     * The most important property of a mineable graph is that all its adjacency lists are sorted according to a
     * strict total order, thus ensuring that a dag built from it will not have cycles.
     *
     * The form without parameters sorts the adjacency lists only when it's required (e.g. it's skipped if the graph
     * was built from a file with its data already sorted) and using Graph::VertexFrequencyComparer.
     * For the second overloaded form, the sorting is always (re)done; the given comparer must define a strict
     * ordering between the vertexes.
     *
     * See Graph::VertexComparer and Graph::VertexFrequencyComparer for some already available comparers.
     */
    void rebuildForMining();

    template<typename Comparer>
//==============================================================================    
    void rebuildForMining(Comparer);

    /*
     * Let to trigger independently the rebuilding for mining except the final sorting of the adjacency lists.
     * It's exposed in order to use some complex comparers with the templatized form of rebuildForMining (e.g.
     * Graph::VertexFrequencyComparer), that requires calling this procedure in the graph before constructing
     * the comparer instance.
     *
     * Internally, all the forms of rebuildForMining() will not repeat the process done by
     * rebuildForMiningExceptSorting() if it was run previously, so no need to worry about waste processing.
     */
//==============================================================================     
    virtual void print_edge_map() const { /*Do Nothing*/ }
    
    virtual void rebuildForMiningExceptSorting();
    
    virtual WedgeMap* get_edge_map() const { return 0; }
//==============================================================================
    /*
     * Some comparers that can be used with the templatized form of rebuildForMining():
     *
     * - VertexComparer lets to sort the adjacency lists by increasing Vertex id.
     * - VertexFrequencyComparer lets to sort the adjacency lists by decreasing frequency of apparition in the
     *   adjacency lists (i.e. the inlinks count of each Vertex).
     *   To use it, Graph::rebuildForMiningExceptSorting() must have been called previously, otherwise it will trigger
     *   undefined behaviour.
     * - RandomVertexPermutationComparer matches each Vertex with another Vertex choosed randomly and then it sorts
     *   by increasing id of the matched Vertex. It lets to simulate that the ids were assigned originally to each
     *   Vertex in a random way, trying with it to provide a kind of 'middle point' to compare the impact of different
     *   sorting strategies in the quantity and quality of the mined dense subgraphs.
     *   To use it, Graph::rebuildForMiningExceptSorting() must have been called previously too.
     */
    struct VertexComparer;
    struct VertexFrequencyComparer;
    struct RandomVertexPermutationComparer;


    //// Inspectors ///////////////////////////////////////////////////////////////////////////////////////////////

    bool empty() const { return innerGraph.empty(); }

    std::size_t listsCount() const { return innerGraph.size(); }
    std::size_t nodesCount() const;     // It can be slow
    unsigned long arcsCount() const;

    bool isMineable() const { return mineability == 2; }
    bool isSortedByVertex() const { return sortedByVertex; }

    /*
     * Streamable brief summary of the graph, made to fit in one line.
     * It can be slow with huge social/web graphs, though.
     */
    friend std::ostream& operator<<(std::ostream&, const Graph&);

    /*
     * Print the graph via standard output. Accepted formats:
     *   - 0: Print one vertex with its respective adjacency list by each line. It's the default, and it's compatible
     *        with the format accepted by the Graph constructor from file.
     *   - 1: Print one arc (the pair of nodes defining it) by line.
     */
    void print(std::ostream&, unsigned int format=0) const;

    /*
     * Export the graph to a file. For now, only dumping as adjacency lists is supported, as it's the unique format
     * accepted back again by the Graph constructors.
     *
     * The extraSummary parameter let to add a short summary of the graph as a comment in the first line of the
     * output file. It uses the streamable summary form.
     */
    void dump(const std::string& fileName, bool extraSummary=false) const;  // It can throw an exception

//==============================================================================
    virtual int get_type(){ return UNDIRECTED_UNWEIGHTED; }
//==============================================================================

protected:

    std::map<Vertex, AdjacencyList> innerGraph;

    /*
     * Knowing than the adjacency lists are currently sorted by Vertex lets us improve the performance of the
     * implementation of many of our procedures. It's critical as this sorting is the only one
     * that will be used with big social/web graphs, who already come with their data sorted by id.
     */
    bool sortedByVertex;

    unsigned int mineability;   // 0: no mineable or unknown; 1: only sorting is pending; 2: mineable


    //// Internal helpers /////////////////////////////////////////////////////////////////////////////////////////

    typedef std::map<Vertex, AdjacencyList>::iterator iterator;

    iterator begin() { return innerGraph.begin(); }
    iterator   end() { return   innerGraph.end(); }

    //// Helpers for print() and dump()
    void printAsAdjacencyLists(std::ostream&) const;
    void printAsArcs(std::ostream&) const;
};


//// Graph constructors ///////////////////////////////////////////////////////////////////////////////////////////////

template<typename ContainerT>
Graph::Graph(const std::map<Vertex, ContainerT>& m): innerGraph(), sortedByVertex(false), mineability(0) {

    for (typename std::map<Vertex, ContainerT>::const_iterator it = m.begin(); it != m.end(); ++it) {
        innerGraph[it->first] = AdjacencyList(it->second.begin(), it->second.end());
    }
}


//// Graph comparers //////////////////////////////////////////////////////////////////////////////////////////////////

struct Graph::VertexComparer: public std::binary_function<Vertex, Vertex, bool> {
    bool operator()(Vertex vx1, Vertex vx2) const { return vx1 < vx2; }
};


struct Graph::VertexFrequencyComparer: public std::binary_function<Vertex, Vertex, bool> {
    explicit VertexFrequencyComparer(const Graph&);
    bool operator()(Vertex, Vertex) const;
private:
    std::map<Vertex, unsigned int> frequencies;
};


struct Graph::RandomVertexPermutationComparer: public std::binary_function<Vertex, Vertex, bool> {
    typedef std::map<Vertex, Vertex> Permutation;

    explicit RandomVertexPermutationComparer(const Graph&);
    bool operator()(Vertex, Vertex) const;
private:
    Permutation permutation;
};


//// Graph mutators ///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Comparer>
inline void
Graph::rebuildForMining(Comparer comparer) {
    rebuildForMiningExceptSorting();

    for (iterator it = begin(); it != end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), comparer);
    }
    sortedByVertex = false;
    mineability = 2;
}

template<>      // Template specialization for Graph::VertexComparer in order to update sortedByVertex
inline void
Graph::rebuildForMining<Graph::VertexComparer>(Graph::VertexComparer comparer) {
    rebuildForMiningExceptSorting();

    for (iterator it = begin(); it != end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), comparer);
    }
    sortedByVertex = true;
    mineability = 2;
}


}       // namespace odsg
#endif  // SRC_GRAPH_HPP_INCLUDED
