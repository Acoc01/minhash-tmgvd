#ifndef SRC_DAG_HPP_INCLUDED
#define SRC_DAG_HPP_INCLUDED

#include <cstddef>      // NULL, std::size_t
#include <vector>
#include <map>
#include <iosfwd>

#include "DagNode.hpp"          // All of our container-like classes include the definition of the contained element
#include "Vertex.hpp"

namespace odsg {


class DenseSubGraphsMaximalSet;
class Graph;
class GraphCluster;

/*
 * A Dag object is a collection of DagNode objects linked between them, from where dense subgraphs are mined.
 * It's an acyclic directed graph, though in this project the word 'graph' alone is always assumed to apply to Graph
 * objects and not to Dag objects. Following the same idea, a different vocabulary for the relations between nodes
 * of the Dag is used here to reduce confusion (e.g. 'node' instead of 'vertex'; 'child' instead of 'outlink'),
 * being it similar to that seen in other contexts dealing with DAGs, as compilers and version control systems.
 */
class Dag {
public:

    //// Types ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector<const DagNode*>::const_iterator const_iterator;


    //// Constructors & Destructor ////////////////////////////////////////////////////////////////////////////////

    /*
     * All Dag constructor take as argument a different type of specification of a graph (the 'source graph', as
     * it's called in some places); it's expected that it's mineable: Graph objects have methods to check/ensure it,
     * and all cluster of a mineable graph is OK too. Currently, Dag objects are created only as part of the
     * process of construction of DagForest objects, where mineability of the input graphs is ensured.
     *
     * Clusters doesn't know much about the graph from whom they were built, so the comeSortedByVertex property
     * needs be indicated by the caller in these cases to get a performance boost during construction.
     * TODO: It sucks, but I'm in a hurry now.
     * In case of doubt, don't set it, otherwise it will trigger undefined behaviour.
     */
    explicit Dag(const Graph&);
    explicit Dag(const GraphCluster& cluster, bool comeSortedByVertex=false);

    ~Dag();


    //// No public mutators: a dag isn't altered outside of the contructors & destructor //////////////////////////


    //// Iterators ////////////////////////////////////////////////////////////////////////////////////////////////

    const_iterator begin() const { return nodeCache.begin(); }
    const_iterator end() const { return nodeCache.end(); }


    //// Inspectors ///////////////////////////////////////////////////////////////////////////////////////////////

    bool empty() const { return nodeCache.empty(); }

    std::size_t nodesCount() const { return nodeCache.size(); }
    unsigned long arcsCount() const;

    /*
     * Mine a collection of maximal dense subgraphs from the dag.
     *
     * The objective argument lets to choose from a list of predefined MinerObjective implementations:
     *   - 0: Maximize the size of the centers set while it is kept as a subset of the sources, favoring the mining
     *        of cliques.
     *   - 1: Maximize the number of arcs.
     *   - 2: Maximize the number of shared elements between sources and centers.
     *
     * If the asCliquesOnly option is set, for all the mined dense subgraphs, their centers can be directly treated
     * as a clique (the vertexes are present in the sources set too). Note that for these cases minArcsCount is not
     * for the arcs of the clique itself. Also it forces the use of the objective number 0, regardless of the value
     * passed as the objective argument, as currently it's the only compatible with it.
     */
    DenseSubGraphsMaximalSet getDenseSubGraphs(unsigned int traveler,
                                               unsigned int objective,
                                               bool asCliquesOnly,
                                               unsigned long minArcsCount) const;


    friend std::ostream& operator<<(std::ostream&, const Dag&);     // Short summary, made to fit in one line
    void print(std::ostream&) const;


private:

    /*
     * All the nodes in the dag, saved to can iterate linearly over them. During initialization, this vector is
     * sorted in such a way that, while iterating over it, no node is visited before all its parents, i.e. a
     * 'topological sorting', as it's known in the context of DAGs.
     */
    std::vector<const DagNode*> nodeCache;

    /*
     * All those nodes without parents. The mining algorithms doesn't get really favored having here these nodes
     * cached, but the overhead is fairly low, so it's not a problem.
     */
    std::vector<const DagNode*> roots;

    /*
     * The maximum value seen of a node in the dag for its 'maxDepth' property. It, then, corresponds to the length
     * of the largest possible path between two nodes (plus one): the so called in the literature 'diameter' property
     * for graphs.
     *
     * It's cached only for informational purposes, as the complexity of our current algorithms for mining dense
     * subgraphs depends, in some way, of this value.
     */
    unsigned int maxNodeMaxDepth;

    //// Internal helpers /////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector<const DagNode*>::iterator iterator;

    iterator begin() { return nodeCache.begin(); }
    iterator end() { return nodeCache.end(); }


    //// Helpers for the constructors

    bool fromGraphSortedByVertex;

    template<typename GraphT>
    void initialize(const GraphT&);

    void insert(Vertex, const std::vector<Vertex>&, std::map<Vertex, DagNode*>&);
    void setTopologicalCacheSorting(const std::map<Vertex, DagNode*>&);
    void updateNodeMaxDepths();

//==============================================================================     
    const Graph* const ptr_graph; 
    /*
     * Keeps reference to "Graph with weight"
    */ 
//==============================================================================    


    //// Others ///////////////////////////////////////////////////////////////////////////////////////////////////

    // The next two are declared and deliberately NOT implemented, to prevent copying objects of this classs
    Dag(const Dag&);
    Dag& operator=(const Dag&);
};


}       // namespace odsg
#endif  // SRC_DAG_HPP_INCLUDED
