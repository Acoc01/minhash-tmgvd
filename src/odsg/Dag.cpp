#include "Dag.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <ostream>
#include <algorithm>
#include <memory>       // std::auto_ptr

#include "utils/algorithms.hpp"
#include "DenseSubGraphsMiner.hpp"
#include "DenseSubGraphsMaximalSet.hpp"
#include "Graph.hpp"
//==============================================================================
#include "WGraph.hpp"
#include <iostream>
//==============================================================================
#include "GraphCluster.hpp"
#include "MinerDagTraveler.hpp"
#include "MinerObjective.hpp"

namespace odsg {


namespace {     // Put here general, global definitions limited to this file

    std::ostream&
    operator<<(std::ostream& os, const std::vector<const DagNode*>& nodes) {
        for (std::vector<const DagNode*>::const_iterator nit = nodes.begin(); nit != nodes.end(); ++nit) {
            if (nit != nodes.begin())
                os << ' ';
            os << (*nit)->label;
        }
        return os;
    }

}   // namespace


//// Dag //////////////////////////////////////////////////////////////////////////////////////////////////////////////

Dag::Dag(const Graph& graph)
: nodeCache(), roots(), maxNodeMaxDepth(0), fromGraphSortedByVertex(graph.isSortedByVertex()),
  ptr_graph(&graph) {
    assert(graph.isMineable());

    initialize(graph);
    //const WGraph* const wgraph = (const WGraph*)(&graph);
    //wgraph->print_edge_map();

}

Dag::Dag(const GraphCluster& cluster, bool comeSortedByVertex)
: nodeCache(), roots(), maxNodeMaxDepth(0), fromGraphSortedByVertex(comeSortedByVertex),
  ptr_graph(cluster.get_ptrGraph()) {

    initialize(cluster);

}


template<typename GraphT>
void
Dag::initialize(const GraphT& graph) {
    // Only for the creation of the nodes, we choose to use a map to store the existent nodes, because the fast
    // retrieval capacity and simpler syntax. After finish it, its contents are passed to nodeCache and then the
    // map is discarted.
    {
        std::map<Vertex, DagNode*> tmpNodeMapCache;

        for (typename GraphT::const_iterator it = graph.begin(); it != graph.end(); ++it) {
            insert(it->first, it->second, tmpNodeMapCache);
        }
        assert(tmpNodeMapCache.size() == graph.nodesCount());

        setTopologicalCacheSorting(tmpNodeMapCache);
    }

    // As the insertion of one node can change the maxDepth property for all the descendent nodes in the dag, setting
    // this property for all the nodes can be done just after ending all the insertions.
    updateNodeMaxDepths();


    assert(roots.size() <= graph.listsCount());
    assert(roots.size() <= nodesCount() - maxNodeMaxDepth + 1);
    assert(std::count_if(roots.begin(), roots.end(),            // roots don't have parents
                         DagNode::countParentsComparer(0)) == roots.size());
    assert(std::count_if(nodeCache.begin(), nodeCache.end(),    // Only roots don't have parents
                         DagNode::countParentsComparer(0)) == roots.size());
}


void
Dag::insert(Vertex vertex, const std::vector<Vertex>& outlinks, std::map<Vertex, DagNode*>& nodeMapCache) {
    assert(algorithms::is_found(outlinks, vertex));     // Self-loops are present
    assert(outlinks.size() > 1);                        // No trivial outlinks

    DagNode* prevNode = NULL;
    for (std::vector<Vertex>::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
        Vertex outlink = *vxit;

        DagNode* node = nodeMapCache[outlink];      // There will not be two nodes with the same label in the dag
        if (!node)
            node = new DagNode(outlink);

        // Update relations between nodes
        node->addVertex(vertex);
        if (prevNode)
            prevNode->addChild(node);

        // Update general components of the dag
        if (!prevNode && !nodeMapCache[outlink])
            roots.push_back(node);  // Any unknown node found at start of a adjacency list is saved as root
        if (prevNode) {
            std::vector<const DagNode*>::iterator rit = std::find(roots.begin(), roots.end(), node);
            if (rit != roots.end())
                roots.erase(rit);   // Any known root found not at start of a adjacency list is not root anymore
        }
        nodeMapCache[outlink] = node;

        // Other preparations for the next iteration
        prevNode = node;
    }
}


void
Dag::setTopologicalCacheSorting(const std::map<Vertex, DagNode*>& nodeMapCache) {
    nodeCache.reserve(nodeMapCache.size());

    if (fromGraphSortedByVertex) {
        // When the adjacency lists of the source graph are sorted by Vertex, *one* possible topological sorting is
        // simply the lineal listing of the nodes by increasing label. It's a great help with huge social/web graphs.
        for (std::map<Vertex, DagNode*>::const_iterator it = nodeMapCache.begin(); it != nodeMapCache.end(); ++it) {
            nodeCache.push_back(it->second);
        }
    } else {
        // The next general implementation to find a topological sorting is fairly inneficient, but very easy to
        // understand too. Given that it is never applied to huge social/web graphs, it can be good enough.
        // To see a better implementation, check Cormen's 'Introduction to Algorithms', 3rd Ed., section 22.4.
        while (nodeCache.size() < nodeMapCache.size()) {
            const DagNode* nextNode = NULL;

            for (std::map<Vertex, DagNode*>::const_iterator it = nodeMapCache.begin(); it != nodeMapCache.end(); ++it) {
                const DagNode* const node = it->second;

                if (!algorithms::is_found(nodeCache, node) && algorithms::are_found(nodeCache, node->parents)) {
                    nextNode = node;
                    break;
                }
            }
            assert(nextNode != NULL);

            nodeCache.push_back(nextNode);
        }
    }

    assert(nodeCache.size() == nodeMapCache.size());
    assert(algorithms::has_unique(nodeCache));
}


void
Dag::updateNodeMaxDepths() {
    // As the maxDepth value of each node depends only of the maxDepth values of its parents, the topological sorting
    // of nodeCache let us to perform this task without resorting to recursion.

    for (iterator it = begin(); it != end(); ++it) {
        DagNode* const node = const_cast<DagNode*>(*it);    // Remove constness

        unsigned int maxParentMaxDepth = 0;
        for (std::vector<const DagNode*>::const_iterator pit = node->parents.begin();
             pit != node->parents.end();
             ++pit) {

            if (maxParentMaxDepth < (*pit)->maxDepth)
                maxParentMaxDepth = (*pit)->maxDepth;
        }
        node->maxDepth = maxParentMaxDepth + 1;

        if (node->isLeaf() && node->maxDepth > maxNodeMaxDepth)
            maxNodeMaxDepth = node->maxDepth;

        assert(node->maxDepth >= 1);
        assert(node->maxDepth <= nodesCount() - roots.size() + 1);      // Theoretical upper limit for maxDepth
    }

    assert(std::count_if(roots.begin(), roots.end(),            // roots have maxDepth == 1
                         DagNode::maxDepthComparer(1)) == roots.size());
    assert(std::count_if(nodeCache.begin(), nodeCache.end(),    // Only roots have maxDepth == 1
                         DagNode::maxDepthComparer(1)) == roots.size());
    assert(std::count_if(nodeCache.begin(), nodeCache.end(),    // At most one node have the theoretical max. maxDepth
                         DagNode::maxDepthComparer(nodesCount() - roots.size() + 1)) <= 1);

    assert((maxNodeMaxDepth == 0) == nodeCache.empty());
}


Dag::~Dag() {
    for (iterator it = begin(); it != end(); ++it) {
        delete *it;
    }
}


unsigned long
Dag::arcsCount() const {    // TODO: Have it cached after first use
    unsigned long arcs = 0;

    for (const_iterator it = begin(); it != end(); ++it) {
        arcs += (*it)->parents.size();      // Counting parents or children is the same
    }

    assert(arcs >= nodesCount() - roots.size());    // Theoretical lower limit
                                                    // Theoretical upper limit can do overflow easily
    return arcs;
}


DenseSubGraphsMaximalSet
Dag::getDenseSubGraphs(unsigned int traveler,
                       unsigned int objective,
                       bool asCliquesOnly,
                       unsigned long minArcsCount) const {

    assert(traveler == 0 || traveler == 1);
    assert(objective <= 7);


    MinerDagTraveler* minerTravelerPtr = NULL;
    switch (traveler) {
        case 0:
            minerTravelerPtr = new DeepestParentTraveler;
            break;
        case 1:
            minerTravelerPtr = new SharingMoreVertexesParentTraveler;
            break;
    }
    std::auto_ptr<MinerDagTraveler> minerTraveler(minerTravelerPtr);    // It takes care doing 'delete' on the object


    MinerObjective* minerObjectivePtr = NULL;
    if (asCliquesOnly) {
        // For the moment, only this MinerObjective is compatible with the option asCliquesOnly.
        // Others will provide wrong results regarding the maximality of the mined collection!
        minerObjectivePtr = new AsCliqueMinerObjective;
    } else {
        WedgeMap *wedge_map = ptr_graph->get_edge_map();
        switch (objective) {
            case 0:
                minerObjectivePtr = new AsCliqueMinerObjective;
                break;
            case 1:
                minerObjectivePtr = new LegacyMinerObjective;
                break;
            case 2:
                minerObjectivePtr = new MaxIntersectionMinerObjective;
                break;

//==============================================================================
            case 3:
                minerObjectivePtr = new SimpleEdgeDensity(wedge_map);
                break;
            case 4:
                minerObjectivePtr = new SimpleDegreeDensity(wedge_map);
                break;
            case 5:
                minerObjectivePtr = new DegreeAndEdgeDensity(wedge_map);
                break;
            case 6:
                minerObjectivePtr = new FullEdgeDensity(wedge_map);
                break;
            case 7:
                minerObjectivePtr = new FullDegreeDensity(wedge_map);
                break;
//==============================================================================
        }
    }
    std::auto_ptr<MinerObjective> minerObjective(minerObjectivePtr);    // It takes care doing 'delete' on the object


    return DenseSubGraphsMiner(this, minerTraveler.get(), minerObjective.get(), asCliquesOnly, minArcsCount)
            .mine();
}


std::ostream&
operator<<(std::ostream& os, const Dag& dag) {
    // Comparing with Dag::print(), it ensures an short output, still with big dags

    os << "Dag with "
       << dag.nodesCount() << " nodes, "
       << dag.arcsCount() << " arcs, "
       << dag.roots.size() << " roots and "
       << "maximum maxDepth = " << dag.maxNodeMaxDepth;

    return os;
}


void
Dag::print(std::ostream& os) const {
    os << "The dag has "
       << nodesCount() << " nodes, "
       << arcsCount() << " arcs, "
       << roots.size() << " roots {" << roots << "} and "
       << "maximum maxDepth = " << maxNodeMaxDepth << ". "
       << "The nodes are:\n";

    for (const_iterator it = begin(); it != end(); ++it) {
        os << "  ";
        (*it)->print(os);
        os << '\n';
    }
}


}   // namespace odsg
