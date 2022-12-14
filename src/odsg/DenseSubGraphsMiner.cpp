#include "DenseSubGraphsMiner.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <cstddef>      // NULL

#include "utils/algorithms.hpp"
#include "Dag.hpp"
#include "DagNode.hpp"
#include "DenseSubGraph.hpp"
#include "DenseSubGraphsMaximalSet.hpp"
#include "MinerDagTraveler.hpp"
#include "MinerObjective.hpp"

namespace odsg {


DenseSubGraphsMiner::DenseSubGraphsMiner(const Dag* d,
                                         MinerDagTraveler* traveler,
                                         MinerObjective* objective,
                                         bool asCliques,
                                         unsigned long arcsCount)
: dag(d), parentTraveler(traveler), minerObjective(objective), asCliquesOnly(asCliques), minArcsCount(arcsCount) {

    assert(d);
    assert(traveler);
    assert(objective);
    assert(!asCliques || dynamic_cast<AsCliqueMinerObjective*>(objective));     // Use proper objective for cliques

    // Set the traveling paths inside the DagNodes objects for a faster retrieval later, during mining.
    // Note that it relies in no other DenseSubGraphsMiner instances being built for the same dag until that the
    // mining entrusted to the first instance is done.
    for (Dag::const_iterator it = dag->begin(); it != dag->end(); ++it) {
        const DagNode* node = *it;

        node->setNextTravelingNode(traveler->next(node));

        assert(node->getNextTravelingNode() != NULL || node->isRoot());
    }
}


DenseSubGraphsMaximalSet
DenseSubGraphsMiner::mine() const {
    DenseSubGraphsMaximalSet dagDSGs(asCliquesOnly);

    for (Dag::const_iterator it = dag->begin(); it != dag->end(); ++it) {
        const DagNode* const node = *it;

        if (willNotProvideEnoughGoodDsg(node))
            continue;

        // The next will find one big dense subgraph that includes to the current node in its centers set
        DenseSubGraph nodeDsg = getDenseSubGraphFrom(node);

        if (isNotGoodEnough(nodeDsg))
            continue;

        // Try to insert the dsg in the ongoing collection; if the dense subgraph is not maximal it will be ignored
        dagDSGs.insert(nodeDsg);
    }

    assert(dagDSGs.size() <= dag->nodesCount());
    return dagDSGs;
}

//==============================================================================
//VER ESTO!!
DenseSubGraph
DenseSubGraphsMiner::getDenseSubGraphFrom(const DagNode* node) const {
    assert(node != NULL);

    DenseSubGraph nodeDsg(node->getVertexes(), node->label);
    
    // For each node in the path, we check if adding it to the centers set provides us a 'better' dense subgraph
    
    const DagNode* pathNode = node;

    while ((pathNode = pathNode->getNextTravelingNode())) {     // Yes, it's an assignment
        DenseSubGraph candidateDsg(pathNode->getVertexes(), pathNode->label);
        candidateDsg.merge(nodeDsg);
        if (minerObjective->better(nodeDsg, candidateDsg)) {
            nodeDsg.swap(candidateDsg);
            if (minerObjective->best(nodeDsg))break;
        }
    }
    assert(algorithms::is_found(nodeDsg.getCenters(), node->label));
    assert(nodeDsg.getCenters().size() <= node->getMaxDepth());
    assert(nodeDsg.getSources().size() <= dag->nodesCount());

    return nodeDsg;

}
//==============================================================================

bool
DenseSubGraphsMiner::willNotProvideEnoughGoodDsg(const DagNode* node) const {
    assert(node != NULL);

    return node->isRoot();
}


bool
DenseSubGraphsMiner::isNotGoodEnough(const DenseSubGraph& dsg) const {
    assert(dsg.getCenters().size() > 0);

    return dsg.getCenters().size() == 1 || dsg.getSources().size() <= 1 || dsg.arcsCount() < minArcsCount;
}


}   // namespace odsg
