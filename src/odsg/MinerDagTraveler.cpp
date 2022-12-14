#include "MinerDagTraveler.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <cstddef>      // NULL, std::size_t
#include <vector>

#include "utils/algorithms.hpp"
#include "DagNode.hpp"

namespace odsg {


const DagNode*
ParentTraveler::next(const DagNode* node) {
    const DagNode* firstParent = NULL;

    if (!node->getParents().empty())
        firstParent = node->getParents().front();

    return firstParent;
}


const DagNode*
DeepestParentTraveler::next(const DagNode* node) {
    const DagNode* deepestParent = NULL;

    for (std::vector<const DagNode*>::const_iterator pit = node->getParents().begin();
         pit != node->getParents().end();
         ++pit) {

        if ((*pit)->getMaxDepth() == node->getMaxDepth() - 1) {
            deepestParent = *pit;    // Many parents could fulfill the condition, we take the first
            break;
        }
    }
    return deepestParent;
}


const DagNode*
SharingMoreVertexesParentTraveler::next(const DagNode* node) {
    const DagNode* parentSharingMoreVertexes = NULL;
    std::size_t maxSharedCount = 0;

    for (std::vector<const DagNode*>::const_iterator pit = node->getParents().begin();
         pit != node->getParents().end();
         ++pit) {

        std::size_t sharedCount = algorithms::set_intersection_count((*pit)->getVertexes(), node->getVertexes());
        assert(sharedCount >= 1);

        if (maxSharedCount < sharedCount) {
            maxSharedCount = sharedCount;
            parentSharingMoreVertexes = *pit;
        }
        if (sharedCount == node->getVertexes().size())
            break;      // Many parents could fulfill the requeriments, we take the first
    }
    return parentSharingMoreVertexes;
}


}   // namespace odsg
