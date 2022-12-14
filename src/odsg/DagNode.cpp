#include "DagNode.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <cstddef>      // NULL, std::size_t
#include <ostream>
#include <algorithm>

#include "utils/algorithms.hpp"

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


//// DagNode //////////////////////////////////////////////////////////////////////////////////////////////////////////

DagNode::DagNode(Vertex vx)
: label(vx), vertexes(), children(), parents(), maxDepth(1), cachedNextTravelingNode(NULL) {}


void
DagNode::addVertex(Vertex vertex) {
    vertexes.insert(vertex);
}


void
DagNode::addChild(DagNode* node) {
    assert(node != NULL);
    assert(isParentOf(node) == node->isChildOf(this));  // Symmetry in the parent/child relationship is satisfied

    if (isParentOf(node)) return;

    children.push_back(node);
    node->parents.push_back(this);

    assert(std::count(children.begin(), children.end(), node) == 1);            // No duplication in children
    assert(std::count(node->parents.begin(), node->parents.end(), this) == 1);  // No duplication in parents
}


bool
DagNode::isChildOf(const DagNode* node) const {
    assert(node != NULL);
    assert(std::count(node->children.begin(), node->children.end(), this)   // Search by ptr or label must be the same
           == std::count_if(node->children.begin(), node->children.end(), labelComparer(this->label)));

    return algorithms::is_found(node->children, this);
}


void
DagNode::print(std::ostream& os) const {
    os << label << " has "
       << vertexes.size() << " vertexes {" << vertexes << "}, "
       << parents.size() << " parents {" << parents << "}, "
       << children.size() << " children {" << children << "} and "
       << "maxDepth = " << maxDepth;
}


}   // namespace odsg
