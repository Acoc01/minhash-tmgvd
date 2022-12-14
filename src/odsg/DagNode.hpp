#ifndef SRC_DAGNODE_HPP_INCLUDED
#define SRC_DAGNODE_HPP_INCLUDED

#include <vector>
#include <iosfwd>

#include "Vertex.hpp"
#include "VertexSet.hpp"

namespace odsg {


/*
 * DagNode objects are expected to be created and manipulated only from the Dag class, and inspected by classes that
 * iterate over these Dag objects. The notable exception are DenseSubGraphsMiner objects that query and manipulate
 * a single caching attribute by performance concerns.
 *
 * DagNode could have been defined as nested in Dag, but dealing with nested classes can be tricky when it's combined
 * with friendship and other things, so, for now, it is an independent class, with Dag as a friend class.
 */
class DagNode {
public:
    explicit DagNode(Vertex);

    friend class Dag;

    // No public mutators

    // Inspectors
    const VertexSet& getVertexes() const { return vertexes; }                       //
    const std::vector<const DagNode*>& getChildren() const { return children; }     // Const-references to internals
    const std::vector<const DagNode*>& getParents() const { return parents; }       //
    unsigned int getMaxDepth() const { return maxDepth; }                           //

    bool isChildOf(const DagNode*) const;
    bool isParentOf(const DagNode* node) const { return node->isChildOf(this); }

    bool isRoot() const { return parents.empty(); }
    bool isLeaf() const { return children.empty(); }

    void print(std::ostream&) const;


    // The next are for use only by DenseSubGraphsMiner objects: getting from the nodes themselves the traveling path
    // let us improve the performance of the whole mining process.
    const DagNode* getNextTravelingNode() const { return cachedNextTravelingNode; }
    void setNextTravelingNode(const DagNode* node) const { cachedNextTravelingNode = node; }


    const Vertex label;

private:
    VertexSet vertexes;             // Labels of the inlinks of this vertex in the original source graph
    std::vector<const DagNode*> children;
    std::vector<const DagNode*> parents;
    unsigned int maxDepth;          // Length of the largest path between any root and this node, plus one


    mutable const DagNode* cachedNextTravelingNode;

    // Mutators
    void addVertex(Vertex);
    void addChild(DagNode*);
    void addParent(DagNode* node) { node->addChild(this); }


    // The next two are declared and deliberately NOT implemented, to prevent copying objects of this class
    DagNode(const DagNode&);
    DagNode& operator=(const DagNode&);


    // The next predicates are used only by assertions
#ifndef NDEBUG
    struct labelComparer {
        explicit labelComparer(Vertex vx): label(vx) {}
        bool operator()(const DagNode* node) const { return label == node->label; }
    private:
        const Vertex label;
    };

    struct countParentsComparer {
        explicit countParentsComparer(unsigned int n): parentsCount(n) {}
        bool operator()(const DagNode* node) const { return parentsCount == node->parents.size(); }
    private:
        const unsigned int parentsCount;
    };

    struct maxDepthComparer {
        explicit maxDepthComparer(unsigned int d): depth(d) {}
        bool operator()(const DagNode* node) const { return depth == node->maxDepth; }
    private:
        const unsigned int depth;
    };
#endif  // NDEBUG

};


}       // namespace odsg
#endif  // SRC_DAGNODE_HPP_INCLUDED
