#ifndef SRC_MINER_DAG_TRAVELER_HPP_INCLUDED
#define SRC_MINER_DAG_TRAVELER_HPP_INCLUDED

namespace odsg {


class DagNode;

/*
 * MinerDagTraveler objects let to customize one aspect of the mining process from dags done by DenseSubGraphsMiner
 * objects: how to move from each node to other in a Dag object according to some criteria.
 *
 * The determination of the next node depends only of the data of the current node (a.k.a. no context-aware of the
 * mining in process), so multiples queries for the same node always gives the same result, letting then to have it
 * cached inside the DagNode objects themselves for the fastest possible retrieval during the mining process.
 *
 * The implicit path resultant from a starting node to some ending root is called in some places a 'mining path'.
 * All implementations must ensure at least no cyclic mining paths.
 *
 * There isn't here an explicit reference to the Dag object itself; it's assumed that all queries to a traveler
 * will be for nodes from the same unique dag.
 *
 * MinerDagTraveler is an abstract class. Derived classes needs to define how each next node is selected.
 */
class MinerDagTraveler {
public:
    virtual ~MinerDagTraveler() {}

    virtual const DagNode* next(const DagNode*) = 0;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Move from each node to 'some' parent: the most basic idea. Not very useful, though.
 * Used as comparison base for some analysis, nothing more.
 */
class ParentTraveler : public MinerDagTraveler {
public:
    /*virtual*/ const DagNode* next(const DagNode*);
};


/*
 * Move from each node to a parent with the maximum value for its maxDepth property (always lower by one that the
 * maxDepth of the node itself, by definition). It allows, from any initial node, take a path to some root so that
 * the length of the path is maximal.
 *
 * Between all travelers here defined, it provides the best results for mining.
 */
class DeepestParentTraveler : public MinerDagTraveler {
public:
    /*virtual*/ const DagNode* next(const DagNode*);
};


/*
 * Move from each node to a parent with whom it shares the largest number of vertexes; it could be useful in
 * conjuntion with MaxIntersectionMinerObjective.
 */
class SharingMoreVertexesParentTraveler : public MinerDagTraveler {
public:
    /*virtual*/ const DagNode* next(const DagNode*);
};


}       // namespace odsg
#endif  // SRC_MINER_DAG_TRAVELER_HPP_INCLUDED
