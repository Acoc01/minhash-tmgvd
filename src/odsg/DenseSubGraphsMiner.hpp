#ifndef SRC_DENSE_SUB_GRAPHS_MINER_HPP_INCLUDED
#define SRC_DENSE_SUB_GRAPHS_MINER_HPP_INCLUDED

namespace odsg {


class Dag;
class DagNode;
class DenseSubGraph;
class DenseSubGraphsMaximalSet;
class MinerDagTraveler;
class MinerObjective;

/*
 * Define the general, basic, approach to mining dense subgraphs from a dag.
 *
 * The lifetime of the arguments passed as pointers must extend past any use of the respective DenseSubGraphsMiner
 * object.
 */
class DenseSubGraphsMiner {
public:
    DenseSubGraphsMiner(const Dag*, MinerDagTraveler*, MinerObjective*, bool asCliques, unsigned long arcsCount);

    DenseSubGraphsMaximalSet mine() const;

private:
    const Dag* const dag;
    MinerDagTraveler* const parentTraveler;     // Currently only used during construction
    MinerObjective* const minerObjective;
    bool asCliquesOnly;
    unsigned long minArcsCount;

    // Helpers for mine()
    DenseSubGraph getDenseSubGraphFrom(const DagNode*) const;

    bool willNotProvideEnoughGoodDsg(const DagNode*) const;
    bool isNotGoodEnough(const DenseSubGraph&) const;
};


}       // namespace odsg
#endif  // SRC_DENSE_SUB_GRAPHS_MINER_HPP_INCLUDED
