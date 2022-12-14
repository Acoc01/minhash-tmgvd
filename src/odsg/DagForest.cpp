#include "DagForest.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <stdexcept>
#include <fstream>
#include <ostream>
#include <memory>       // std::auto_ptr
//==============================================================================
#include<iostream>
//==============================================================================



#include "Graph.hpp"
#include "GraphCluster.hpp"
#include "GraphPartitioner.hpp"

namespace odsg {


DagForest::DagForest(const          Graph& graph,
                     int        clusteringScheme,
                     unsigned int minClusterSize,
                     bool sortClusterByFrequency)
: forest() {

    assert(clusteringScheme == 0 || clusteringScheme == 1 || clusteringScheme == 2);

    if (graph.empty())
        return;     // A empty graph will lead to a empty forest: hardly useful but allowed

    if (!graph.isMineable()) {
        throw std::logic_error("DagForest::DagForest(): graph must be mineable");
    }

    std::vector<GraphCluster> partition;

    GraphPartitioner* partitionerPtr = NULL;
    switch (clusteringScheme) {
        case 0:     // No partitioning
            break;
        case 1:     // Partitioning by common initial outlink
            partitionerPtr = new GraphPartitionerByInitialOutlink(&graph);
            break;
        case 2:     // Partitioning by common hashing signature
            partitionerPtr = new GraphPartitionerBySignature(&graph);
            break;
    }
    if (partitionerPtr) {
        std::auto_ptr<GraphPartitioner> partitioner(partitionerPtr);    // It takes care of delete the object

        GraphCluster cluster = partitioner->next(minClusterSize);
        while (!cluster.empty()) {
            partition.push_back(cluster);
            cluster = partitioner->next(minClusterSize);
        }
    }

    // The partition is empty when no partitioning was set
    if (partition.empty()) {
        forest.push_back(new Dag(graph));
    } else {
        for (std::vector<GraphCluster>::const_iterator it = partition.begin(); it != partition.end(); ++it) {

            if (sortClusterByFrequency) {
                const GraphCluster& cluster = *it;

                // Due to the current overall workflow to build dags from graphs (sorting followed by clustering),
                // to support without much pain this added-in-final-stages sortClusterByFrequency option (that
                // requires clustering followed by sorting) it's inevitable to duplicate some data, maybe doing it
                // not very suitable for huge social-web graphs.
                Graph clusterGraph(cluster);        // It create a copy of the data
                clusterGraph.rebuildForMiningExceptSorting();   // Required by VertexFrequencyComparer
                clusterGraph.rebuildForMining(Graph::VertexFrequencyComparer(clusterGraph));

                forest.push_back(new Dag(clusterGraph));
            } else {
                forest.push_back(new Dag(*it, graph.isSortedByVertex()));
            }
        }
    }

    assert(size() <= graph.listsCount());
}

DagForest::~DagForest() {
    for (iterator tit = begin(); tit != end(); ++tit) {
        delete *tit;
    }
}


std::ostream&
operator<<(std::ostream& os, const DagForest& forest) {
    os << "The dag forest has " << forest.size() << " dags";
    return os;
}


void
DagForest::print(std::ostream& os, bool onlySummaries) const {
    os << *this;

    if (!empty()) {
        os << ":\n";

        for (const_iterator tit = begin(); tit != end(); ++tit) {
            if (tit != begin())
                os << '\n';

            if (onlySummaries)
                os << **tit;
            else
                (*tit)->print(os);
        }
    }
}


void
DagForest::dump(const std::string& fileName, bool onlySummaries) const {
    assert(!fileName.empty());

    std::ofstream outfile(fileName.c_str());
    if (!outfile) {
        throw std::runtime_error("DagForest::dump(): can not open output dump file");
    }

    print(outfile, onlySummaries);
}


}   // namespace odsg
