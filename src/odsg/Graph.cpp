#include "Graph.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <cstdlib>      // std::rand, std::srand
#include <set>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iterator>     // std::advance
#include <ctime>        // std::time

#include "utils/algorithms.hpp"
#include "utils/strings.hpp"
#include "GraphCluster.hpp"

namespace odsg {


namespace {     // Put here general, global definitions limited to this file

    std::ostream&
    operator<<(std::ostream& os, const Graph::AdjacencyList& outlinks) {
        for (Graph::AdjacencyList::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
            if (vxit != outlinks.begin())
                os << ' ';
            os << *vxit;
        }
        return os;
    }


    namespace graphFileFormat {     // Specific stuff to our text file format for graphs. TODO: Move to a class
        const char ADJACENCY_LIST_DELIMITER = ':';
        const char COMMENTS_TOKEN = '#';
    }


    inline void
    warnAboutInputFile(unsigned int line, const std::string& description, const std::string& resolution) {
        std::cerr << "warning: line " << line << " in input file: " << description << ": " << resolution
                  << std::endl;
    }

}   // namespace


//// Graph ////////////////////////////////////////////////////////////////////////////////////////////////////////////

Graph::Graph(): innerGraph(), sortedByVertex(true), mineability(2) {}


Graph::Graph(const std::string& fileName, bool comeSortedByVertex)
: innerGraph(), sortedByVertex(comeSortedByVertex), mineability(0) {

    assert(!fileName.empty());

    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("Graph::Graph(): can not open input graph file");
    }

    unsigned int lineCount = 0;
    std::string line;
    while (std::getline(infile, line)) {
        lineCount++;

        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r") == std::string::npos)
            continue;
        if (line[0] == graphFileFormat::COMMENTS_TOKEN)  // TODO: trim start of the line before checking
            continue;

        std::istringstream iss(line);

        Vertex vertex;
        if (!(iss >> vertex)) {
            warnAboutInputFile(lineCount, "no a vertex at the line start", "ignoring line");
            continue;
        }
        assert(innerGraph.find(vertex) == end());    // Not an already known vertex

        char delimiter = '\0';
        if (!(iss >> delimiter) || delimiter != graphFileFormat::ADJACENCY_LIST_DELIMITER) {
            warnAboutInputFile(lineCount, "missing or wrong delimiter", "ignoring line");
            continue;
        }
        AdjacencyList outlinks;
        for (Vertex vx; iss >> vx; ) {
            outlinks.push_back(vx);
        }
        assert(!comeSortedByVertex || algorithms::is_sorted(outlinks));
        assert(algorithms::has_unique(outlinks));

        // The remainings in the line are simply ignored

        innerGraph[vertex] = outlinks;
    }
}


Graph::Graph(const GraphCluster& cluster): innerGraph(), sortedByVertex(false), mineability(2) {

    for (GraphCluster::const_iterator it = cluster.begin(); it != cluster.end(); ++it) {
        innerGraph.insert(*it);
    }

    assert(listsCount() == cluster.listsCount());
}


void
Graph::rebuildForMining() {
    if (isMineable())
        return;

    // Unlike the templatized form of rebuildForMining(), here we try to skip the sorting whenever possible
    if (sortedByVertex) {
        rebuildForMiningExceptSorting();
        mineability = 2;
    } else {
        rebuildForMiningExceptSorting();    // Required by VertexFrequencyComparer
        rebuildForMining(VertexFrequencyComparer(*this));
    }

    assert(isMineable());
}


void
Graph::rebuildForMiningExceptSorting() {
    if (mineability >= 1) return;

    // Add self-loops and remove vertexes with 'trivial' (we need a better word) adjacency lists
    for (iterator it = begin(); it != end();  ) {
        Vertex vertex = it->first;
        AdjacencyList& outlinks = it->second;

        if (outlinks.size() >= 1 && !algorithms::is_found(outlinks, vertex)) {
            // If the graph was specified as already sorted by vertex, we don't want lose that property inserting
            // the self-loop without care.
            if (sortedByVertex) {
                outlinks.insert(std::upper_bound(outlinks.begin(), outlinks.end(), vertex), vertex);

                assert(algorithms::is_sorted(outlinks));
            } else {
                outlinks.push_back(vertex);
            }
        }
        if (outlinks.size() <= 1) {
            iterator tmp = it;
            ++it;
            innerGraph.erase(tmp);
        } else {
            ++it;
        }
    }

    // Finally, ensure that this procedure can't be done twice needlessly
    mineability = 1;
}


std::size_t
Graph::nodesCount() const {
    std::set<Vertex> nodes;

    for (const_iterator it = begin(); it != end(); ++it) {
        // The set nodes doesn't store duplicates, so we can to insert all without need to care
        nodes.insert(it->first);
        nodes.insert(it->second.begin(), it->second.end());
    }

    assert(nodes.size() >= listsCount());
    return nodes.size();
}


unsigned long
Graph::arcsCount() const {
    unsigned long arcs = 0;

    for (const_iterator it = begin(); it != end(); ++it) {
        arcs += it->second.size();      // Add the total of outlinks of each vertex
    }

    assert(!isMineable() || arcs >= listsCount() * 2);
    return arcs;
}


std::ostream&
operator<<(std::ostream& os, const Graph& graph) {
    // Comparing with Graph::print(), it ensures an short output, still with big graphs

    os << "Graph with "
       << graph.listsCount() << " adjacency lists, "
       << graph.nodesCount() << " nodes and "
       << graph.arcsCount() << " arcs"
       ;

    if (graph.isMineable())
        os << "; it was rebuilt for mining";

    return os;
}


void
Graph::print(std::ostream& os, unsigned int format) const {
    assert(format == 0 || format == 1);

    switch (format) {       // TODO: manage print formats as a class enum
        case 0: printAsAdjacencyLists(os);
                break;
        case 1: printAsArcs(os);
                break;
    }
}


void
Graph::dump(const std::string& fileName, bool extraSummary) const {
    assert(!fileName.empty());

    std::ofstream outfile(fileName.c_str());
    if (!outfile) {
        throw std::runtime_error("Graph::dump(): can not open output dump file");
    }

    if (extraSummary)
        outfile << graphFileFormat::COMMENTS_TOKEN << " " << *this << "\n\n";

    printAsAdjacencyLists(outfile);    // TODO: Support other formats in a good way
}


void
Graph::printAsAdjacencyLists(std::ostream& os) const {
    for (const_iterator it = begin(); it != end(); ++it) {
        os << it->first << ": " << it->second << '\n';
    }
}

void
Graph::printAsArcs(std::ostream& os) const {
    for (const_iterator it = begin(); it != end(); ++it) {
        Vertex vertex = it->first;
        const AdjacencyList& outlinks = it->second;

        for (Graph::AdjacencyList::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
            os << vertex << " " << *vxit << '\n';
        }
    }
}


//// VertexFrequencyComparer //////////////////////////////////////////////////////////////////////////////////////////

Graph::VertexFrequencyComparer::VertexFrequencyComparer(const Graph& graph): frequencies() {
    assert(graph.mineability >= 1);

    // Count and save the number of apparitions of each Vertex in all the adjacency lists
    for (Graph::const_iterator it = graph.begin(); it != graph.end(); ++it) {
        const Graph::AdjacencyList& outlinks = it->second;

        for (Graph::AdjacencyList::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
            frequencies[*vxit] += 1;
        }
    }
}


bool
Graph::VertexFrequencyComparer::operator()(Vertex vx1, Vertex vx2) const {
    assert(frequencies.find(vx1) != frequencies.end());
    assert(frequencies.find(vx2) != frequencies.end());

    unsigned int freq1 = frequencies.find(vx1)->second;
    unsigned int freq2 = frequencies.find(vx2)->second;

    if (freq1 != freq2)
        return freq1 > freq2;   // Decrescent orden by frequency...
    else
        return vx1 < vx2;   // ...then crescent orden by label
}


//// RandomVertexPermutationComparer //////////////////////////////////////////////////////////////////////////////////

Graph::RandomVertexPermutationComparer::RandomVertexPermutationComparer(const Graph& graph): permutation() {

    // First, save all vertexes that are part of the adjacency lists of the graph
    for (Graph::const_iterator it = graph.begin(); it != graph.end(); ++it) {
        const Graph::AdjacencyList& outlinks = it->second;

        for (Graph::AdjacencyList::const_iterator vxit = outlinks.begin(); vxit != outlinks.end(); ++vxit) {
            permutation[*vxit] = *vxit;
        }
    }
    assert(permutation.size() == graph.nodesCount());

    std::srand(std::time(0));   // Initialize seed

    // Random shuffling.
    // Loosely based in the 'possible implementation' for std::random_shuffle, from:
    //   http://en.cppreference.com/w/cpp/algorithm/random_shuffle
    // See also: Fisher-Yates shuffle
    Permutation::reverse_iterator it = permutation.rbegin();
    for (int i = permutation.size() - 1; i > 0; --i) {
        Permutation::reverse_iterator randit = it;
        std::advance(randit, std::rand() % (i + 1));

        std::swap(it->second, randit->second);
        ++it;
    }
}


bool
Graph::RandomVertexPermutationComparer::operator()(Vertex vx1, Vertex vx2) const {
    assert(permutation.find(vx1) != permutation.end());
    assert(permutation.find(vx2) != permutation.end());

    // The next is logically the same as:
    //     permutation[vx1] < permutation[vx2]
    return permutation.find(vx1)->second < permutation.find(vx2)->second;
}

}   // namespace odsg
