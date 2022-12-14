#include "DenseSubGraph.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <ostream>
#include <istream>

#include "utils/algorithms.hpp"

namespace odsg {


DenseSubGraph::DenseSubGraph(const VertexSet& ss, const VertexSet& cc): sources(ss), centers(cc), density_value(0.0f) {}

DenseSubGraph::DenseSubGraph(const VertexSet& ss, Vertex c): sources(ss), centers(), density_value(0.0f) {
    centers.insert(c);
}


bool
DenseSubGraph::operator==(const DenseSubGraph& dsg) const {
    return sources == dsg.sources && centers == dsg.centers;
}
bool
DenseSubGraph::operator!=(const DenseSubGraph& dsg) const {
    return !(*this == dsg);
}


DenseSubGraph&
DenseSubGraph::merge(const DenseSubGraph& dsg) {
    // A merge of two dense subgraphs is built intersecting both sources sets and joining both centers sets

    VertexSet tmp;
    algorithms::set_intersection(sources, dsg.sources, tmp);
    sources.swap(tmp);

    centers.insert(dsg.centers.begin(), dsg.centers.end());

    assert(sources.size() <= dsg.sources.size() && sources.size() <= tmp.size());
    assert(centers.size() >= dsg.centers.size());

    return *this;   // Let concatenate merge() calls
}


void
DenseSubGraph::swap(DenseSubGraph& dsg) {
    sources.swap(dsg.sources);
    centers.swap(dsg.centers);
//==============================================================================    
    float temp_value  =     density_value;
    density_value     = dsg.density_value;
    dsg.density_value =        temp_value;
//==============================================================================
}

//==============================================================================
void
DenseSubGraph::set_density_value(float value) {
   density_value = value;
}
//==============================================================================
bool
DenseSubGraph::includes(const DenseSubGraph& subDsg) const {
    return algorithms::set_includes(sources, subDsg.sources) &&
           algorithms::set_includes(centers, subDsg.centers);
}


unsigned long
DenseSubGraph::arcsCount() const {
    return sources.size() * centers.size();
}

unsigned long
DenseSubGraph::nodesCount() const {
    VertexSet nodes = sources;
    nodes.insert(centers.begin(), centers.end());
    return nodes.size();
}

bool
DenseSubGraph::clique() const {
    return sources == centers;
}

bool
DenseSubGraph::asClique() const {
    // The next is the same as:
    //      return !clique() && algorithms::set_includes(sources, centers));
    // ...but a bit faster.
    return sources.size() != centers.size() && algorithms::set_includes(sources, centers);
}

bool
DenseSubGraph::biClique() const {
    return algorithms::set_intersection_count(sources, centers) == 0;
}


std::string
DenseSubGraph::description() const {
    if (clique())
        return "clique";
    if (asClique())
        return "as-clique";
    if (biClique())
        return "biclique";

    assert(generic());      // Keep this method consistent with related inspector methods
    return "generic";
}


//==============================================================================
float
DenseSubGraph::get_density_value() const { 
   return density_value;
}
//==============================================================================

std::istream&
operator>>(std::istream& is, DenseSubGraph& dsg) {
    VertexSet cc, ss;

    is >> cc;

    std::string delimiter;
    is.width(std::string("<---").length());     // Max. number of characters to read in the next call
    if (is >> delimiter && delimiter == "<---") {
        is >> ss;
    } else {    // A clique, that is represented as a set
        is.clear();     // TODO: put back the string extracted in the stream too if it was something unknown!
        ss = cc;
    }

    // TODO: The passed-by-reference dsg parameter should be modified only if everything went well
    dsg.sources.swap(ss);
    dsg.centers.swap(cc);

    return is;
}


std::ostream&
operator<<(std::ostream& os, const DenseSubGraph& dsg) {
    os << dsg.centers;

    if (dsg.sources != dsg.centers || dsg.empty()) {
        if (!dsg.centers.empty())
            os << ' ';

        os << "<---";

        if (!dsg.sources.empty())
            os << ' ' << dsg.sources;
    }

    return os;
}


}   // namespace odsg
