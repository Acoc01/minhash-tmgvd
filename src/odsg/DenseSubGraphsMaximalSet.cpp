#include "DenseSubGraphsMaximalSet.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <stdexcept>
#include <fstream>
#include <ios>          // std::ios::app, std::ios::out

#include "utils/algorithms.hpp"

namespace odsg {


/*
 * TODO: This procedure is slow. We need to investigate in the literature how to improve it, probably using a custom
 * data structure to store the dense subgraphs. By example, maybe having indexed by the sizes of sources and centers
 * could reduce both the number of dense subgraphs with which the dsg is evaluated via supersede(), and vice versa.
 */
bool
DenseSubGraphsMaximalSet::insert(const DenseSubGraph& dsg) {

    // Do nothing if dsg is already superseded
    for (const_iterator it = begin(); it != end(); ++it) {
        if (supersede(*it, dsg))    // It take into account when both are equals too
            return false;
    }

    // Ok, dsg is maximal. Before adding it, we need to remove any other dsg superseded by it
    for (iterator it = begin(); it != end();  ) {
        assert(dsg != *it);

        if (supersede(dsg,*it))
            it = dsgs.erase(it);    // erase() is the one who can give us the proper iterator to continue
        else
            ++it;
    }

    dsgs.push_back(dsg);
    return true;
}


bool
DenseSubGraphsMaximalSet::supersede(const DenseSubGraph& candidate, const DenseSubGraph& dsg) const {
    if (onlyCentersMaximality)
        return algorithms::set_includes(candidate.getCenters(), dsg.getCenters());
    else
        return candidate.includes(dsg);
}


void
DenseSubGraphsMaximalSet::dump(const std::string& fileName, bool append, bool includeDescriptions) const {
    assert(!fileName.empty());

    std::ofstream outfile(fileName.c_str(), (append ? std::ios::app : std::ios::out));
    if (!outfile) {
        throw std::runtime_error("DenseSubGraphsMaximalSet::dump(): can not open output dump file");
    }

    if (append)
        outfile << '\n';

    // Dump one dense subgraph by line
    for (const_iterator it = begin(); it != end(); ++it) {
        if (onlyCentersMaximality) {
            outfile << it->getCenters();
        } else {
            outfile << *it;
            if (includeDescriptions && !it->clique() && !it->generic()) {
                // Any non-digit character at the start is Ok to mark the end of the already streamed dense subgraph
                outfile << "  ## " << it->description();
            }
        }

        outfile << '\n';
    }
}


}   // namespace odsg
