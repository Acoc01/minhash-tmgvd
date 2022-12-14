#ifndef SRC_VERTEXSET_HPP_INCLUDED
#define SRC_VERTEXSET_HPP_INCLUDED

#include <set>
#include <ostream>
#include <istream>

#include "Vertex.hpp"

namespace odsg {


typedef std::set<Vertex> VertexSet;


inline std::ostream&
operator<<(std::ostream& os, const VertexSet& vxset) {
    for (VertexSet::const_iterator vxit = vxset.begin(); vxit != vxset.end(); ++vxit) {
        if (vxit != vxset.begin())
            os << ' ';
        os << *vxit;
    }
    return os;
}

inline std::istream&
operator>>(std::istream& is, VertexSet& vxset) {
    Vertex vx;
    while (is >> vx) {
        vxset.insert(vx);
    }
    is.clear();     // Getting even an empty VertexSet is Ok

    return is;
}


}       // namespace odsg
#endif  // SRC_VERTEXSET_HPP_INCLUDED
