#ifndef SRC_DENSE_SUB_GRAPH_HPP_INCLUDED
#define SRC_DENSE_SUB_GRAPH_HPP_INCLUDED

#include <iosfwd>

#include "VertexSet.hpp"

namespace odsg {


class DenseSubGraph {
public:
//==============================================================================
    DenseSubGraph(): sources(), centers(), density_value(0.0f) {}
//==============================================================================
    DenseSubGraph(const VertexSet&, const VertexSet&);
    DenseSubGraph(const VertexSet&, Vertex);

    // Operators
    bool operator==(const DenseSubGraph&) const;
    bool operator!=(const DenseSubGraph&) const;

    // Mutators
    DenseSubGraph& merge(const DenseSubGraph&);
    void swap(DenseSubGraph&);      // WARN: use it directly only; not relies in lookup for general swap

//==============================================================================
    void set_density_value( float value );
//==============================================================================

    // Inspectors
    const VertexSet& getSources() const { return sources; }     // Const-references to internals
    const VertexSet& getCenters() const { return centers; }     //

    bool includes(const DenseSubGraph&) const;

    bool empty() const { return sources.empty() && centers.empty(); }
    unsigned long arcsCount() const;
    unsigned long nodesCount() const;

    bool clique() const;
    bool asClique() const;      // Aka its centers set can be treated as a clique
    bool biClique() const;
    bool generic() const { return !clique() && !asClique() && !biClique(); }

    std::string description() const;    // An human-readable, short description
    
//==============================================================================
    float get_density_value() const;
//==============================================================================    
    
    friend std::istream& operator>>(std::istream&, DenseSubGraph&);
    friend std::ostream& operator<<(std::ostream&, const DenseSubGraph&);

private:
    VertexSet sources;
    VertexSet centers;
//==============================================================================    
    float density_value;
//==============================================================================    
};


}       // namespace odsg
#endif  // SRC_DENSE_SUB_GRAPH_HPP_INCLUDED
