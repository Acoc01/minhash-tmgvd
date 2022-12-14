#ifndef SRC_WGRAPH_HPP_INCLUDED
#define SRC_WGRAPH_HPP_INCLUDED

#include    "Graph.hpp"
#include "WedgeMap.hpp"

namespace odsg {

class WGraph : public Graph {

public:

   WGraph();

   WGraph(const WGraph& WG);

   WGraph(const std::map<Vertex, VertexSet>& dataset,
          const UndirectedWedgeMap&            emap);

   //WGraph(const std::map<Vertex, VertexSet>& dataset,
   //       const DirectedWedgeMap&              emap);

  ~WGraph();

   WGraph& operator = (const WGraph& WG);

   int  get_type()                const;

   void print_edge_map()          const;
   //void rebuildForMiningExceptSorting();

   WedgeMap* get_edge_map() const;

private:

   WedgeMap* edge_map;
   void clone_edge_map_ptr( const WedgeMap *ptr );
};

}// namespace odsg

#endif  // SRC_WGRAPH_HPP_INCLUDED
