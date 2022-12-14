#include<cassert> // Support run-time assertions. They can be disabled defining the NDEBUG macro

#include           "WGraph.hpp"
#include      "WGraphTypes.hpp"
#include "utils/algorithms.hpp"

namespace odsg {

void
WGraph::clone_edge_map_ptr( const WedgeMap *ptr ) {
   if( !ptr ) return;
   switch( ptr->get_type() ) {
      case UNDIRECTED_WITH_SYMETRIC_WEIGHT:
           edge_map = new UndirectedWedgeMap
                      (*(dynamic_cast<const UndirectedWedgeMap*>(ptr)));
           break;
      /*     
      case DIRECTED_WEIGHTED:
           edge_map = new DirectedWedgeMap
                      (*(dynamic_cast<const DirectedWedgeMap*>(ptr)));
           break;
      */     
      //add here more types
   }
}

WGraph::WGraph(){ edge_map = 0; }

WGraph::WGraph(const WGraph& WG) : Graph(WG) {
	clone_edge_map_ptr( WG.edge_map );   
}

WGraph::WGraph( const std::map<Vertex, VertexSet>& dataset,
                const UndirectedWedgeMap& emap )
              : Graph(dataset),
                edge_map(new UndirectedWedgeMap(emap)) { }

/*
WGraph::WGraph( const std::map<Vertex, VertexSet>& dataset,
                const DirectedWedgeMap& emap)
              : Graph(dataset),
                edge_map(new DirectedWedgeMap(emap)) { }
*/

WGraph::~WGraph(){ delete edge_map; }

WGraph&
WGraph::operator = (const WGraph& WG) {
   static_cast<Graph&>(*this) = WG;
   clone_edge_map_ptr( WG.edge_map );
   return *this;
}

int
WGraph::get_type() const {
   if( edge_map )return edge_map->get_type();
   return EMPTY_WEIGHTED_GRAPH;
}

void
WGraph::print_edge_map() const {
   if( edge_map ) edge_map->print_map();
}

WedgeMap*
WGraph::get_edge_map() const {
   return edge_map;
}

/*
void
WGraph::rebuildForMiningExceptSorting() {
   if (mineability >= 1) return;
   Vertex vertex;
   for (iterator it = begin(); it != end(); ) {
      vertex = it->first;
      
      //Adding self loops to map
      //edge_map->add_edge( vertex, vertex, 1.0f );
      
      AdjacencyList& outlinks = it->second;

     if (outlinks.size() >= 1 && !algorithms::is_found(outlinks, vertex)) {
         if (sortedByVertex) {
             outlinks.insert(std::upper_bound(outlinks.begin(), outlinks.end(), vertex), vertex);
             assert(algorithms::is_sorted(outlinks));
         }
         else outlinks.push_back(vertex);
     }
     if (outlinks.size() <= 1) {
         iterator tmp = it;
         ++it;
         innerGraph.erase(tmp);
     }
     else ++it;
   }
   mineability = 1;
}
*/
}// namespace
