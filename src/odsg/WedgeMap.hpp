#ifndef SRC_WEDGEMAP_HPP_INCLUDED
#define SRC_WEDGEMAP_HPP_INCLUDED

#include <map>
#include <functional>   // std::binary_function
#include <utility>      // std::pair, std::make_pair

#include   "WGraphTypes.hpp"
#include "DenseSubGraph.hpp"

namespace odsg {

class WedgeMap {

public:

   WedgeMap() {}

   virtual ~WedgeMap() { }
   virtual void   add_edge( Vertex v1, Vertex v2, float value ) = 0;
   virtual int    get_type()                              const = 0;
   virtual void  print_map()                                  const;

   //float simple_average_weight      ( const DenseSubGraph& dsg ) const;
   //float get_weight_edge_average( const DenseSubGraph& dsg ) const;
   
   float get_simple_edge_density    ( const DenseSubGraph& dsg ) const;
   float get_full_edge_density      ( const DenseSubGraph& dsg ) const;
   float get_simple_degree_density  ( const DenseSubGraph& dsg ) const;
   float get_full_degree_density    ( const DenseSubGraph& dsg ) const;
   
   float get_u_simple_degree_density( const DenseSubGraph& dsg ) const;
   float get_u_full_degree_density  ( const DenseSubGraph& dsg ) const;
   
protected:

   struct Edge {
      const Vertex  first;
      const Vertex second;
      Edge( Vertex v1, Vertex v2 ) : first(v1), second(v2) { }
   };

   struct EdgeComparer: public std::binary_function<Edge, Edge, bool> {
      inline bool operator()(Edge e1, Edge e2) const {
         Vertex v1 = e1.first;
         Vertex v2 = e2.first;   
         return v1 != v2 ? v1 < v2 : e1.second < e2.second;
      }
   };
   
   std::map<Edge, float, EdgeComparer>  edge_map;
   
   virtual float weight_sum
                 (const VertexSet& centers, 
                  const VertexSet& sources,
                  const bool average = false) const = 0;

   virtual float full_weight_sum
                 (const VertexSet& centers,
                  const VertexSet& sources,
                  const bool average = false) const = 0;
                  
   virtual float edge_count(const VertexSet& vertex_set) const = 0;
};

class UndirectedWedgeMap : public WedgeMap {

public:
   
   UndirectedWedgeMap() { }

   int get_type() const { return UNDIRECTED_WITH_SYMETRIC_WEIGHT; }
   
   void    add_edge( Vertex v1, Vertex v2, float value );
   float get_weight( Vertex v1, Vertex v2 )        const;

private:

   float weight_sum(const VertexSet& centers, const VertexSet& sources,
                    const bool average = false) const;

   float full_weight_sum(const VertexSet& centers, const VertexSet& sources,
                         const bool average = false) const;
   
   float edge_count( const VertexSet& vertex_set ) const;
};

/*
class DirectedWedgeMap : public WedgeMap {

public:   
   
   DirectedWedgeMap(){ }

   int     get_type() const { return DIRECTED_WEIGHTED; }
   
   void    add_edge( Vertex v1, Vertex v2, float value );
   float get_weight( Vertex v1, Vertex v2 )        const;

private:   
   float weight_sum(const DenseSubGraph& dsg) const;
};
*/

} // namespace

#endif  // SRC_WEDGEMAP_HPP_INCLUDED
