#include<iostream>
#include "WedgeMap.hpp"

namespace odsg {

//_WedgeMap_Methods_____________________________________________________________

void
WedgeMap::print_map() const {
   std::map<Edge, float>::const_iterator it;
   for( it = edge_map.begin(); it != edge_map.end(); ++it ) {
      std::cout<<((Edge)(it->first)).first <<"\t"
               <<((Edge)(it->first)).second<<"\t"
               <<(float)(it->second)       <<"\n";
   }   
}
/*
float 
WedgeMap::simple_average_weight(const DenseSubGraph& dsg) {
   const VertexSet c = dsg.getCenters();
   const VertexSet s = dsg.getSources();
   return c.size() + s.size() < 5 ? 0.0f : weight_sum( c, s, true );
}
*/
float 
WedgeMap::get_simple_edge_density(const DenseSubGraph& dsg) const {
   const VertexSet c = dsg.getCenters();
   const VertexSet s = dsg.getSources();
   return c.size() + s.size() < 5 ? 0.0f : weight_sum( c, s, true );
}

float 
WedgeMap::get_full_edge_density(const DenseSubGraph& dsg) const {
   const VertexSet c = dsg.getCenters();
   const VertexSet s = dsg.getSources();
   return (c.size() + s.size()) < 5 ? 0.0f : full_weight_sum( c, s, true );
}

float
WedgeMap::get_simple_degree_density( const DenseSubGraph& dsg ) const {
   const VertexSet centers = dsg.getCenters();
   const VertexSet sources = dsg.getSources();
   if(centers.size() + sources.size() < 5)return 0.0f;
   VertexSet set_union;
   set_union.insert( centers.begin(), centers.end() );
   set_union.insert( sources.begin(), sources.end() );
   return weight_sum(centers, sources) / set_union.size();
}

float
WedgeMap::get_full_degree_density( const DenseSubGraph& dsg ) const {
   const VertexSet c = dsg.getCenters();
   const VertexSet s = dsg.getSources();
   return (c.size() + s.size()) < 5 ? 0.0f : full_weight_sum(c, s);
}

float
WedgeMap::get_u_simple_degree_density( const DenseSubGraph& dsg ) const {
   const VertexSet centers = dsg.getCenters();
   const VertexSet sources = dsg.getSources();
   const int size_c = centers.size();
   const int size_s = sources.size();
   if( size_c + size_s < 5)return 0.0f;
   VertexSet set_union;
   set_union.insert( centers.begin(), centers.end() );
   set_union.insert( sources.begin(), sources.end() );
   return ((float)(size_c*size_s)) / ((float)set_union.size());
}

float
WedgeMap::get_u_full_degree_density( const DenseSubGraph& dsg ) const {
   const VertexSet centers = dsg.getCenters();
   const VertexSet sources = dsg.getSources();
   if( centers.size() + sources.size() < 5)return 0.0f;
   VertexSet set_union;
   set_union.insert( centers.begin(), centers.end() );
   set_union.insert( sources.begin(), sources.end() );  
   return edge_count(set_union) / set_union.size();
}

//_________________________________________________________End_WedgeMap_Methods_


//_UndirectedWedgeMap_Methods___________________________________________________

void
UndirectedWedgeMap::add_edge( Vertex v1, Vertex v2, float value ) {
   if(v1 < v2) edge_map.insert(std::make_pair( Edge(v1, v2), value ));
   else        edge_map.insert(std::make_pair( Edge(v2, v1), value ));
}

float
UndirectedWedgeMap::weight_sum (const VertexSet& centers,
                                const VertexSet& sources, const bool average)
                                const {
   std::set<Vertex>::const_iterator it_c;
   std::set<Vertex>::const_iterator it_s;
   std::map<Edge, float>::const_iterator it;
   int edge_counter =    0;
   float sum        = 0.0f;
   Vertex            v1,v2;
   std::set<Edge, EdgeComparer> edge_set;
   for( it_c = centers.begin(); it_c != centers.end(); ++it_c )
      for( it_s = sources.begin(); it_s != sources.end(); ++it_s ) {
         if( (v1 = *it_c) == (v2 = *it_s) ) continue;
         if( v1 < v2 ) edge_set.insert(Edge(v1, v2));
         else          edge_set.insert(Edge(v2, v1));
      }
   std::set<Edge, EdgeComparer>::const_iterator it_set;
   
   for( it_set = edge_set.begin(); it_set != edge_set.end(); ++it_set )
      if( (it = edge_map.find(*it_set)) != edge_map.end() ) {
         sum += (float)it->second;
         edge_counter++;
      }
   return average ? sum/edge_counter : sum;
}

float
UndirectedWedgeMap::full_weight_sum (const VertexSet& centers,
                                     const VertexSet& sources,
                                     const bool average) const {
   VertexSet set_union;
   set_union.insert( centers.begin(), centers.end() );
   set_union.insert( sources.begin(), sources.end() );
   const unsigned int union_size = set_union.size();
   if( union_size < 2 )return 0.0f;
   const unsigned int max_index  =   union_size - 1;
   Vertex *v = new Vertex[union_size];
   std::copy(set_union.begin(), set_union.end(), v);
   std::map<Edge, float>::const_iterator it;
   int edge_counter =    0;
   float sum        = 0.0f;
   for( int i = 0; i < max_index; i++ )
      for( int j = max_index; j > i; j-- )
         if( (it = edge_map.find(Edge(v[i],v[j]))) != edge_map.end() ) {
            sum += (float)it->second;
            edge_counter++;
         }
   delete []v;
   return average ? sum/edge_counter : sum/union_size;
}

float
UndirectedWedgeMap::edge_count( const VertexSet& vertex_set ) const {
   const unsigned int vertex_set_size = vertex_set.size();
   if( vertex_set_size < 2 )return 0.0f;
   const unsigned int max_index  =   vertex_set_size - 1;
   Vertex *v = new Vertex[vertex_set_size];
   std::copy(vertex_set.begin(), vertex_set.end(), v);
   std::map<Edge, float>::const_iterator it;
   int sum = 0;
   for( int i = 0; i < max_index; i++ )
      for( int j = max_index; j > i; j-- )
         if( (it = edge_map.find(Edge(v[i],v[j]))) != edge_map.end() ) sum++;
   delete []v;
   return (float)sum;
}

//_______________________________________________End_UndirectedWedgeMap_Methods_

//_DirectedWedgeMap_Methods_____________________________________________________
/*
void
DirectedWedgeMap::add_edge( Vertex v1, Vertex v2, float value ) {
   edge_map.insert (std::make_pair(Edge(v1, v2), value));
}


float
DirectedWedgeMap::get_weight( Vertex v1, Vertex v2 ) const {
   std::map<Edge, float>::const_iterator it;
   it = edge_map.find( Edge(v1,v2) );
   return it != edge_map.end() ? (float)(it->second) : 0.0f;
}


std::pair<float,int>
DirectedWedgeMap::weight_sum(const DenseSubGraph& dsg) const {
   
   const VertexSet centers = dsg.getCenters();
   const VertexSet sources = dsg.getSources();
   
   std::set<Vertex>::const_iterator it_c;
   std::set<Vertex>::const_iterator it_s;

   std::map<Edge, float>::const_iterator it;
   
   float sum = 0.0f;
   Vertex v1, v2;
   
   std::set<Edge, EdgeComparer> edge_set;
   
   for( it_c = centers.begin(); it_c != centers.end(); ++it_c )
      for( it_s = sources.begin(); it_s != sources.end(); ++it_s ) {
         if( (v1 = *it_c) == (v2 = *it_s) ) continue;
         edge_set.insert(Edge(v1, v2));
      }
   std::set<Edge, EdgeComparer>::const_iterator it_set;
   
   for( it_set = edge_set.begin(); it_set != edge_set.end(); ++it_set )
      if( (it = edge_map.find(*it_set)) != edge_map.end() )
         sum += (float)it->second;
   
   return sum;
   
   
  // return 0.0f;
}
*/
//_________________________________________________End_DirectedWedgeMap_Methods_
   
} // namespace   
