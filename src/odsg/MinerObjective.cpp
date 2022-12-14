#include "MinerObjective.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro

#include "utils/algorithms.hpp"
#include "DenseSubGraph.hpp"

//==============================================================================
#include "WedgeMap.hpp"
//==============================================================================

namespace odsg {


//// AsCliqueMinerObjective ///////////////////////////////////////////////////////////////////////////////////////////

/*
 * Note that the 'current' parameter isn't used in non-debug builds, so the ((__unused__)) attribute prevents
 * the respective warning from gcc (or clang) for that cases:
 *   https://gcc.gnu.org/onlinedocs/gcc/Common-Variable-Attributes.html
 */
bool
AsCliqueMinerObjective::better(const DenseSubGraph& current __attribute__ ((__unused__)),
                               const DenseSubGraph& candidate) const {

    assert(current.getCenters().size() < candidate.getCenters().size());
    return algorithms::set_includes(candidate.getSources(), candidate.getCenters());
}


bool
AsCliqueMinerObjective::best(const DenseSubGraph& candidate) const {
    return candidate.getCenters().size() == candidate.getSources().size();
}


//// LegacyMinerObjective /////////////////////////////////////////////////////////////////////////////////////////////

bool
LegacyMinerObjective::better(const DenseSubGraph& current,
                             const DenseSubGraph& candidate) const {
    return current.arcsCount() < candidate.arcsCount();
}


//// MaxIntersectionMinerObjective ////////////////////////////////////////////////////////////////////////////////////

bool
MaxIntersectionMinerObjective::better(const DenseSubGraph& current,
                                      const DenseSubGraph& candidate) const {
    return algorithms::set_intersection_count(current.getCenters(), current.getSources()) <
           algorithms::set_intersection_count(candidate.getCenters(), candidate.getSources());
}

//_DensityMinerObjective________________________________________________________

SimpleEdgeDensity::SimpleEdgeDensity(const WedgeMap* wedgemap_ptr)
: threshold(0.0f), wedgemap(wedgemap_ptr) { }

bool
SimpleEdgeDensity::better(const DenseSubGraph& current,
                          const DenseSubGraph& candidate) const {

    return wedgemap->get_simple_edge_density( current ) <
           wedgemap->get_simple_edge_density( candidate );
}

FullEdgeDensity::FullEdgeDensity(const WedgeMap* wedgemap_ptr)
: threshold(0.0f), wedgemap(wedgemap_ptr) { }

bool
FullEdgeDensity::better(const DenseSubGraph& current,
                        const DenseSubGraph& candidate) const {

    return wedgemap->get_full_edge_density( current ) <
           wedgemap->get_full_edge_density( candidate );
}


SimpleDegreeDensity::SimpleDegreeDensity(const WedgeMap* wedgemap_ptr)
: wedgemap(wedgemap_ptr) { }

bool
SimpleDegreeDensity::better(const DenseSubGraph& current,
                            const DenseSubGraph& candidate) const {

   float a = wedgemap->get_simple_degree_density( current );
   float b = wedgemap->get_simple_degree_density( candidate );
   if( a < b )
      return true;//wedgemap->get_simple_edge_density( candidate ) > 0.43f;

   return false;

}

FullDegreeDensity::FullDegreeDensity(const WedgeMap* wedgemap_ptr)
: wedgemap(wedgemap_ptr) { }


bool
FullDegreeDensity::better(const DenseSubGraph& current,
                            const DenseSubGraph& candidate) const {
   return wedgemap->get_full_degree_density( current ) <
          wedgemap->get_full_degree_density( candidate );
}


DegreeAndEdgeDensity::DegreeAndEdgeDensity(const WedgeMap* wedgemap_ptr)
: wedgemap(wedgemap_ptr) { }

bool
DegreeAndEdgeDensity::better(const DenseSubGraph& current,
                             const DenseSubGraph& candidate) const {

   float   current_value;
   float candidate_value;
   current_value   = wedgemap->get_u_simple_degree_density( current );
   candidate_value = wedgemap->get_u_simple_degree_density( candidate );
   return current_value < candidate_value;

}

//_____________________________________________________End_DensityMinerObjective
}   // namespace odsg
