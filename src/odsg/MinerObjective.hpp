#ifndef SRC_MINER_OBJECTIVE_HPP_INCLUDED
#define SRC_MINER_OBJECTIVE_HPP_INCLUDED

namespace odsg {


class DenseSubGraph;
//==============================================================================
class      WedgeMap;
//==============================================================================
/*
 * MinerObjective objects let to customize one aspect of the mining process done by DenseSubGraphsMiner objects:
 * Which factor to favor to try to mine the more and 'bigger' (whatever how this word is interpreted) dense subgraphs.

 * It's an abstract class; derived implementations define different strategies.
 */
class MinerObjective {
public:
    virtual ~MinerObjective() {}

    /*
     * Decide if a newly built dense subgraph (a 'candidate') is preferable over the currently kept dense subgraph,
     * while traveling though a mining path.
     */
    virtual bool better(const DenseSubGraph&, const DenseSubGraph&) const = 0;

    /*
     * For some cases, it's possible to cut the traveling though a mining path because can be ensured that no better
     * dense subgraphs will be found. It helps to reduce the processing time.
     */
    virtual bool best(const DenseSubGraph&) const { return false; }
//==============================================================================
    virtual bool has_best() const { return false; }
//==============================================================================
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Maximize the size of the centers set, while keeping it as a subset of the sources set, favoring (but not limited
 * to) the mining of cliques.
 *
 * Also, it's the unique currently possible to use in conjuntion with cliques mining.
 */
class AsCliqueMinerObjective : public MinerObjective {
public:
    /*virtual*/ bool better(const DenseSubGraph&, const DenseSubGraph&) const;
    /*virtual*/ bool best(const DenseSubGraph&) const;
    /*virtual*/ bool has_best() const { return true; }

};


/*
 * Maximize the number of arcs of the dense subgraphs mined.
 *
 * It's the same idea as the legacy implementation by chernand.
 */
class LegacyMinerObjective : public MinerObjective {
public:
    /*virtual*/ bool better(const DenseSubGraph&, const DenseSubGraph&) const;
};


/*
 * Maximize the number of shared elements between sources and centers.
 *
 * With little examples, it provides very similar results with AsCliqueMinerObjective, but with some huge social/web
 * graphs it seems to provide larger collections of dense subgraphs.
 */
class MaxIntersectionMinerObjective : public MinerObjective {
public:
    /*virtual*/ bool better(const DenseSubGraph&, const DenseSubGraph&) const;
};


//==============================================================================

class SimpleEdgeDensity : public MinerObjective {

public:
   SimpleEdgeDensity( const WedgeMap* wedgemap_ptr );
   bool better(const DenseSubGraph&, const DenseSubGraph&) const;

private:
   const float threshold;
   const WedgeMap* wedgemap;
};

class FullEdgeDensity : public MinerObjective {

public:
   FullEdgeDensity( const WedgeMap* wedgemap_ptr );
   bool better(const DenseSubGraph&, const DenseSubGraph&) const;

private:
   const float threshold;
   const WedgeMap* wedgemap;
};


class SimpleDegreeDensity : public MinerObjective {
public:
   SimpleDegreeDensity( const WedgeMap* wedgemap_ptr );
   bool better(const DenseSubGraph&, const DenseSubGraph&) const;

private:
   const WedgeMap* wedgemap;
};

class FullDegreeDensity : public MinerObjective {
public:
   FullDegreeDensity( const WedgeMap* wedgemap_ptr );
   bool better(const DenseSubGraph&, const DenseSubGraph&) const;

private:
   const WedgeMap* wedgemap;
};


class DegreeAndEdgeDensity : public MinerObjective {
public:
   DegreeAndEdgeDensity( const WedgeMap* wedgemap_ptr );
   bool better(const DenseSubGraph&, const DenseSubGraph&) const;

private:
   const WedgeMap* wedgemap;
};

//==============================================================================

}       // namespace odsg
#endif  // SRC_MINER_OBJECTIVE_HPP_INCLUDED
