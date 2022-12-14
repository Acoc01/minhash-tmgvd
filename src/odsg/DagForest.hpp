#ifndef SRC_DAG_FOREST_HPP_INCLUDED
#define SRC_DAG_FOREST_HPP_INCLUDED

#include <cstddef>      // std::size_t
#include <string>
#include <vector>
#include <iosfwd>

#include "Dag.hpp"              // All of our container-like classes include the definition of the contained element

namespace odsg {


class Graph;
/*
 * The main motivations to have a collection of dags as an class (versus passing std::vector<Dag>, by example, all
 * around the place) were:
 *   1) To have a single place where to group different strategies about how dags can be built from a Graph object:
 *      sometimes we might want a single dag covering the whole graph, if processing memory is not a problem, or
 *      to generate many dags, according to different partitioning/clustering schemes, each one implemented by
 *      a GraphPartitioner object.
 *   2) As copying dag objects is discouraged (disabled, in fact), the dags are saved and passed via pointers, so
 *      the respective RAII class was desired.
 *
 * Note that 'partitioning' and 'clustering' terms are treated as synonyms, though 'partition' and 'cluster' aren't.
 * Sorry about that.
 */
class DagForest {
public:
    // Types
    typedef std::vector<const Dag*>::const_iterator const_iterator;

    // Constructors
    explicit DagForest(const Graph&,
                       int clusteringScheme=0,
                       unsigned int minClusterSize=1,       // With 'size' we refers to the number of arcs
                       bool sortClusterByFrequency=false);  // It can throw an exception
    ~DagForest();

    // No public mutators: a dag forest isn't altered outside of the constructor & destructor

    // Iterators
    const_iterator begin() const { return forest.begin(); }
    const_iterator end() const { return forest.end(); }

    // Inspectors
    bool empty() const { return forest.empty(); }
    std::size_t size() const { return forest.size(); }

    friend std::ostream& operator<<(std::ostream&, const DagForest&);   // Short summary, made to fit in one line
    void print(std::ostream&, bool onlySummaries=false) const;

    void dump(const std::string& fileName, bool onlySummaries=false) const;   // It can throw an exception

private:
    std::vector<const Dag*> forest;

    // General helpers
    typedef std::vector<const Dag*>::iterator iterator;

    iterator begin() { return forest.begin(); }
    iterator end() { return forest.end(); }


    // The next two are declared and deliberately NOT implemented, to prevent copying objects of this class
    DagForest(const DagForest&);
    DagForest& operator=(const DagForest&);
};


}       // namespace odsg
#endif  // SRC_DAG_FOREST_HPP_INCLUDED
