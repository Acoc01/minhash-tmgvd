#ifndef BIO_TYPEDEFS_HPP_INCLUDED
#define BIO_TYPEDEFS_HPP_INCLUDED

#include <string>
#include <map>
#include <set>

namespace bio_odsg {


typedef std::string ProteinName;
typedef unsigned int ProteinId;

typedef std::map<ProteinName, ProteinId> ProteinsMap;   // TODO: Move to boost::bimap
typedef std::set<ProteinId> Complex;                    // TODO? Move to vector<>



/*
 * The next two, dirty, auxiliary utils will be removed after moving ProteinsMap to boost::bimap
 */
inline bool
isProteinIdInProteinsMap(const ProteinsMap& mapping, ProteinId id) {
    for (ProteinsMap::const_iterator it = mapping.begin(); it != mapping.end(); ++it) {
        if (it->second == id)
            return true;
    }
    return false;
}
inline ProteinName
getProteinName(const ProteinsMap& mapping, ProteinId id) {
    for (ProteinsMap::const_iterator it = mapping.begin(); it != mapping.end(); ++it) {
        if (it->second == id)
            return it->first;
    }
    return "";
}


}       // namespace bio_odsg
#endif  // BIO_TYPEDEFS_HPP_INCLUDED
