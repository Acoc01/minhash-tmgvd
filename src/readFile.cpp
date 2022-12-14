#include "readFile.hpp"

#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <stdexcept>
#include <cstdlib>      // std::atof
#include <map>
#include <fstream>
#include <sstream>

#include <odsg/utils/strings.hpp>

namespace bio_odsg {


namespace {     // Put here general, global definitions limited to this file

    ProteinId
    maxProteinIdIn(const ProteinsMap& mapping) {
        ProteinId maxId = 0;
        for (ProteinsMap::const_iterator it = mapping.begin(); it != mapping.end(); ++it) {
            if (it->second > maxId)
                maxId = it->second;
        }
        return maxId;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ProteinsMap
readDatasetMappingFromFile(const std::string& fileName) {

    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readDatasetMappingFromFile(): can not open input file");
    }

    ProteinsMap mapping;

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')       // TODO: drop lines with empty spaces too
            continue;

        std::istringstream iss(line);

        ProteinName protein;
        if (!(iss >> protein)) {
            throw std::runtime_error("bad line format");
        }
        ProteinId id;
        if (!(iss >> id)) {
            throw std::runtime_error("bad line format");
        }

        // Basic consistency checks
        if (mapping.find(protein) != mapping.end()) {
            throw std::invalid_argument(std::string("bad mapping: found duplicated protein: ") + protein);
        }
        if (isProteinIdInProteinsMap(mapping, id)) {
            throw std::invalid_argument(std::string("bad mapping: found duplicated id: ") + odsg::strings::to_str(id));
        }

        mapping.insert(std::make_pair(protein, id));
    }
    return mapping;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<odsg::WGraph>
readDatasetFromFileWW(const std::string& fileName, bool weighted, ProteinsMap& mapping) {

    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readDatasetFromFileWW(): can not open input file");
    }

    std::map<ProteinId, std::set<ProteinId> > dataset;
    odsg::UndirectedWedgeMap ppi_dataset;
    std::vector<odsg::WGraph> Clusters;

    unsigned int nextProteinId = maxProteinIdIn(mapping) + 1;

    std::string line;
    std::string interaction_value;
    odsg::Vertex left_vertex;
    odsg::Vertex right_vertex;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#'){
            Clusters.push_back(odsg::WGraph(dataset,ppi_dataset));
            dataset.clear();
            odsg::UndirectedWedgeMap ppi_nuevo;
            ppi_dataset = ppi_nuevo;
            continue;
        }

        std::istringstream iss(line);

        ProteinName lprotein;
        if (!(iss >> lprotein)) {
            throw std::runtime_error("bad line format");
        }
        ProteinName rprotein;
        if (!(iss >> rprotein)) {
            throw std::runtime_error("bad line format");
        }

        // Existent weights can be ignored if choosed, and missing weights can have a default value
        if (!weighted || !(iss >> interaction_value)) {
            interaction_value = "1.0";
        }

        // If the proteins are unknown, save them in the mapping
        if (mapping.find(lprotein) == mapping.end()) {
            mapping.insert(std::make_pair(lprotein, nextProteinId));
            nextProteinId++;
        }
        if (mapping.find(rprotein) == mapping.end()) {
            mapping.insert(std::make_pair(rprotein, nextProteinId));
            nextProteinId++;
        }

        // Save the interaction in the adjacency lists of both proteins
        left_vertex  = mapping[lprotein];
        right_vertex = mapping[rprotein];

        dataset[left_vertex].insert(right_vertex);
        dataset[right_vertex].insert(left_vertex);

        ppi_dataset.add_edge(left_vertex, right_vertex, std::atof(interaction_value.c_str()));
    }
    return Clusters;
    //return odsg::WGraph(dataset, ppi_dataset);
}

}   // namespace bio_odsg
