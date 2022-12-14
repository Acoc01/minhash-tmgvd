#ifndef BIO_READ_FILE_HPP_INCLUDED
#define BIO_READ_FILE_HPP_INCLUDED

#include <string>

//==============================================================================
//#include <odsg/Graph.hpp>
#include <odsg/WGraph.hpp>
//==============================================================================

#include "typedefs.hpp"


namespace bio_odsg {

/*
 * Read a text file containing a custom mapping from proteins from a dataset to internal integer identifiers.
 *
 * This file must contain one pair by line, a protein specified by name and an integer id, both separated by spaces.
 * The function will throw if the file contains an id previously seen for other protein, or a protein previously
 * seen with other id.
 *
 * The mapping don't need to cover all the proteins in the dataset, but all unknown proteins will get an id greater
 * than the greatest in this custom mapping. Also, the order of the lines in the mapping file is irrelevant: later,
 * the graph for the dataset is built with the adjacency lists sorted according to the id.
 */
ProteinsMap
readDatasetMappingFromFile(const std::string& fileName);


/*
 * Read a text file containing a network of protein-protein interactions (PPI): a dataset.
 *
 * This file must contain one interaction by line, both proteins specified by name and separated by spaces.
 * A reliability score can be found accompanying each pair in most cases.
 *
 * Returns the graph and the mapping updated with all the new proteins seen in the dataset.
 */
std::vector<odsg::WGraph>
readDatasetFromFileWW(const std::string& fileName,
                      bool weighted,
                      ProteinsMap& mapping);          // Out-parameter


}       // namespace bio_odsg
#endif  // BIO_READ_FILE_HPP_INCLUDED
