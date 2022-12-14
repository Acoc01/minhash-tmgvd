#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <string>
#include <exception>
#include <stdexcept>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <utility>

#include <tclap/CmdLine.h>

#include <odsg/utils/algorithms.hpp>


struct CmdLineArgs {    // The definition of processCmdLine() constains descriptions for each option
    // Input files
    std::string referenceComplexesFileName;

    // Other options
    unsigned int minComplexSize;
    double minOverlapScore;
};
CmdLineArgs processCmdLine(int argc, char* argv[]);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * The functions of the bio_odsg namespace works with proteins being treated as numeric ids, something not required
 * here, so, instead, we re-implement here some stuff from there.
 */
typedef std::string ProteinName;
typedef unsigned int ComplexId;
typedef std::set<ProteinName> ComplexData;


std::map<ComplexId, ComplexData>
readComplexesFromFile(const std::string& fileName, unsigned int minSize) {

    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readComplexesFromFile(): can not open input file");
    }

    std::map<ComplexId, ComplexData> complexesData;
    ComplexId currentComplexId = 0;

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')       // TODO: drop lines with empty spaces too
            continue;

        std::istringstream iss(line);

        ComplexData complexData;
        for (ProteinName protein; iss >> protein; ) {
            complexData.insert(protein);
        }

        if (complexData.size() < minSize)
            continue;

        complexesData[currentComplexId] = complexData;
        currentComplexId++;
    }

    return complexesData;
}


inline double
overlapScore(const ComplexData& complexData, const ComplexData& otherComplexData) {
    unsigned int overlap = odsg::algorithms::set_intersection_count(complexData, otherComplexData);

    return double(overlap * overlap) / (complexData.size() * otherComplexData.size());
}


inline std::ostream&
operator<<(std::ostream& os, const ComplexData& complexData) {
    for (ComplexData::const_iterator it = complexData.begin(); it != complexData.end(); ++it) {
        if (it != complexData.begin())
            os << " ";
        os << *it;
    }
    return os;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(int argc, char* argv[]) {

    CmdLineArgs args;
    try {
        args = processCmdLine(argc, argv);
    } catch (TCLAP::ArgException& e) {
        std::cerr << "error: " << e.argId() << '\n'
                  << "       " << e.error() << std::endl;
        return 1;
    }

    std::cout << "#Considering protein complexes with a minimum size of " << args.minComplexSize << std::endl;


    // Read the references from file
    std::cout << "#Reading reference complexes file... " << std::flush;
    std::map<ComplexId, ComplexData> references;
    try {
        references = readComplexesFromFile(args.referenceComplexesFileName, args.minComplexSize);
    } catch (std::exception& e) {
        std::cerr << "ERROR\n"
                  << e.what() << std::endl;
        return 1;
    }
    std::cout << "#OK\n"
              << "#\t" << references.size() << " complexes found in '" << args.referenceComplexesFileName << "'\n";


    // Save the OS between all the pair of reference complexes over the threhold
    std::map<std::pair<ComplexId, ComplexId>, double> ooss;

    for (std::map<ComplexId, ComplexData>::const_iterator it = references.begin(); it != references.end(); ++it) {
        for (std::map<ComplexId, ComplexData>::const_iterator it2 = it; it2 != references.end(); ++it2) {

            if (it == it2)
                continue;

            const ComplexId& complexId1 = it->first;
            const ComplexId& complexId2 = it2->first;
            const ComplexData& complexData1 = it->second;
            const ComplexData& complexData2 = it2->second;

            double os = overlapScore(complexData1, complexData2);

            if (os >= args.minOverlapScore)
                ooss[std::make_pair(complexId1, complexId2)] = os;
        }
    }

    // Print the pairs of complexes with a high OS.
    // The sorting is by crecient id, a.k.a. by apparition in the original input file.
    std::cout << "#Printing all the pairs of protein complexes with an overlap score (OS) over "
              << args.minOverlapScore << "\n" << std::endl;

    for (std::map<std::pair<ComplexId, ComplexId>, double>::const_iterator it = ooss.begin(); it != ooss.end(); ++it) {
        const std::pair<ComplexId, ComplexId>& complexIdPair = it->first;
        double os = it->second;

        std::cout << "With OS=" << os << " :\n";
        std::cout << "\tComplex_" << complexIdPair.first << " : " << references[complexIdPair.first] << "\n";
        std::cout << "\tComplex_" << complexIdPair.second << " : " << references[complexIdPair.second] << "\n";
    }
    std::cout << "#\t" << ooss.size() << " pairs of complexes found" << std::endl;


    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * The next does use of the Templatized C++ Command Line Parser (TCLAP) library, in include/ directory.
 *   http://tclap.sourceforge.net/manual.html
 */
CmdLineArgs
processCmdLine(int argc, char* argv[]) {

    //// Define the main command line object //////////////////////////////////////////////////////////////////////
    TCLAP::CmdLine cmd("Print all pairs of overlapping protein complexes in a reference collection"
                       /*+*/ " -- Carlos Mella, Cecilia Hernandez",
                                    // Message to be displayed in the USAGE output
                       ' ',         // Character used to separate the argument flag/name from the value
                       "1",         // Version number to be displayed by the --version switch
                       false);      // Whether or not to create the automatic --help and --version switches


    //// Arguments are separate objects, added to the CmdLine object one at a time ////////////////////////////////

    // Value args defines a flag and a type of value that it expects
	TCLAP::ValueArg<unsigned int> minComplexSizeArg(
        "s",
        "min-size",
        "Minimum size for the protein complexes that will be considered; defaults to 2.",
        false,
        2,
        "MINIMUM_COMPLEX_SIZE",
        cmd);
	TCLAP::ValueArg<double> minOverlapScoreArg(
        "",
        "os",
        "Minimum overlap score (OS) considered to match two reference complexes; defaults to 0.2.",
        false,
        0.2,
        "MINIMUM_OVERLAP_SCORE",
        cmd);

    // Unlabeled value args aren't identified by a flag, instead they are identified by their position in the argv
    // array. Note that:
    //  - the order that they are added here to the cmd object is the order that they will be parsed.
    //  - only one optional UnlabeledValueArg is possible; it must be the last listed.
    TCLAP::UnlabeledValueArg<std::string> referenceComplexesFileNameArg(
        "REFERENCE_COMPLEXES_FILE", // A one word name for the argument, used only for identification
        "Path to an input text file defining a collection of protein complexes as reference.",
        true,                       // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "REFERENCE_COMPLEXES_FILE", // A short description of the value type, displayed in the USAGE output
        cmd);                       // The parser object to add this argument to

    //// Parse the argv array /////////////////////////////////////////////////////////////////////////////////////
    cmd.parse(argc, argv);

    // Extra validation checks
    if (referenceComplexesFileNameArg.getValue().empty())
        throw TCLAP::CmdLineParseException("Empty argument!", referenceComplexesFileNameArg.longID());
    if (minOverlapScoreArg.getValue() == 0.0 || minOverlapScoreArg.getValue() > 1.0)
        throw TCLAP::CmdLineParseException("Value out of range!", minOverlapScoreArg.longID());


    //// Get the value parsed by each argument ////////////////////////////////////////////////////////////////////
    CmdLineArgs args;

    args.referenceComplexesFileName = referenceComplexesFileNameArg.getValue();
    args.minComplexSize = minComplexSizeArg.getValue();
    args.minOverlapScore = minOverlapScoreArg.getValue();

    return args;
}
