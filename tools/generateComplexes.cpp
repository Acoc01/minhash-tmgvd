#include <cassert>      // Support run-time assertions. They can be disabled defining the NDEBUG macro
#include <exception>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <ctime>                           // for timing


#include <tclap/CmdLine.h>

#include <odsg/utils/algorithms.hpp>
#include <odsg/DagForest.hpp>
#include <odsg/DenseSubGraphsMaximalSet.hpp>
#include <odsg/Graph.hpp>
#include <odsg/Vertex.hpp>
#include <odsg/VertexSet.hpp>

#include <typedefs.hpp>
#include <readFile.hpp>
#include <metrics.hpp>
#include <printing.hpp>
#include <min.hpp>

using namespace odsg;
using namespace bio_odsg;


struct CmdLineArgs {    // The definition of processCmdLine() constains descriptions for each option
    // Input files
    std::string datasetMappingFileName;
    std::string datasetFileName;

    // Options related to the way that complexes are generated
    int partitioning;
    std::string outlinksSorting;
    bool cliquesOnly;

    bool weightedDataset;

    std::string weightDensityMetric;
    unsigned int objective;     // Not exposed, dependent of weightDensityMetric

    std::string extendedLogFileName;

    // Options related to the way that generated complexes are treated
    unsigned int minComplexSize;
    int similarityFiltering;
};
CmdLineArgs processCmdLine(int argc, char* argv[]);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Build a single proteins set from a pair <Sources, Centers> from a dense subgraph.
 */
VertexSet
unifyDenseSubGraph(const DenseSubGraph& dsg, bool treatAsClique) {
    VertexSet complex = dsg.getCenters();
    if (!treatAsClique)
        complex.insert(dsg.getSources().begin(), dsg.getSources().end());

    return complex;
}

VertexSet centers_nodes(const DenseSubGraph& dsg){
    VertexSet cent = dsg.getCenters();
    return cent;
}

VertexSet source_nodes(const DenseSubGraph& dsg){
    VertexSet src = dsg.getSources();
    return src;
}

bool isbiclique(const DenseSubGraph& dsg){
    return dsg.biClique();
}


/*
 * Add a new generated (o 'predicted') complex to a collection of these.
 * In some cases, many generated complexes shares a great similarity, so these can be filtered/merged,
 * according to the filteringBySimilarity flag; the accepted values are:
 *   - 0: No filtering/merging.
 *   - 1: Having two complexes very similar, the biggest is kept and the other is forgotten/removed. This is the
 *        implementation done by ClusterBFS algorithm; see its paper for further details.
 *   - 2: Having two complexes very similar, both are merged and kept as one; This is done at least by clusterONE.
 *
 * The 'big picture' implementation is from ClusterBFS.
 */

void
//addToFilteredPredictedComplexes(std::vector<Complex>& complexes,
addToFilteredPredictedComplexes(std::set<Complex>& complexes,
                                const Complex& complex,
                                int filteringBySimilarity,
                                double filteringBySimilarityThreshold=0.8) {
    assert(!complex.empty());
    assert(0 <= filteringBySimilarity && filteringBySimilarity <= 2);
    assert(0.0 <= filteringBySimilarityThreshold && filteringBySimilarityThreshold <= 1.0);

    if (filteringBySimilarity == 0) {   // No filtering
        //complexes.push_back(complex);
        complexes.insert(complex);
        return;
    }

    // Search in the collection the most similar to complex
    //std::vector<Complex>::iterator mostSimilarComplex = complexes.end();
    std::set<Complex>::iterator mostSimilarComplex = complexes.end();
    double maxSimilarScore = -1.0;

    //for (std::vector<Complex>::iterator it = complexes.begin(); it != complexes.end(); ++it) {
    for (std::set<Complex>::iterator it = complexes.begin(); it != complexes.end(); ++it) {
        double score = overlapScore(*it, complex);
        if (score > maxSimilarScore) {
            mostSimilarComplex = it;
            maxSimilarScore = score;
        }
    }

    if (maxSimilarScore < filteringBySimilarityThreshold) {
        // The complex is enought different of all the existent complexes. OK adding it.
        //complexes.push_back(complex);
        complexes.insert(complex);
    } else {
        if (filteringBySimilarity == 1) {
            if (complex.size() > mostSimilarComplex->size()) {
                complexes.erase(mostSimilarComplex);
                //complexes.push_back(complex);
                complexes.insert(complex);
            }
        } else if (filteringBySimilarity == 2) {
            Complex merged;
            std::set_union(complex.begin(), complex.end(),
                           mostSimilarComplex->begin(), mostSimilarComplex->end(),
                           std::inserter(merged, merged.end()));

            complexes.erase(mostSimilarComplex);
            //complexes.push_back(merged);
            complexes.insert(merged);
        }
    }
}


/*
 * Output messages
 */
std::string
outlinksSortingMsg(std::string sorting) {
    if (sorting == "ID") {
        return "Sorting of the adjacency lists of the graph by protein id";
    } else if (sorting == "FREQUENCY") {
        return "Sorting of the adjacency lists of the graph by frequency";
    } else {    // Unknown
        return "Sorting of the adjacency lists of the graph by... Uh?";
    }
}

std::string
graphTypeMsg(bool weighted) {
    if (weighted)
        return "Dataset treated as an weighted graph";
    else
        return "Dataset treated as an unweighted graph";
}

std::string
weightDensityMsg(std::string densitymetric) {
    if (densitymetric == "WEDGE")
        return "Mining of dense subgraphs optimizing by weighted-edge-density (WEDGE)";
    if (densitymetric == "WDEGREE")
        return "Mining of dense subgraphs optimizing by weighted-degree-density (WDEGREE)";
    if (densitymetric == "FWDEGREE")
        return "Mining of dense subgraphs optimizing by weighted-degree-density considering all possibles edges (FWDEGREE)";
    if (densitymetric == "FWEDGE")
        return "Mining of dense subgraphs optimizing by weighted-edge-density considering all possibles edges (FWEDGE)";
    if (densitymetric == "DEGREE_WEDGE")
        return "Mining of dense subgraphs optimizing by degree-density and weighted-edge-density (DEGREE_WEDGE)";

    return "Mining of dense subgraphs ignoring weights (if present)";
}

std::string
partitioningMsg(int partitioning) {
    switch (partitioning) {
        case 0:     return "No graph partitioning";
        case 1:     return "Graph partitioning according to common initial outlink";
        case 2:     return "Graph partitioning according to common signature via 'shingles'";
        default:    return "Graph partitioning according to... Uh?";
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(int argc, char* argv[]) {

    clock_t start, finish;
    double time;
    start = clock();


    CmdLineArgs args;
    try {
        args = processCmdLine(argc, argv);
    } catch (TCLAP::ArgException& e) {
        std::cerr << "error: " << e.error() << " " << e.argId() << std::endl;
        return 1;
    }

    // All messages are sent to cerr; cout is reserved to the list of predicted complexes
    std::cerr << "Considering protein complexes with a minimum size of " << args.minComplexSize << "\n";


    // Initialize the (optional) extended log file
    bool extendedLogging = !args.extendedLogFileName.empty();

    std::ofstream extendedLogFile(args.extendedLogFileName.c_str());
    if (extendedLogging && !extendedLogFile) {
        std::cerr << "error: can not open file for extended log";
        return 1;
    }
    //Vector que contiene los clusters obtenidos con el minhash
    std::vector<std::vector<int>> v1 = minhash::min(args.datasetFileName);
    //Grafo con la listas de cada nodo para poder mapear los nodos de un cluster con todas sus aristas
    std::vector<std::vector<long long int>> gadj = minhash::graph;
    //Vector de WGraph, uno por cada cluster
    std::vector<WGraph> datasetWGraph;
    //Vector de punteros de WGraph para poder armar los dagForest
    std::vector<Graph*> datasetGraph_ptr;
    //Mapeo de cluster a grafo en un archivo que se crea
    //Imprime un cluster con sus nodos y listas, luego agrega un # para indicar el final del cluster
    std::ofstream myfile;
    myfile.open("clusters.txt");
    for(int i = 0; i < v1.size(); ++i){
        int flag = 0;
        for(int j = 0; j < v1[i].size(); ++j){
            int node = v1[i][j];
            int flag_sloop = 0;
            if(gadj[node-1].size() < 3) flag = 1; //Esto es para considerar el minimo tamanio especificado en el algoritmo original
            for(int k = 0; k < gadj[node-1].size(); ++k){
                if(node == gadj[node-1][k])flag_sloop =1; //Si el nodo tenia un self loop, no lo agregamos
                if(k < gadj[node-1].size())myfile<<node<<' '<<gadj[node-1][k]<<' '<<1.0<<'\n';
            }
            if(!flag_sloop)myfile<<node<<' '<<node<<' '<<1.0<<'\n';//Si el nodo no tenia self loop, se lo agregamos.
        }if(flag==0)myfile<<"#\n";
    }
    myfile.close();
    // This mapping will apply to all the proteins seen from now
    ProteinsMap proteinMapping;


    // Get the (optional) custom mapping, previous to the dataset parsing
   if (!args.datasetMappingFileName.empty()) {
        std::cerr << "\nReading dataset mapping file... " << std::flush;
        try {
            proteinMapping = readDatasetMappingFromFile(args.datasetMappingFileName);
        } catch (std::exception& e) {
            std::cerr << "ERROR\n" << e.what() << std::endl;
            return 1;
        }
        std::cerr << "OK\n"
                  << '\t' << proteinMapping.size() << " different proteins\n"
                  ;
    }
    // Aqui modifique la funcion readDatasetFromFileWW ubicada en dapg_complexes/src/readFile.cpp
    datasetWGraph = readDatasetFromFileWW("clusters.txt", args.weightedDataset, proteinMapping);
    std::cout<<args.outlinksSorting<<std::endl;
    //Introducimos los WGraph al vector de punteros
    for(int i = 0; i < datasetWGraph.size(); ++i){
        datasetGraph_ptr.push_back(&datasetWGraph[i]);
    }
    //Preparamos los grafos para ser minados
    if(args.outlinksSorting == "ID"){
        for(int i = 0; i < datasetGraph_ptr.size();++i){
            datasetGraph_ptr[i]->rebuildForMining(Graph::VertexComparer());
        }
    }else{
        for(int i = 0; i < datasetGraph_ptr.size();++i){
            datasetGraph_ptr[i]->rebuildForMining();
        }
    }
    for(int i = 0; i < datasetGraph_ptr.size(); ++i){
        if(datasetGraph_ptr[i]->empty()){
            std::cerr << i <<"One of the dataset graph after rebuilding for mining is empty." << std::endl;
            return 1;
        }
    }
    //Definimos contadores y vectores para guardar cantidad y elementos.
    int cliques = 0; // Contador de cliques
    int biclique_r = 0; // Contador de Bicliques Rigurosos (Interseccion de S y C vacia)
    int biclique_nr = 0; // Contador de Bicliques no Rigurosos (Interseccion no vacia pero S != C)
    std::vector<DenseSubGraph> vector_cliques;
    std::vector<DenseSubGraph> vector_bicliques;
    std::vector<DenseSubGraph> vector_bicliques_no_riguroso;

    std::cout<<"Generating DagForests"<<std::endl;
    for(int i = 0; i < datasetGraph_ptr.size();++i){
        const DagForest myForest(*datasetGraph_ptr[i],args.partitioning,1);
        //std::cout<<"DagForest for cluster "<<i+1<<" ready"<<std::endl;
        unsigned int totalDagsBuilt = 0;
        for (DagForest::const_iterator dit = myForest.begin(); dit != myForest.end(); ++dit) {
            const Dag& dag = **dit;
            const DenseSubGraphsMaximalSet dagDSGs = dag.getDenseSubGraphs(0,
                                                                        args.objective,
                                                                        args.cliquesOnly,
                                                                            1);
            if (dagDSGs.empty())
                continue;
            // Indica el numero de dense sub graphs por cluster
            std::cout<<"Number of Dense Sub Graphs for Cluster "<<i+1<<' '<<dagDSGs.size()<<'\n';
            int cont = 1;
            for (std::vector<DenseSubGraph>::const_iterator it = dagDSGs.begin(); it != dagDSGs.end(); ++it) {
                if(it->biClique()){
                    vector_bicliques.push_back(*it);
                    biclique_r++;
                }
                if(it->clique()){
                    vector_cliques.push_back(*it);
                    cliques++;
                }
                if(!it->biClique() && !it->clique()){
                    vector_bicliques_no_riguroso.push_back(*it);
                    biclique_nr++;
                }
                cont++;
            }

            totalDagsBuilt += dagDSGs.size();
        }
    }
    /*La seccion que sigue se encarga de escribir en un archivo por separado los sets S y C dependiendo
    de si es un Biclique, Clique o Biclique no Riguroso, tarda demasiado en escribir todo porque debe
    hacer el mapeo inverso con getProteinName que se encuentra en typedefs.hpp
    
    Así que es recomendable omitir esto a la hora de testear o limitar la cantidad de elementos a escribir.*/
    std::ofstream res_file;
    res_file.open("results.txt");
    res_file << "BICLIQUES\n";
    for(int i = 0; i < vector_bicliques.size(); ++i){
        VertexSet S,C;
        S = vector_bicliques[i].getSources();
        C = vector_bicliques[i].getCenters();
        for(VertexSet::const_iterator it = S.begin(); it != S.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
        for(VertexSet::const_iterator it = C.begin(); it != C.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
    }res_file<<'\n';
    res_file << "BICLIQUE NO RIGUROSO\n";
    for( int i = 0; i < vector_bicliques_no_riguroso.size(); ++i){
        VertexSet S,C;
        S = vector_bicliques_no_riguroso[i].getSources();
        C = vector_bicliques_no_riguroso[i].getCenters();
        for(VertexSet::const_iterator it = S.begin(); it != S.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
        for(VertexSet::const_iterator it = C.begin(); it != C.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
    }res_file << '\n';
    for(int i = 0; i < vector_cliques.size(); ++i){
        VertexSet S,C;
        S = vector_cliques[i].getSources();
        C = vector_cliques[i].getCenters();
        for(VertexSet::const_iterator it = S.begin(); it != S.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
        for(VertexSet::const_iterator it = C.begin(); it != C.end(); ++it){
            res_file << getProteinName(proteinMapping,*it)<< ' ';
        }res_file << '\n';
    }res_file << '\n';
    res_file.close();
    // Print data to cout
    //printComplexes(predictedComplexes, proteinMapping, std::cout);

    finish = clock();
    std::cout<< "Cliques: " << cliques << '\n';
    std::cout<< "Bicliques Rigurosos: " << biclique_r << '\n';
    std::cout<< "Bicliques no Rigurosos: " << biclique_nr << '\n';
    time = double(finish - start) / CLOCKS_PER_SEC;
    std::cerr << "\tDAPG execution time : " << time << '\n';

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
    TCLAP::CmdLine cmd("Predict protein complexes in PPI networks based on the mining of overlapped dense subgraphs"
                       /*+*/ " using prefix dags -- Carlos Mella, Cecilia Hernandez, Jaime Araya",
                                    // Message to be displayed in the USAGE output
                       ' ',         // Character used to separate the argument flag/name from the value
                       "1",         // Version number to be displayed by the --version switch
                       false);      // Whether or not to create the automatic --help and --version switches


    //// Arguments are separate objects, added to the CmdLine object one at a time ////////////////////////////////

    // Switch args are boolean arguments to define a flag; its presence negates its default value
    TCLAP::SwitchArg unifiedArg(
        "u",
        "unified",
        "<internal> Build a unique ('unified') prefix dag from where to extract all the dense subgraphs."
            " It's a deprecated option: use -p NONE  instead.",
        cmd,
        false);
    TCLAP::SwitchArg denseSubgraphsArg(
        "d",
        "dense-subgraphs",
        "<internal> Expand the mining to general dense subgraphs. IGNORED option: it's enabled anyway."
            "Use the option -c  instead to change it.",
        cmd,
        true);
    TCLAP::SwitchArg cliquesOnlyArg(
        "c",
        "cliques-only",
        "<internal> Limit the mining to dense subgraphs with maximal centers sets between themselves, i.e. cliques.",
        cmd,
        false);

    // Value args defines a flag and a type of value that it expects
    TCLAP::ValueArg<std::string> datasetMappingFileNameArg(
        "m",                        // The one character flag that identifies this argument on the command line
        "mapping",                  // One word name for the argument. Can be used as a long flag
        "<internal> Path to a input text file containing a custom mapping from protein names to internal ids;"
            " if it's not specified, a default internal mapping will assign ids incrementally by apparition in the"
            " dataset input file.",
        false,                      // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "DATASET_MAPPING_FILE",     // A short description of the value type, displayed in the USAGE output
        cmd);
    TCLAP::ValueArg<std::string> extendedLogFileNameArg(
        "e",
        "extended-log-file",
        "<internal> Path to a output text file that will contain a more detailed log about the generation process.",
        false,
        "",
        "EXTENDED_LOG_FILE",
        cmd);
    TCLAP::ValueArg<unsigned int> minComplexSizeArg(
        "s",
        "min-size",
        "Minimum size for the generated protein complexes; defaults to 3."
            " Anyway, complexes with size bigger than 100 are ignored too; it's not configurable (yet)",
        false,
        3,
        "MINIMUM_COMPLEX_SIZE",
        cmd);

    std::vector<std::string> outlinksSortingValues;
    outlinksSortingValues.push_back("ID");
    outlinksSortingValues.push_back("FREQUENCY");
    TCLAP::ValuesConstraint<std::string> outlinksSortingConstraint(outlinksSortingValues);
    TCLAP::ValueArg<std::string> outlinksSortingArg(
        "r",
        "outlinks-sorting",
        "<internal> Select how the adjacency lists from the PPI graph are sorted previous to its partitioning"
            "in clusters."
            " ID sorts by the internal protein id, assigned by the internal mapping (see -m option);"
            " FREQUENCY sorts by apparition frequency in the adjacency lists (i.e. inlinks count)."
            " Defaults to FREQUENCY.",
        false,
        "FREQUENCY",
        &outlinksSortingConstraint,
        cmd);

    std::vector<std::string> graphTypeValues;
    graphTypeValues.push_back("UNONE");
    graphTypeValues.push_back("USYM");
    TCLAP::ValuesConstraint<std::string> graphTypeConstraint(graphTypeValues);
    TCLAP::ValueArg<std::string> graphTypeArg(
        "g",
        "graph-types",
        "Select how the weights for the protein interactions, present the dataset given as input file, are treated."
            " UNONE will ignore the weights (if present);"
            " USYM will include the weights (if missing, they are assumed as 1.0)."
            " Defaults to UNONE.",
        false,
        "UNONE",
        &graphTypeConstraint,
        cmd);

    std::vector<std::string> weightDensityValues;
    weightDensityValues.push_back("WEDGE");
    weightDensityValues.push_back("FWEDGE");
    weightDensityValues.push_back("WDEGREE");
    weightDensityValues.push_back("FWDEGREE");
    weightDensityValues.push_back("DEGREE_WEDGE");
    TCLAP::ValuesConstraint<std::string> weightDensityConstraint(weightDensityValues);
    TCLAP::ValueArg<std::string> weightDensityArg(
        "w",
        "weight-density-objective-function",
        "<internal> Select the weight-density objective function:"
            " WEDGE optimizes by greater weighted-edge-density;"
            " FWEDGE optimizes by greater weighted-edge-density considering all possible edges;"
            " WDEGREE optimizes by greater weighted-degree-density;"
            " FWDEGREE optimizes by greater weighted-degree-density considering all possible edges;"
            " DEGREE_WEDGE optimizes first by greater degree-density and then by weighted-edge-density."
            " Defaults to WDEGREE. It's ignored if -g NONE  or -c  are used.",
        false,
        "WDEGREE",
        &weightDensityConstraint,
        cmd);

    std::vector<std::string> partitioningValues;
    partitioningValues.push_back("NONE");
    partitioningValues.push_back("HASHING");
    partitioningValues.push_back("INITIAL_OUTLINK");
    TCLAP::ValuesConstraint<std::string> partitioningConstraint(partitioningValues);
    TCLAP::ValueArg<std::string> partitioningArg(
        "p",
        "partitioning-scheme",
        "<internal> Select how the set of the adjacency lists from the PPI graph is partitioned in 'clusters',"
            " from which each prefix dags is built:"
            " NONE prevents all partitioning (same as the -u option);"
            " INITIAL_OUTLINK groups adjacency lists sharing the same initial outlink after sorting (see -r option);"
            " HASHING groups adjacency lists with the same 'shingling' signature. Defaults to INITIAL_OUTLINK;"
            " it's ignored if the -u option is used.",
        false,
        "INITIAL_OUTLINK",
        &partitioningConstraint,
        cmd);

    std::vector<std::string> similarityFilteringValues;
    similarityFilteringValues.push_back("NONE");
    similarityFilteringValues.push_back("BIGGEST");
    similarityFilteringValues.push_back("UNION");
    TCLAP::ValuesConstraint<std::string> similarityFilteringConstraint(similarityFilteringValues);
    TCLAP::ValueArg<std::string> similarityFilteringArg(
        "f",
        "filter-similar-predicted",
        "Select how pairs of generated protein complexes that are very similar to each other are treated, with"
            " such 'very similar' property granted if the overlap score (OS) is bigger than the (fixed) value 0.8:"
            " NONE does nothing, keeping both (it will not remove duplicated predicted complexes);"
            " BIGGEST kept only the bigger complex;"
            " UNION kept only a complex made from merging both complexes. Defaults to UNION.",
        false,
        "UNION",
        &similarityFilteringConstraint,
        cmd);


    // Unlabeled value args aren't identified by a flag, instead they are identified by their position in the argv
    // array. Note that:
    //  - the order that they are added here to the cmd object is the order that they will be parsed.
    //  - only one optional UnlabeledValueArg is possible; it must be the last listed.
    TCLAP::UnlabeledValueArg<std::string> datasetFileNameArg(
        "DATASET_FILE",             // A one word name for the argument, used only for identification
        "Path to an input text file defining (probably weighted) protein-protein interactions: a PPI network.",
        true,                       // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "DATASET_FILE",             // A short description of the value type, displayed in the USAGE output
        cmd);                       // The parser object to add this argument to

    //// Parse the argv array /////////////////////////////////////////////////////////////////////////////////////
    cmd.parse(argc, argv);


    // Extra validation checks
    if (datasetFileNameArg.getValue().empty())
        throw TCLAP::CmdLineParseException("Empty argument!", datasetFileNameArg.longID());

    //// Get the value parsed by each argument ////////////////////////////////////////////////////////////////////
    CmdLineArgs args;

    args.datasetMappingFileName = datasetMappingFileNameArg.getValue();
    args.datasetFileName        =        datasetFileNameArg.getValue();
    args.outlinksSorting        =        outlinksSortingArg.getValue();
    args.cliquesOnly            =            cliquesOnlyArg.getValue();
    args.extendedLogFileName    =    extendedLogFileNameArg.getValue();
    args.minComplexSize         =         minComplexSizeArg.getValue();

    args.weightedDataset = graphTypeArg.getValue() == "USYM";
    args.weightDensityMetric = (graphTypeArg.getValue() == "USYM") ? weightDensityArg.getValue() : "";

    args.objective = 2;      // Best objetive function for un-weighted graphs
    if (args.weightDensityMetric == "WEDGE")
        args.objective = 3;
    else if (args.weightDensityMetric == "WDEGREE")
        args.objective = 4;
    else if (args.weightDensityMetric == "DEGREE_WEDGE")
        args.objective = 5;
    else if (args.weightDensityMetric == "FWEDGE")
        args.objective = 6;
    else if (args.weightDensityMetric == "FWDEGREE")
        args.objective = 7;

    args.similarityFiltering = 0;
    if (similarityFilteringArg.getValue() == "BIGGEST")
        args.similarityFiltering = 1;
    else if (similarityFilteringArg.getValue() == "UNION")
        args.similarityFiltering = 2;

    args.partitioning = 0;
    if (!unifiedArg.getValue() && partitioningArg.getValue() == "INITIAL_OUTLINK")
        args.partitioning = 1;
    else if (!unifiedArg.getValue() && partitioningArg.getValue() == "HASHING")
        args.partitioning = 2;

    return args;
}
