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
#include <odsg/utils/strings.hpp>


struct CmdLineArgs {    // The definition of processCmdLine() constains descriptions for each option
    // Input files
    std::string referenceComplexesFileName;
    std::string refpredComplexesFileName;
    std::string referenceComplexesNames;

    // Other options
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
typedef std::set<ComplexId> complexOverlap;
typedef std::map<int, std::set<int> >::iterator cIterator;
typedef std::map<int, ComplexData >::iterator iterData;
std::string fileNameInput="";
std::ofstream summary;
std::ofstream allOverlapData;
std::ofstream namesCrCp;
std::ofstream namesCrCr;
std::ofstream namesCpCp;
std::ofstream namesCpCpMul;

double sum_re =0.0;
double numcpcp=0.0;
double sum_ssr =0.0;
double num_ssr=0.0;
double sum_ssr_ref =0.0;
double num_ssr_ref=0.0;

void printHisto(std::map<float, int> &histo, std::string type, std::string file){
	std::ofstream outfile(file.c_str());
	std::map<float, int>::iterator it;
	outfile<<type<<std::endl;
	for(it=histo.begin(); it!= histo.end(); it++){
		outfile<<it->first<<" "<<it->second<<"\n";
	}
	outfile.close();
}

void updateHisto(std::map<float, int> &histo, float os){
	float pf1 = 0.1, pf2=0.2, pf3=0.3, pf4=0.4, pf5=0.5, pf6=0.6, pf7=0.7, pf8=0.8, pf9=0.9, pf10=1.0;
	if(os >= pf1)
		histo[pf1]++;
	if(os >= pf2)
		histo[pf2]++;
	if(os >= pf3)
		histo[pf3]++;
	if(os >= pf4)
		histo[pf4]++;
	if(os >= pf5)
		histo[pf5]++;
	if(os >= pf6)
		histo[pf6]++;
	if(os >= pf7)
		histo[pf7]++;
	if(os >= pf8)
		histo[pf8]++;
	if(os >= pf9)
		histo[pf9]++;
	if(os == pf10)
		histo[pf10]++;
}

void addToOverlaps(int c1, int c2, std::map<int, std::set<int> > &complexesOverlap){
    cIterator ct = complexesOverlap.find(c1);
    std::set<int> oc;
    if(ct != complexesOverlap.end()){
	oc = ct->second;
    }
    oc.insert(c2);
    complexesOverlap[c1] = oc;
}

void printMap(std::map<int, std::set<int> > & overlap){
    std::set<int>::iterator sit;
    cIterator ct;
    std::set<int> oc;
    for(ct = overlap.begin(); ct!= overlap.end(); ct++){
	std::set<int> oc = ct->second;
	if(oc.size() < 2 )continue;
	std::cout<<"cr "<<ct->first<<" : ";
	for(sit=oc.begin(); sit!=oc.end(); sit++)
		std::cout<<*sit<<" ";

	std::cout<<std::endl;
    }
}

void printData(std::map<int, ComplexData > & data){
    iterData ct;
    ComplexData oc;
    ComplexData::iterator sit;
    std::cout<<" number of complexes "<<data.size()<<"\n";
    for(ct = data.begin(); ct!= data.end(); ct++){
	oc = ct->second;
	std::cout<<"cr "<<ct->first<<" : ";
	for(sit=oc.begin(); sit!=oc.end(); sit++)
		std::cout<<*sit<<" ";

	std::cout<<std::endl;
    }
}

void printCData(ComplexData c){
	ComplexData::iterator s;
	for(s=c.begin(); s!=c.end(); s++)
		std::cout<<*s<<" ";
	std::cout<<'\n';
}

void printCDataFile(std::ofstream & x, ComplexData c){
	ComplexData::iterator s;
	for(s=c.begin(); s!=c.end(); s++)
		x<<*s<<" ";
	x<<'\n';
}

// compute all pcs ids that overlaps with rcs and creates the map pcrcmap that has pcid mapped to rcs ComplexData (used for getting
// complex name to write in namesCpCpMul file done in compOver function
std::vector<int> buildAllpcs(std::vector<int> allcs, std::map<int, std::pair<int, ComplexData> > &refpreds, std::map<int, ComplexData> &predData, std::map<std::pair<int, int>, float> &rpOverlap, std::map<int, ComplexData> &pcrcmap, std::map<int, ComplexData> &rdata){
    std::vector<int> allpcs;
    std::map<int, std::pair<int, ComplexData> >::iterator it;
    for(unsigned int i=0; i<allcs.size(); i++){
	summary<<"\nRC "<<allcs[i]<<" ";
	allOverlapData<<"\nRC "<<allcs[i]<<" ";
	it = refpreds.find(allcs[i]);
	if(it != refpreds.end()){
		std::pair<int, ComplexData> p = it->second;
		allpcs.push_back(p.first);
		pcrcmap[p.first] = rdata[allcs[i]];
		predData[p.first] = p.second;
		summary<<"PC "<<p.first<<" OS "<<rpOverlap[std::make_pair(allcs[i], p.first)];
		allOverlapData<<"PC "<<p.first<<" OS "<<rpOverlap[std::make_pair(allcs[i], p.first)];
	}
    }
    return allpcs;
}

//compute intersection among all complexes in allcs
void compOver(std::vector<int> &allcs, std::map<int, ComplexData> &data, std::map<std::pair<int,int>, float> &os,
	std::map<ComplexData, std::string> &refNames, std::map<int, ComplexData> &pcrcmap){
    ComplexData r1, r2;
    ComplexData::iterator sit;
    float sum_rcrc =0.0, sum_pcpc = 0.0;
    float num_rcs = 0.0, num_pcs = 0.0;
    allOverlapData<<" overlap\n";
    for(unsigned int i=0; i<allcs.size()-1; i++){
	for(unsigned int j=i+1; j<allcs.size(); j++){
		r1 = data[allcs[i]];
		r2 = data[allcs[j]];
		for(sit = r1.begin(); sit != r1.end(); sit++)
			allOverlapData<<*sit<<" ";
		allOverlapData<<std::endl;
		for(sit = r2.begin(); sit != r2.end(); sit++)
			allOverlapData<<*sit<<" ";
		allOverlapData<<std::endl;
		std::vector<std::string> common_data;
		set_intersection(r1.begin(),r1.end(),r2.begin(),r2.end(), std::back_inserter(common_data));
		allOverlapData<<"set intersection \n";
		for(unsigned int k=0; k<common_data.size(); k++)
			allOverlapData<<common_data[k]<<" ";
		allOverlapData<<std::endl;
		// this is for rcs
		if(os.size() > 0){
		    summary<<"RC_RC_OS("<<allcs[i]<<", "<<allcs[j]<<") = "<<os[std::make_pair(allcs[i], allcs[j])]<<"\n";
		    allOverlapData<<"RC_RC_OS("<<allcs[i]<<", "<<allcs[j]<<") = "<<os[std::make_pair(allcs[i], allcs[j])]<<"\n";
		    sum_rcrc += os[std::make_pair(allcs[i], allcs[j])];
		    num_rcs = num_rcs + 1.0;
		    if(os[std::make_pair(allcs[i], allcs[j])] > 0.0 && allcs.size() > 2)
			namesCpCpMul<<std::setprecision(3)<<os[std::make_pair(allcs[i], allcs[j])]<<" #RC1# "<<refNames[r1]<<" #RC2# "<<refNames[r2]<<"\n";
		} else {
			// this is for rcs
		    float os_score = common_data.size()*common_data.size();
		    os_score = os_score/(r1.size()*r2.size());
		    summary<<"PC_PC_OS("<<allcs[i]<<", "<<allcs[j]<<") = "<<os_score<<"\n";
		    allOverlapData<<"PC_PC_OS("<<allcs[i]<<", "<<allcs[j]<<") = "<<os_score<<"\n";
		    sum_pcpc += os_score;
		    num_pcs = num_pcs + 1.0;
		    if(os_score > 0.0 && allcs.size() > 2)
			namesCpCpMul<<std::setprecision(3)<<os_score<<" #PC1# "<<refNames[pcrcmap[allcs[i]]]<<" #PC2# "<<refNames[pcrcmap[allcs[j]]]<<"\n";

		}
	}
    }
    if(os.size() > 0){
      sum_ssr_ref += sum_rcrc/num_rcs;
      num_ssr_ref += 1.0;
      summary<<"\nsum_rcrc "<<sum_rcrc<<" num_rcs "<<num_rcs<<" RC_RC_OS_AVG = "<<sum_rcrc/num_rcs<<"\n";
      allOverlapData<<"\nsum_rcrc "<<sum_rcrc<<" num_rcs "<<num_rcs<<" RC_RC_OS_AVG = "<<sum_rcrc/num_rcs<<"\n";
    } else { // for PCs os map is empty
      sum_ssr += sum_pcpc/num_pcs;
      num_ssr += 1.0;
      summary<<"\nsum_pcpc "<<sum_pcpc<<" num_pcs "<<num_pcs<<" PC_PC_OS_AVG = "<<sum_pcpc/num_pcs<<"\n";
      allOverlapData<<"\nsum_pcpc "<<sum_pcpc<<" num_pcs "<<num_pcs<<" PC_PC_OS_AVG = "<<sum_pcpc/num_pcs<<"\n";
    }
}

void printVector(std::vector<int> v){
	for(unsigned int i=0; i<v.size(); i++)
		std::cout<<v[i]<<" ";
	std::cout<<"\n";
}

void addMaximalSet(std::vector<int> candidate, std::vector< std::vector<int> > &allcsVector){
    bool found = false;
    sort(candidate.begin(), candidate.end());
    for(unsigned int i=0; i<allcsVector.size(); i++){
	if(allcsVector[i].size() >= candidate.size()){
	    if(includes(allcsVector[i].begin(), allcsVector[i].end(), candidate.begin(), candidate.end())){
		//std::cout<<" found \n";
		//printVector(allcsVector[i]);
		found = true;
		return;
            }
	} else {
	    if(includes(candidate.begin(), candidate.end(), allcsVector[i].begin(), allcsVector[i].end())){
		allcsVector.erase(allcsVector.begin()+i);
		allcsVector.push_back(candidate);
		return;
            }
	}
    }
    if(!found){
	allcsVector.push_back(candidate);
    }
}

void computeOverlaps(std::map<int, std::set<int> > & overlap, std::map<int, ComplexData> & data, std::map<int, std::pair<int, ComplexData > > &refpreds, std::map<std::pair<int, int>, float> &refOverlap, std::map<std::pair<int,int>, float> &rpOverlap,
	std::map<ComplexData, std::string> &refNames){
    std::set<int>::iterator sit;
    cIterator ct;
    std::set<int> oc;
    std::vector< std::vector<int> > allcsVector;
    for(ct = overlap.begin(); ct!= overlap.end(); ct++){
	std::set<int> oc = ct->second;
        std::vector<int> allcs;
	allcs.push_back(ct->first);
	for(sit=oc.begin(); sit!=oc.end(); sit++){
		allcs.push_back(*sit);
	}
	addMaximalSet(allcs, allcsVector);
    }

    summary<<"\nComplexes in Overlap Analysis using Reference Complexes total "<<allcsVector.size()<<"\n";
    allOverlapData<<"\nComplexes in Overlap Analysis using Reference Complexes total "<<allcsVector.size()<<"\n";
    unsigned int matchings = 0;
    unsigned int largestPcs = 0;

    for(unsigned int x=0; x<allcsVector.size(); x++){
	std::vector<int> allcs = allcsVector[x];
	//printVector(allcsVector[x]);
	summary<<"\nReference Complexes in Overlap Set "<<x<<" size "<<allcs.size()<<"\n";
	if(allcs.size() > 2){
		namesCpCpMul<<"------------------------------------------------------\n";
		namesCpCpMul<<"\nReference Complexes in Overlap Set "<<x<<" size "<<allcs.size()<<"\n";
		namesCpCpMul<<"------------------------------------------------------\n";
	}
	allOverlapData<<"\nReference Complexes in Overlap Set "<<x<<" size "<<allcs.size()<<"\n";
	allOverlapData<<"\nCOMPARISON AMONG Reference Complexes (RC)\n";
	summary<<"\nCOMPARISON AMONG Reference Complexes (RC)\n";
	std::map<int, ComplexData> pcrcmap;
	compOver(allcs, data, refOverlap, refNames, pcrcmap);
	std::map<int, ComplexData> maprefpreds;
	summary<<"\nReference Complexes (RC) MMR mapping at OS=0.2 Predicted Complexes (PC)\n";
	allOverlapData<<"\nReference Complexes (RC) MMR mapping at OS=0.2 to Predicted Complexes (PC)\n";
	std::vector<int> allpcs = buildAllpcs(allcs,refpreds,maprefpreds, rpOverlap, pcrcmap, data);
	std::map<std::pair<int,int>, float> rpO;
	if(allcs.size() > 2)
		namesCpCpMul<<"------------------------------------------------------\n";
	if(allpcs.size() > 1){
	    allOverlapData<<"\nCOMPARISON AMONG Predicted Complexes (PC)\n";
	    summary<<"\nCOMPARISON AMONG Predicted Complexes in (PC)\n";
	    compOver(allpcs, maprefpreds, rpO, refNames, pcrcmap);
	    matchings++;
	    if(allpcs.size() > largestPcs) largestPcs = allpcs.size();
	} else {
	    allOverlapData<<"\nReference complexes do not have at least 2 matching in Predicted Complexes\n";
	    summary<<"\nReference complexes do not have at least 2 matching in Predicted Complexes\n";
	}
    }
    summary<<"\nTotal PCs in overlap = "<<matchings<<" largest set of PCs in overlap = "<<largestPcs<<"\n";
    allOverlapData<<"\nTotal PCs in overlap = "<<matchings<<" largest set of PCs in overlap = "<<largestPcs<<"\n";
}


std::map<ComplexData, std::string> readRefNamesFile(const std::string& fileName){

    ComplexData::iterator sit;
    std::map<ComplexData, std::string > refComplexNames;
    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readRefComplexesNameFile(): can not open input file");
    }

    std::cout<<"\n reading file "<<fileName<<std::endl;
    std::string line, prot;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')       // TODO: drop lines with empty spaces too
            continue;

	std::string nameComplex;
        std::istringstream iss(line);

        ComplexData complexData;
	int flag = 0;
        for (std::string token; iss >> token; ) {
                std::size_t found = token.find("#");
                if(found == std::string::npos){
		    if(flag)
			complexData.insert(token);
		    else
                    	nameComplex += token + " ";
                } else {
		    flag=1;
                    prot= token.substr(found+1);
		    std::string str = token.substr(0,token.length()-(token.length()-found));
		    nameComplex += str;
		    complexData.insert(prot);
		    continue;
		}
        }
        refComplexNames[complexData] = nameComplex;
    }
    return refComplexNames;
}

std::map<int, std::pair<int, ComplexData> >
readRefPredComplexesFile(const std::string& fileName, double minOverlapScore, std::map<std::pair<int,int>, float> & rpoverlap,
	std::map< ComplexData, std::string> &refNames) {
    std::map<float, int> histoCrCp;
    histoCrCp[0.1] = histoCrCp[0.2] = histoCrCp[0.3] = histoCrCp[0.4] = histoCrCp[0.5] = 0;
    histoCrCp[0.6] = histoCrCp[0.7] = histoCrCp[0.8] = histoCrCp[0.9] = histoCrCp[1.0] = 0;
    std::map<int, std::pair<int, ComplexData> > refPredComplexes;
    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readRefPredComplexesFile(): can not open input file");
    }

    std::cout<<"\n reading file "<<fileName<<std::endl;

    std::map<int, ComplexData> complexesData;
    std::set<int>::iterator sit;

    unsigned int currentLine = 0;
    std::string line, osstr;
    int complex1_id = -1, complex2_id = -1;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')       // TODO: drop lines with empty spaces too
            continue;

        std::istringstream iss(line);

        ComplexData complexData1, complexData2;
        if(currentLine % 3 == 0 ){
            for (std::string token; iss >> token; ) {
                std::size_t found = token.find("=");
                if(found != std::string::npos){
                    osstr = token.substr(found+1);
                }
            }
        } else if(currentLine % 3 == 1){ // complex 1
            for (std::string token; iss >> token; ) {
                std::size_t found = token.find("_");
                if(found != std::string::npos){
                    complex1_id = atoi(token.substr(found+1).c_str());
                } else if(token.find(":")){
                        complexData1.insert(token);
                }
            }
            complexesData[complex1_id] = complexData1;
       } else if(currentLine % 3 == 2){ // complex 2
            for (std::string token; iss >> token; ) {
                std::size_t found = token.find("_");
                if(found != std::string::npos){
                    complex2_id = atoi(token.substr(found+1).c_str());
                } else if(token.find(":")){
                        complexData2.insert(token);
                }
            }
            //complexesData[complex2_id] = complexData2;
            if(osstr != "" && complex1_id != -1 && complex2_id != -1){
		float osf = atof(osstr.c_str());
                //std::cout<<" os "<<osstr<<" cr1 "<<complex1_id<<" cr2 "<<complex2_id<<" osf "<<osf<<std::endl;
		//namesCrCp<<" & "<<std::setprecision(3)<<osf<<" & "<<refNames[complexesData[complex1_id]]<<"\\"<<"\\"<<'\n';
		namesCrCp<<" & "<<std::setprecision(3)<<osf<<" & "<<refNames[complexesData[complex1_id]];
		namesCrCp<<" & ";
		printCDataFile(namesCrCp,complexData2);	
		namesCrCp<<"\\"<<"\\"<<'\n';
		updateHisto(histoCrCp,osf);
		refPredComplexes[complex1_id] = std::make_pair(complex2_id, complexData2);
		rpoverlap[std::make_pair(complex1_id, complex2_id)] = atof(osstr.c_str());
            }
        }
        currentLine++;
    }


    std::ostringstream ss;
    ss << minOverlapScore;
    std::string osinput = ss.str();
    std::string histoname = fileName + ".refOS-" + osinput +".histoCrCp";
    fileNameInput = fileName + ".refOS-" + osinput;
    printHisto(histoCrCp, "CrCp monotic crecient histogram", histoname);
    return refPredComplexes;
}

void computeHistoCpCp(std::map<std::pair<int, int>, float> &refOverlap, std::map<int, std::pair<int, ComplexData> > &refpreds,
		std::map<int, ComplexData> &data, std::map<ComplexData, std::string> & refNames, double minOS){

    std::map<float, int> histoCpCp;
    std::map<std::pair<int, int>, float>::iterator ref;
    std::map<int, std::pair<int, ComplexData> >::iterator pred1, pred2;
    for(ref = refOverlap.begin(); ref!= refOverlap.end(); ref++){
	int cr1 = ref->first.first;
	int cr2 = ref->first.second;
	pred1 = refpreds.find(cr1);
	ComplexData d1, d2;
	if(pred1 != refpreds.end()){
		std::pair<int, ComplexData> p = pred1->second;
		d1 = p.second;
		pred2 = refpreds.find(cr2);
		if(pred2 != refpreds.end()){
			std::pair<int, ComplexData> p2 = pred2->second;
			d2 = p2.second;
			ComplexData common_data, common_refs, diff;
			odsg::algorithms::set_intersection(d1,d2,common_data);

			float os = common_data.size()*common_data.size();
			os = os/(d1.size()*d2.size());
			updateHisto(histoCpCp, os);
			if(os >= minOS){
				odsg::algorithms::set_intersection(data[cr1],data[cr2],common_refs);
				std::set_symmetric_difference(common_refs.begin(), common_refs.end(), common_data.begin(), common_data.end(), std::inserter(diff, diff.begin()));

				if(common_refs.size() > 0)
					sum_re += 1.0*(diff.size())/common_refs.size();
				numcpcp++;
				namesCpCp<<"& "<<std::setprecision(3)<<os<<" & "<<refNames[data[cr1]]<<" & "<<refNames[data[cr2]]<<"\\"<<"\\"<<"\n";
				namesCpCp<<"PC1 : ";
				printCDataFile(namesCpCp, d1);
				namesCpCp<<"PC2 : ";
				printCDataFile(namesCpCp, d2);
				namesCpCp<<"Common PC1-PC2 : ";
				printCDataFile(namesCpCp, common_data);
				namesCpCp<<"Common RC1-RC2 : ";
				printCDataFile(namesCpCp, common_refs);
				namesCpCp<<"|Co sym_diff Po| : "<<diff.size()<<"\n";
			}
		}
	}
    }

    std::string histoname = fileNameInput + ".histoCpCp";
    printHisto(histoCpCp, "CpCp monotonic crecient histogram", histoname);

}

void readComplexesFromFile(const std::string& fileName, double minOverlapScore, std::map<int, std::pair<int, ComplexData> > &refpreds, std::map<std::pair<int, int>, float> &rpOverlap, std::map< ComplexData, std::string> & refNames) {

    std::map<float, int> histoCrCr;
    std::map<int, std::set<int> > complexesOverlap1;
    std::map<int, std::set<int> > complexesOverlap2;
    std::map<std::pair<int, int>, float> refOverlap;
    cIterator ct;
    std::ifstream infile(fileName.c_str());
    if (!infile) {
        throw std::runtime_error("readComplexesFromFile(): can not open input file");
    }

    std::cout<<"\n reading file "<<fileName<<std::endl;
    summary<<"\n reading file "<<fileName<<std::endl;
    allOverlapData<<"\n reading file "<<fileName<<std::endl;

    std::map<int, ComplexData> complexesData;
    std::set<int>::iterator sit;

    unsigned int currentLine = 0;
    std::string line, osstr;
    int complex1_id = -1, complex2_id = -1;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')       // TODO: drop lines with empty spaces too
            continue;

        std::istringstream iss(line);

        ComplexData complexData1, complexData2;
 	if(currentLine % 3 == 0 ){
            for (std::string token; iss >> token; ) {
	        std::size_t found = token.find("=");
	        if(found != std::string::npos){
		    osstr = token.substr(found+1);
	        }
            }
	} else if(currentLine % 3 == 1){ // complex 1
            for (std::string token; iss >> token; ) {
	        std::size_t found = token.find("_");
	        if(found != std::string::npos){
		    complex1_id = atoi(token.substr(found+1).c_str());
	        } else if(token.find(":")){
			complexData1.insert(token);
		}
            }
	    complexesData[complex1_id] = complexData1;

	} else if(currentLine % 3 == 2){ // complex 2
            for (std::string token; iss >> token; ) {
	        std::size_t found = token.find("_");
	        if(found != std::string::npos){
		    complex2_id = atoi(token.substr(found+1).c_str());
	        } else if(token.find(":")){
			complexData2.insert(token);
		}
            }
	    complexesData[complex2_id] = complexData2;
	    if(osstr != "" && complex1_id != -1 && complex2_id != -1){
		//std::cout<<" os "<<osstr<<" cr1 "<<complex1_id<<" cr2 "<<complex2_id<<std::endl;
		addToOverlaps(complex1_id, complex2_id, complexesOverlap1);
		addToOverlaps(complex2_id, complex1_id, complexesOverlap1);
		refOverlap[std::make_pair(complex1_id, complex2_id)] = atof(osstr.c_str());
		updateHisto(histoCrCr,atof(osstr.c_str()));
		namesCrCr<<" & "<<std::setprecision(3)<<atof(osstr.c_str())<<" & "<<refNames[complexesData[complex1_id]]<<" & "<<refNames[complexesData[complex2_id]]<<"\\"<<"\\"<<"\n";
		ComplexData common;
		odsg::algorithms::set_intersection(complexesData[complex1_id], complexesData[complex2_id], common);
		namesCrCr<<"RC1 : ";
		printCDataFile(namesCrCr, complexesData[complex1_id]);
		namesCrCr<<"RC2 : ";
		printCDataFile(namesCrCr, complexesData[complex2_id]);
		namesCrCr<<"Common : ";
		printCDataFile(namesCrCr, common);
	    }
	}
	currentLine++;
    }

    computeOverlaps(complexesOverlap1, complexesData, refpreds, refOverlap, rpOverlap, refNames);
    std::ostringstream ss;
    ss << minOverlapScore;
    std::string osinput = ss.str();
    std::string histoname = fileName + ".refOS-"+osinput+".histoCrCr";
    printHisto(histoCrCr, "CrCr monotonic crecient histogram", histoname);
    computeHistoCpCp(refOverlap, refpreds, complexesData, refNames, minOverlapScore);

    return;
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

    std::cout << "Considering protein complexes with a OS of " << args.minOverlapScore<< std::endl;
    summary << "Considering protein complexes with a OS of " << args.minOverlapScore<< std::endl;
    allOverlapData << "Considering protein complexes with a OS of " << args.minOverlapScore<< std::endl;


    // Read the references from file
    std::cout << "Reading reference complexes file... " << std::flush;
    summary << "Reading reference complexes file... " << std::flush;
    allOverlapData << "Reading reference complexes file... " << std::flush;
    std::map<ComplexId, ComplexData> references;
    std::map<int, std::pair<int, ComplexData > > refpreds;
    std::map<ComplexData, std::string> refNames;
    std::map<std::pair<int, int>, float > rpOverlap;
    try {
	std::string namesCrCpFile = args.refpredComplexesFileName + ".refOS-" + odsg::strings::to_str(args.minOverlapScore) + ".namesCrCp";
	namesCrCp.open(namesCrCpFile.c_str(), std::ofstream::out);
	std::string namesCpCpFile = args.refpredComplexesFileName + ".refOS-" + odsg::strings::to_str(args.minOverlapScore) + ".namesCpCp";
	namesCpCp.open(namesCpCpFile.c_str(), std::ofstream::out);
	std::string namesCrCrFile = args.refpredComplexesFileName + ".refOS-" + odsg::strings::to_str(args.minOverlapScore) + ".namesCrCr";
	namesCrCr.open(namesCrCrFile.c_str(), std::ofstream::out);
	std::string namesCpCpMulFile = args.refpredComplexesFileName + ".refOS-" + odsg::strings::to_str(args.minOverlapScore) + ".namesCpCpMul";
	namesCpCpMul.open(namesCpCpMulFile.c_str(), std::ofstream::out);
        refNames = readRefNamesFile(args.referenceComplexesNames);
        refpreds = readRefPredComplexesFile(args.refpredComplexesFileName, args.minOverlapScore, rpOverlap, refNames);
	std::string summaryFile = fileNameInput + ".summary";
	summary.open(summaryFile.c_str(), std::ofstream::out);
	std::string allOverlapDataFile = fileNameInput + ".allOverlapData";
	allOverlapData.open(allOverlapDataFile.c_str(), std::ofstream::out);
        readComplexesFromFile(args.referenceComplexesFileName, args.minOverlapScore, refpreds, rpOverlap, refNames);
    } catch (std::exception& e) {
        std::cerr << "ERROR\n"
                  << e.what() << std::endl;
        return 1;
    }
    std::cout << "OK\n"
              << '\t' << "sum_re "<<sum_re<<" T ="<<numcpcp<<" Average Relative Error = '" << sum_re/numcpcp << "'\n";
    std::cout << "|SSR|_ref = "<<num_ssr_ref<<" AvgOSOverlap_ref = "<<sum_ssr_ref/num_ssr_ref<<"\n";
    std::cout << "|SSR| = "<<num_ssr<<" AvgOSOverlap = "<<sum_ssr/num_ssr_ref<<"\n";

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
    TCLAP::CmdLine cmd("Generate cummulative overlap histogram and others overlap metrics"
                       /*+*/ " -- Carlos Mella, Cecilia Hernandez",
                                    // Message to be displayed in the USAGE output
                       ' ',         // Character used to separate the argument flag/name from the value
                       "1",         // Version number to be displayed by the --version switch
                       false);      // Whether or not to create the automatic --help and --version switches


    //// Arguments are separate objects, added to the CmdLine object one at a time ////////////////////////////////

    TCLAP::ValueArg<double> minOverlapScoreArg(
    "",
    "os",
    "Minimum overlap score (OS) considered to match two reference complexes; defaults to 0.2.",
    false,
    0.2,
    "MINIMUM_OVERLAP_SCORE",
    cmd);

    TCLAP::UnlabeledValueArg<std::string> refpredComplexesFileNameArg(
        "REFERENCE_PREDICTED_COMPLEX_FILE", // A one word name for the argument, used only for identification
        "Path to an input text file defining MMR matching protein complexes reference to predicted complexes.",
        true,                       // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "REFERENCE_PREDICTED_COMPLEXES_FILE", // A short description of the value type, displayed in the USAGE output
        cmd);                       // The parser object to add this argument to

    TCLAP::UnlabeledValueArg<std::string> referenceComplexesFileNameArg(
        "REFERENCE_COMPLEXES_FILE", // A one word name for the argument, used only for identification
        "Path to an input text file defining a collection of protein complexes as reference.",
        true,                       // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "REFERENCE_COMPLEXES_FILE", // A short description of the value type, displayed in the USAGE output
        cmd);                       // The parser object to add this argument to

    TCLAP::UnlabeledValueArg<std::string> referenceComplexesNamesNameArg(
        "REFERENCE_COMPLEXES_FILE_NAMES", // A one word name for the argument, used only for identification
        "Path to an input text file defining a collection of protein complexes as reference with their names.",
        true,                       // Whether the argument is required on the command line
        "",                         // Default value of this argument; unused if the presence of the arg is required
        "REFERENCE_COMPLEXES_FILE_NAMES", // A short description of the value type, displayed in the USAGE output
        cmd);                       // The parser object to add this argument to

    //// Parse the argv array /////////////////////////////////////////////////////////////////////////////////////
    cmd.parse(argc, argv);

    // Extra validation checks
    if (referenceComplexesFileNameArg.getValue().empty())
        throw TCLAP::CmdLineParseException("Empty argument!", referenceComplexesFileNameArg.longID());
    if (refpredComplexesFileNameArg.getValue().empty())
        throw TCLAP::CmdLineParseException("Empty argument!", refpredComplexesFileNameArg.longID());
    if (referenceComplexesNamesNameArg.getValue().empty())
        throw TCLAP::CmdLineParseException("Empty argument!", referenceComplexesNamesNameArg.longID());
    if (minOverlapScoreArg.getValue() == 0.0 || minOverlapScoreArg.getValue() > 1.0)
        throw TCLAP::CmdLineParseException("Value out of range!", minOverlapScoreArg.longID());


    //// Get the value parsed by each argument ////////////////////////////////////////////////////////////////////
    CmdLineArgs args;

    args.refpredComplexesFileName = refpredComplexesFileNameArg.getValue();
    args.referenceComplexesFileName = referenceComplexesFileNameArg.getValue();
    args.referenceComplexesNames = referenceComplexesNamesNameArg.getValue();
    args.minOverlapScore = minOverlapScoreArg.getValue();

    return args;
}
