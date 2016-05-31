// main.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <algorithm> 
#include "LdInfo.h"
#include "SnpStruct.h"
#include "ParseText.h"

using namespace std;

const int DEFAULT_WINDOW = 1000000;
const float DEFAULT_R2 = 0;
const float DEFAULT_DP = 0;

void printUsage()
{
	cerr << endl;
	cerr << "LdInfo" << endl;
	cerr << endl;
	cerr << "##############################################################################################################" << endl;
	cerr << "Output LD measures R2 and D' for a set of SNPs." << endl;
	cerr << "Given a set of SNPs, it will search the VCF file for SNPs in LD, and will output R2 and D'" << endl;
	cerr << "VCF file must be bgziped and tabix-indexed." << endl;
	cerr << "SNP file must be in tab-delimited format: CHR	POS	ID" << endl;
	cerr << "Comment lines can start with #." << endl;
	cerr << "Output file will contain SNP in first three columns, proxy in 4th-6th columns, distance, R2, and D'." << endl;
	cerr << "##############################################################################################################" << endl;
	cerr << endl;
	cerr << "Usage:" << endl;
	cerr << "./LdInfo --vcf [VCF] --snps [SNPS] --out [OUT]" << endl;
	cerr << "Optional:" << endl;
	cerr << "        " << "--restrict [RESTRICT] Can be a list of SNP IDs, one per line, or in the same format as the input snps file." << endl;
	cerr << "        " << "                      Will only check LD against SNPs in this list (instead of all others in the window)" << endl;
	cerr << "        " << "--min-r2 [float]      Minimum R2 value for SNP to be output." << endl;
	cerr << "        " << "--min-dprime [float]  Minimum D' value for SNP to be output. If R2 and D' thresholds are given, both will be used as criteria." << endl;
	cerr << "        " << "--max-dist [int]      Maximum distance for SNP to be analyzed." << endl;
	cerr << "        " << "--2-snps [flag]       If given, will find proxy for 2 snps present in snps file" << endl;
	cerr << "        " << "--snp1-col [int]      Required if --2-snps is given. Column where SNP 1 starts" << endl;
	cerr << "        " << "--snp2-col [int]      Required if --2-snps is given. Column where SNP 2 starts" << endl;
	cerr << endl;
}

int main(int argc, char* argv[])
{
	string vcfFileName, snpsFileName, restrictFileName, outFileName, minR2Str, minDpStr, maxDistStr, snp1colStr, snp2colStr;
	bool twoSnps = false;
	if (argc < 7){
		printUsage();
		exit (EXIT_FAILURE);
	}
	
	// Parse arguments
	for (int i = 1; i < argc; i+=2)
	{
		string str = string(argv[i]);
		string arg = string(argv[i+1]);
		if (str == "--vcf")
			vcfFileName = arg;
		else if (str == "--snps")
			snpsFileName = arg;
		else if (str == "--out")
			outFileName = arg;
		else if (str == "--restrict")
			restrictFileName = arg;
		else if (str == "--min-r2")
			minR2Str = arg;
		else if (str == "--min-dprime")
			minDpStr = arg;
		else if (str == "--max-dist")
			maxDistStr = arg;
		else if (str == "--2-snps"){
			twoSnps = true;
			i -= 1;
		}
		else if (str == "--snp1-col")
			snp1colStr = arg;
		else if (str == "--snp2-col")
			snp2colStr = arg;
		else{
			cerr << "Unrecognized command line argument " << str << endl;
			printUsage();
			exit (EXIT_FAILURE);
		}
	}
	
	int maxDist;
	if (!maxDistStr.empty()){
		for (size_t s = 0; s < maxDistStr.length(); s++){
			if (!isdigit(maxDistStr[s])){
				cerr << "--max-dist must be an integer >= 1, not " << maxDistStr << endl;
				exit(EXIT_FAILURE);
			}
		}
		maxDist = atoi(maxDistStr.c_str());
		if (maxDist < 1){
			cerr << "--max-dist must be an integer >= 1, not " << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{
		maxDist = DEFAULT_WINDOW;
	}

	float minDp;
	if (!minDpStr.empty()){
		bool pointFlag = false;
		for (size_t s = 0; s < minDpStr.length(); s++){
			if (pointFlag && !isdigit(minDpStr[s])){
				cerr << "--min-dprime must be a floating point number between 0 and 1, not " << minDpStr << endl;
				exit(EXIT_FAILURE);
			}
			else if (minDpStr[s] == '.')
				pointFlag = true;
			else if (!isdigit(minDpStr[s])){
				cerr << "--min-dprime must be a floating point number between 0 and 1, not " << minDpStr << endl;
				exit(EXIT_FAILURE);
			}
		}
		minDp = atof(minDpStr.c_str());
		if (minDp < 0 || minDp > 1){
			cerr << "--min-dprime must be a floating point number between 0 and 1, not " << minDpStr << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{
		minDp = DEFAULT_DP;
	}
	
	float minR2;
	if (!minR2Str.empty()){
		bool pointFlag = false;
		for (size_t s = 0; s < minR2Str.length(); s++){
			if (pointFlag && !isdigit(minR2Str[s])){
				cerr << "--min-r2 must be a floating point number between 0 and 1, not " << minR2Str << endl;
				exit(EXIT_FAILURE);
			}
			else if (minR2Str[s] == '.')
				pointFlag = true;
			else if (!isdigit(minR2Str[s])){
				cerr << "--min-r2 must be a floating point number between 0 and 1, not " << minR2Str << endl;
				exit(EXIT_FAILURE);
			}
		}
		minR2 = atof(minR2Str.c_str());
		if (minR2 < 0 || minR2 > 1){
			cerr << "--min-r2 must be a floating point number between 0 and 1, not " << minR2Str << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{
		minR2 = DEFAULT_R2;
	}

	int snp1col, snp2col;
	if (twoSnps){
		for (size_t s = 0; s < snp1colStr.length(); s++){
			if (!isdigit(snp1colStr[s])){
				cerr << "--snp1-col must be an integer greater than 0" << endl;
				exit(EXIT_FAILURE);
			}
		}
		snp1col = atoi(snp1colStr.c_str());
		if (snp1col < 1){
			cerr << "--snp1-col must be an integer greater than 0" << endl;
			exit(EXIT_FAILURE);
		}
		for (size_t s = 0; s < snp2colStr.length(); s++){
			if (!isdigit(snp2colStr[s])){
				cerr << "--snp2-col must be an integer greater than 0" << endl;
				exit(EXIT_FAILURE);
			}
		}
		snp2col = atoi(snp2colStr.c_str());
		if (snp2col < 1){
			cerr << "--snp2-col must be an integer greater than 0" << endl;
			exit(EXIT_FAILURE);
		}
		// 0-based
		snp1col -= 1;
		snp2col -= 1;
	}

	
	string inLine;
	vector<string> inVector;

	cout << "Restricting output of proxy SNPs to a minimum R2 of " << minR2 << ", a minimum D' of " << minDpStr << " and a window size of " << maxDistStr << endl;
	
	std::unordered_set<std::string> restrictSnps;
	if (!restrictFileName.empty())
	{
		ifstream rfs(restrictFileName, std::ifstream::in);
		if (!rfs.is_open()){
			cerr << "Could not open restrictions file " << restrictFileName << endl;
			exit(EXIT_FAILURE);
		}
		while (rfs.good())
		{
			getline(rfs, inLine);
			if (inLine.empty() || inLine[0] == '#')
				continue;
			parse(inLine, '\t', inVector);
			if (inVector.size() == 1)
				restrictSnps.insert(inVector.at(0));
			else
				restrictSnps.insert(inVector.at(2));
		}
		cout << "Restricting output of proxy SNPs to these " << restrictSnps.size() << " in " << restrictFileName << endl;
	}
	
	LinkageDisequilibrium li(vcfFileName);
	
	list<SnpLd> outProxies;
	ifstream inSnpFs(snpsFileName, std::ifstream::in);
	if (!inSnpFs.is_open()){
		cerr << "Could not open SNP file " << snpsFileName << endl;
		exit(EXIT_FAILURE);
	}
	ofstream outSnpFs(outFileName, std::ofstream::out);
	if (!outSnpFs.is_open()){
		cerr << "Could not open output file " << outFileName << endl;
		exit(EXIT_FAILURE);
	}

	// Store input SNPs into vector
	vector<string> ils;
	while (inSnpFs.good())
	{
		string istring;
		getline(inSnpFs, istring);
		if (istring[0] != '#' && !istring.empty())
			ils.push_back(istring);
	}
	cout << endl << "Analyzing " << ils.size() << " input tag SNPs." << endl;
	
	// This is to track progress of the SNPs input
	bool ten = false, twenty = false, thirty = false, forty = false, fifty = false, sixty = false, seventy = false, eighty = false, ninety = false;
	
	if (! twoSnps )
		outSnpFs << "#ProxyChr	ProxyPos	ProxyID	TagChr	TagPos	TagID	R2	D'	Distance" << endl;;
	for (size_t snpN = 0; snpN < ils.size(); ++snpN)
	{
		/////////////////////////////////////////////////////////////
		// Track progress
		/////////////////////////////////////////////////////////////
		float snpN_float = snpN, nSnps = ils.size();
		float pd = snpN_float / nSnps;
		if (!ten && pd > .1){
			cout << "10\% of input snps done..." << endl;
			ten = true;
		}
		else if (!twenty && pd > .2){
			cout << "20\% of input snps done..." << endl;
			twenty = true;
		}
		else if (!thirty && pd > .3){
			cout << "30\% of input snps done..." << endl;
			thirty = true;
		}
		else if (!forty && pd > .4){
			cout << "40\% of input snps done..." << endl;
			forty = true;
		}
		else if (!fifty && pd > .5){
			cout << "50\% of input snps done..." << endl;
			fifty = true;
		}
		else if (!sixty && pd > .6){
			cout << "60\% of input snps done..." << endl;
			sixty = true;
		}
		else if (!seventy && pd > .7){
			cout << "70\% of input snps done..." << endl;
			seventy = true;
		}
		else if (!eighty && pd > .8){
			cout << "80\% of input snps done..." << endl;
			eighty = true;
		}
		else if (!ninety && pd > .9){
			cout << "90\% of input snps done..." << endl;
			ninety = true;
		}

		/////////////////////////////////////////////////////////////
		// Finished tracking progress
		/////////////////////////////////////////////////////////////


		inVector.clear();
		string inLine = ils.at(snpN);
		parse(inLine, '\t', inVector);
		if (inVector.size() < 3)
			continue;
		
		if (!twoSnps){
			SnpBase thisSnp;
			thisSnp.chr = inVector.at(0);
			// Not analyzing SNP if chromosome is not present in the VCF file.
			if (!li.hasChr(thisSnp.chr))
				continue;
			for (size_t s = 0; s < inVector.at(1).length(); s++){
				if (!isdigit(inVector.at(1)[s])){
					cerr << "2nd column of SNP file must be an integer genomic position, not " << inVector.at(1) << endl;
					exit(EXIT_FAILURE);
				}
			}
			thisSnp.pos = atoi(inVector.at(1).c_str());
			thisSnp.name = inVector.at(2);
	
			outProxies.clear();
			if (!restrictFileName.empty())
				li.GetProxies(thisSnp, &outProxies, restrictSnps, minR2, minDp, maxDist);
			else
				li.GetProxies(thisSnp, &outProxies, minR2, minDp, maxDist);
			
			list<SnpLd>::iterator it;
			for (it = outProxies.begin(); it != outProxies.end(); ++it)
			{
				outSnpFs << thisSnp.chr << '\t' << thisSnp.pos << '\t' << thisSnp.name
					<< '\t' << it->chr << '\t' << it->pos << '\t' << it->name
					<< '\t' << it->r2 << '\t' << it->Dprime << '\t' << it->distance << endl;
			}
		}
		else{
			SnpBase thisSnp1;
			if (inVector.size() <= (size_t)(snp2col+2)){
				cerr << "Input line has less columns than given from snp2-col" << endl;
				exit(EXIT_FAILURE);
			}
			thisSnp1.chr = inVector.at(snp1col);
			thisSnp1.pos = atoi(inVector.at(snp1col+1).c_str());
			thisSnp1.name = inVector.at(snp1col+2);
			SnpBase thisSnp2;
			thisSnp2.chr = inVector.at(snp2col);
			thisSnp2.pos = atoi(inVector.at(snp2col+1).c_str());
			thisSnp2.name = inVector.at(snp2col+2);

			outProxies.clear();

			li.GetProxies(thisSnp1, thisSnp2, &outProxies, minR2, minDp, maxDist);

			list<SnpLd>::iterator it;
			for (it = outProxies.begin(); it != outProxies.end(); ++it)
			{
				outSnpFs << it->chr << '\t' << it->pos << '\t' << it->name
					<< '\t' << inLine << endl;
			}
		}
	}
	cout << "All done!" << endl;
	return 0;
}
