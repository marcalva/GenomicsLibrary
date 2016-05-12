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

void printUsage()
{
	cerr << endl;
	cerr << "LdInfo" << endl;
	cerr << endl;
	cerr << "Output LD measures R2 and D' for a set of SNPs." << endl;
	cerr << "Given a set of SNPs, it will search the VCF file for SNPs in LD, and will output R2 and D'" << endl;
	cerr << "VCF file must be bgziped and tabix-indexed." << endl;
	cerr << "SNP file must be in tab-delimited format: CHR	POS	ID" << endl;
	cerr << "Comment lines can start with #." << endl;
	cerr << "Output file will contain SNP in first three columns, proxy in 4th-6th columns, distance, R2, and D'." << endl;
	cerr << endl;
	cerr << "Usage:" << endl;
	cerr << "./LdInfo --vcf [VCF] --snps [SNPS] --out [OUT]" << endl;
	cerr << "Optional:" << endl;
	cerr << "        " << "--restrict [RESTRICT] Can be a list of SNP IDs, one per line, or in the same format as the input snps file." << endl;
	cerr << "        " << "                      Will only check LD against SNPs in this list (instead of all others in the window)" << endl;
	cerr << "        " << "--min-r2 [float]      Minimum R2 value for SNP to be output." << endl;
	cerr << "        " << "--max-dist [int]      Maximum distance for SNP to be analyzed." << endl;
}

int main(int argc, char* argv[])
{
	string vcfFileName, snpsFileName, restrictFileName, outFileName, minR2Str, maxDistStr;
	
	if (argc < 7){
		printUsage();
		exit (EXIT_FAILURE);
	}
	
	// Parse arguments
	for (int i = 1; i < argc; i++)
	{
		string str = string(argv[i]);
		if (str == "--vcf")
			vcfFileName = str;
		else if (str == "--snps")
			snpsFileName = str;
		else if (str == "--out")
			outFileName = str;
		else if (str == "--restrict")
			restrictFileName = str;
		else if (str == "--min-r2")
			minR2Str = str.c_str();
		else if (str == "--max-dist")
			maxDistStr = str;
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
	
	string inLine;
	vector<string> inVector;
	
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
			else if (inVector.size() == 3)
				restrictSnps.insert(inVector.at(2));
		}
	}
	
	LinkageDisequilibrium li(vcfFileName);
	
	list<SnpLd> outProxies;
	ifstream inSnpFs(snpsFileName, std::ifstream::in);
	if (!inSnpFs.is_open()){
		cerr << "Could not open SNP file " << snpsFileName << endl;
		exit(EXIT_FAILURE);
	}
	ofstream outSnpFs(outFileName, std::ofstream::out);
	if (outSnpFs.is_open()){
		cerr << "Could not open output file " << outFileName << endl;
		exit(EXIT_FAILURE);
	}
	
	while (inSnpFs.good())
	{
		getline(inSnpFs, inLine);
		if (inLine.empty() || inLine[0] == '#')
			continue;
		parse(inLine, '\t', inVector);
		if (inVector.size() < 3)
			continue;
		
		SnpBase thisSnp;
		thisSnp.chr = inVector.at(0);
		// Not analyzing and outputting SNPs whose chromosome is not present in the VCF file.
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
		
		if (!restrictFileName.empty())
			li.GetProxies(thisSnp, &outProxies, restrictSnps, minR2, maxDist);
		else
			li.GetProxies(thisSnp, &outProxies, minR2, maxDist);
		
		list<SnpLd>::iterator it = outProxies.begin();
		for (; it != outProxies.end(); ++it)
		{
			outSnpFs << thisSnp.chr << '\t' << thisSnp.pos << '\t' << thisSnp.name
				<< '\t' << it->chr << '\t' << it->pos << '\t' << it->name
				<< '\t' << it->r2 << '\t' << it->Dprime << endl;
		}
	}
}