// Create VCF file from input VCF and ASQ output
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>

using namespace std;

int ASQ_NAME_COL = 2;
int ASQ_GT = 5;
int ASQ_A1_PLUS = 6;
int ASQ_A1_MINUS = 7;
int ASQ_A2_PLUS = 8;
int ASQ_A2_MINUS = 9;

struct AsCounts
{
public:
	unsigned int ref;
	unsigned int alt;
};

std::vector<std::string>& splitLine(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


int main(int argc, char* argv[])
{
	if ( argc < 6 )
	{
		cerr << "You must give the variant, bam, and output files" << endl;
		printUsage();
		exit(1);
	}
	
	string vcfName, dirName, ofName, idfName, asqName;
	
	int k = 1;
    while ( k < argc )
    {
		if ( strcmp(argv[k], "--vcf") == 0 )
			vcfName = argv[k+1];
		else if ( strcmp(argv[k], "--dir") == 0 )
			dirName = argv[k+1];
		else if ( strcmp(argv[k], "--ids") == 0 )
			idfName = argv[k+1];
		else if ( strcmp(argv[k], "--asq" ) == 0 )
			asqName = argv[k+1];
		else if ( strcmp(argv[k], "--out" ) == 0 )
			ofName = argv[k+1];
		else
		{
			cout << "Unknown argument " << argv[k] << endl;
			printUsage();
			exit(1);
		}
		k += 2;
	}
	
	/** Read in allelic Counts information into a matrix. Store 2 separate vectors of sample IDs and SNP names. **/
	vector<string> sampleNames;
	map<string, int> si;
	vector<string> snpNames;
	vector< vector< AsCounts> > counts; // Samples x SNPs
	
	ifstream sampleIfs(idfName.c_str(), ios_base::openmode mode = ios_base::in);
	string line;
	int index = 0;
	while (getline(sampleIfs, line))
	{
		sampleNames.push_back(line);
		si[line] = index;
		index++;
	}
	sampleIfs.close();
	counts.resize(sampleNames.size());
	for (int i = 0; i < sampleNames.size(); i++)
	{
		asCountsFn = dirName + "/" + sampleNames.at(i) + "/" + asqName;
		ifstream asCountsFs(asCountsFn, ios_base::openmode mode = ios_base::in);
		getline(asCountsFs, line); // Read out header
		line = "";
		while (getline(asCountsFs, line))
		{
			vector<string> lineV;
			splitLine(line, '\t', lineV);
			string name = lineV.at(ASQ_NAME_COL);
			string gt = lineV.at(ASQ_GT);
			unsigned int a1plus = atoi(lineV.at(ASQ_A1_PLUS).c_str());
			unsigned int a1minus = atoi(lineV.at(ASQ_A1_MINUS).c_str());
			unsigned int a2plus = atoi(lineV.at(ASQ_A2_PLUS).c_str());
			unsigned int a2minus = atoi(lineV.at(ASQ_A2_MINUS).c_str());
			
			unsigned int ref;
			unsigned int alt;
			
			// Homozygous reference
			if (gt == "0|0" || gt == "0/0"){
				ref = a1plus + a1minus;
				alt = a2plus + a2minus;
			}
			else if (gt == "0|1" || gt == "0/1"){
				ref = a1plus + a1minus;
				alt = a2plus + a2minus;
			}
			else if (gt == "1|0"){
				ref = a2plus + a2minus;
				alt = a1plus + a1minus;
			}
			else if (gt == "1|1" || gt == "1/1"){
				ref = a1plus + a1minus;
				alt = a2plus + a2minus;
			}
			else {
				cout << "Unkown genotype " << gt << " at SNP " << name << endl;
				exit(1);
			}
			
			AsCounts c;
			c.ref = ref;
			c.alt = alt;
			if (i == 0){
				snpNames.push_back(name);
			}
			counts.at(i).push_back(c);
		}
		asCountsFs.close();
	} // Stored allele specific counts
	
	string vcf_asHeader = "##FORMAT=<ID=AS,Number=.,Type=String,Description=\"Allele-specific expression counts from RNA-seq\">";
	
	ifstream vcfFs(vcfName, ios_base::openmode mode = ios_base::in);
	ofstream outFs(ofName, ios_base::openmode mode = ios_base::out);
	line = "";
	
	vector<string> vcfSamples;
	int vcfLineNumber = 0;
	while (getline(vcfFs, line)){
		if ( line.at(0) == '#' && line.at(1) == '#' ){
			outFs << line << endl;
		}
		else if ( line.at(0) == '#' && line.at(1) == 'C' ){
			outFs << vcf_asHeader << endl;
			outFs << line << endl;
			vector<string> lineV;
			lineSplit(line, '\t', lineV);
			for (int i = 9; i < lineV.size(); i++)
			{
				vcfSamples.push_back(lineV.at(i));
			}
		}
		else{
			vector<string> lineV;
			lineSplit(line, '\t', lineV);
			
			outFs << lineV.at(0) << '\t' << lineV.at(1) << '\t' << lineV.at(2) << '\t' << lineV.at(3) <<
				'\t' << lineV.at(4) << '\t' << lineV.at(5) << '\t' << lineV.at(6) << '\t' << lineV.at(7) << 
				'\t' << lineV.at(8) << ":AS";
			
			if (lineV.at(2) != snpNames.at(vcfLineNumber)){
				for (int i = 9; i < lineV.size(); i++)
					outFs << '\t' << lineV.at(i) << ":0,0";
				outFs << endl;
				continue;
			}
			else{
				for (int i = 9; i < lineV.size(); i++)
				{
					int vcfI = i - 9;
					int countsIndex = si[vcfSamples.at(vcfI)];
					unsigned int refC = counts.at(countsIndex).at(vcfLineNumber).ref;
					unsigned int altC = counts.at(countsIndex).at(vcfLineNumber).alt;
					outFs << '\t' << lineV.at(i) << ':' << refC << ',' << altC;
				}
				outFs << endl;
				vcfLineNumber++;
			}
		}
	}
	vcfFs.close();
	outFs.close();
}
