
#include <iostream>
#include <string>
#include <cstdlib>
#include "ASQuantifier.h"

using namespace std;

void printUsage()
{
	cout << endl << endl;
	cout << "#############################################" << endl;
	cout << "ASQuantifier" << endl;
	cout << "#############################################" << endl << endl << endl;
	cout << "Usage:" << endl;
	cout << "./ASQuantifier [Options] --vcf VCFFILENAME --bam BAMFILENAME --out OUTFILENAME" << endl;
	cout << "Optional arguments:" << endl;
	cout << "\t--bed GENEFILENAME\tBED file with genomic coordinates. If given, will only output allele counts from SNPs that fall within these coordinates." << endl;
	cout << "\t--use-exonic\t\tIf this flag is given, only reads that fall in gene bed coordinates will be used. Default false" << endl;
	cout << "\t--5prime-cut INT\tSNP must be more than this number of bases away from the 5' end of the alignment. Default 2." << endl;
	cout << "\t--3prime-cut INT\tSNP must be more than this number of bases away from the 3' end of the alignment. Default 2." << endl;
	cout << "\t--ref-max-mismatch INT\tRead with the reference allele at this snp must have no more than this number of mismatches in NM field. Default 10." << endl;
	cout << "\t--alt-max-mismatch INT\tRead with the alternate allele at this snp must have no more than this number of mismatches in NM field. Default 11." << endl;
	cout << endl << endl;
	
}

int main(int argc, char* argv[])
{
	if ( argc < 6 )
	{
		cerr << "You must give the variant, bam, and output files" << endl;
		printUsage();
		exit(1);
	}
	
	string bamFileName, vcfFileName, outName, bedName;
	int fivePrime = 2;
	int threePrime = 2;
	int refMaxMismatch = 10;
	int altMaxMismatch = 11;
	bool useExonic = false;
	
	int k = 1;
	while ( k < argc )
	{
		if ( strcmp(argv[k], "--vcf") == 0 )
			vcfFileName = argv[k+1];
		else if ( strcmp(argv[k], "--bam") == 0 )
			bamFileName = argv[k+1];
		else if ( strcmp(argv[k], "--out" ) == 0 )
			outName = argv[k+1];
		else if ( strcmp(argv[k], "--bed" ) == 0 )
			bedName = argv[k+1];
		else if ( strcmp(argv[k], "--use-exonic" ) == 0 )
		{
			useExonic = true;
			k--;
		}
		else if ( strcmp(argv[k], "--5prime-cut" ) == 0 )
			fivePrime = atoi(argv[k+1]);
		else if ( strcmp(argv[k], "--3prime-cut" ) == 0 )
			threePrime = atoi(argv[k+1]);
		else if ( strcmp(argv[k], "--ref-max-mismatch" ) == 0 )
			refMaxMismatch = atoi(argv[k+1]);
		else if ( strcmp(argv[k], "--alt-max-mismatch" ) == 0 )
			altMaxMismatch = atoi(argv[k+1]);
		else
		{
			cout << "Unknown argument " << argv[k] << endl;
			printUsage();
			exit(1);
		}
		k += 2;
	}
	
	bool check = true;
	
	ASQuantifier::ASQuantifier asq;
	
	check = asq.openBam(bamFileName);
	if ( !check )
	{
		cout << "Could not open Bam file." << endl;
		exit(1);
	}
	
	check = asq.openVcf(vcfFileName);
	if ( !check )
	{
		cout << "Could not open VCF file." << endl;
		exit(1);
	}
	
	check = asq.openOut(outName);
	if ( !check )
	{
		cout << "Could not open output file." << endl;
		exit(1);
	}
	
	if ( !bedName.empty() )
	{
		check = asq.readGenes(bedName);
		if ( !check )
		{
			cout << "Could not open BED file." << endl;
			exit(1);
		}
	}
	
	asq.setExonic(useExonic);
	asq.set5primecutoff(fivePrime);
	asq.set3primecutoff(threePrime);
	asq.setRefMaxMismatch(refMaxMismatch);
	asq.setAltMaxMismatch(altMaxMismatch);

	check = asq.countReads();
	if ( !check )
	{
		cout << "Failed to count allelic reads." << endl;
		exit(1);
	}
	
	asq.closeVcf();
	asq.closeBam();
	asq.closeOut();
	return 0;
}
