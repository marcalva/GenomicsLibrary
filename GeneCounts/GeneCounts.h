//
//  GeneCounts.h
//  
//
//  Created by Marcus Alvarez on 01/12/16.
//  Pajukanta Lab
//
#ifndef __GENE_COUNTS__
#define __GENE_COUNTS__

#include <stdio.h>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <sstream>
#include "Interval.h"
#include "GenomicInterval.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "tabix.hpp"

struct Gene
{
	// Constructor. Sets Count to 0 and ReadNames empty.
	Gene(){};
	Gene(std::string Name);
	
	std::string Name;
	std::string Contig;
	int Start;
	int End;
	std::set<std::string> ReadNames;
	int Count;
	
	bool operator==(const Gene g) const{
		return (Name == g.Name);
	}
	
	bool operator<(const Gene g) const{
		return (Name < g.Name);
	}
};

class GeneCounts
{
public:
	// Default Constructor
	GeneCounts(){};
	
	// Open files
	bool ReadBED(std::string BEDFileName); // Must be tabix indexed
	bool ReadGTF(std::string GTFFileName); // Must be tabix indexed
	bool OpenBAM(std::string BAMFileName);
	
	void CountReads();
	
	bool WriteTable(std::string OutFileName);
	
private:
	std::string _bedFileName;
	BamTools::BamReader _bamfs;
	std::ofstream _ofs;
	
	map<std::string,Gene> _countsTable;
	vector<std::string> _geneOrder; // Preserves initial order when writing genes to file
	
	// Alignment functions
	bool _readOverlap(const BamTools::BamAlignment& al);
	int _getLocIndex(const BamTools::BamAlignment& al, std::string chr, int loc);
};

#endif // __GENE_COUNTS__