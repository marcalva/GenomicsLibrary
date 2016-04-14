//
//  ASQuantifier.h
//  
//
//  Created by Marcus Alvarez on 12/10/15.
//
//

#ifndef _ASQuantifier_
#define _ASQuantifier_

#include <stdio.h>
#include <fstream>
#include <string>
#include "vcf.h"
#include "GenomicInterval.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

namespace ASQuantifier {
	class ASQuantifier
	{
	public:
		// Default Constructor
		// Set default paramters
		ASQuantifier();
		
		// Open input VCF file stream from ifs
		// Return true if able to open, false if not
		bool openVcf(const std::string& ifs);
		void closeVcf();
		
		bool openBam(const std::string& bamfn);
		void closeBam();
		
		// Store gene information from BED file ifs.
		bool readGenes(const std::string& ifs);
		
		// Output file functions
		bool openOut(const std::string& ofs);
		void closeOut();
		
		
		// Go through snp file, count allele specific reads, and output to file.
		bool countReads();
		
		// Set parameters for counting reads.
		void set5primecutoff(const int& c);
		void set3primecutoff(const int& c);
		void setRefMaxMismatch(const int& c);
		void setAltMaxMismatch(const int& c);
		void setExonic(bool c);
		void setExonMinReadFrac(const double& c);
		
	private:
		// File streams
		vcfFileInStream _vcfifs;
		BamTools::BamReader _bamfs;
		std::ifstream _bedifs;
		std::ofstream _ofs;
		
		// Store gene information
		GenomicInterval _gi;
		
		// Parameters
		int _5primecutoff;
		int _3primecutoff;
		int _refMaxMismatch;
		int _altMaxMisMatch;
		bool _useExonic;
		double _exonMinReadFraction;
		
		// Alignment functions
		bool _readOverlap(const BamTools::BamAlignment& al);
		int _getLocIndex(const BamTools::BamAlignment& al, std::string chr, int loc);
	};
	
}

#endif /* defined(_ASQuantifier_) */
