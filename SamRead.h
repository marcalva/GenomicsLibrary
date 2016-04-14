//
//  SamRead.h
//  
//
//  Created by Marcus Alvarez on 12/10/15.
//
//

#ifndef _SamRead_h
#define _SamRead_h

#include <vector>
#include "api/BamAlignment.h"
#include "GenomicInterval.h"
#include "Interval.h"

// Wrapper for Bamtools BamAlignment
class SamRead
{
public:
	Bamtools::BamAlignment ba;
	
	
	// Constructor
	SamRead(BamAlignment b);
	
	// Determine whether read is fully contained within the set of genomic intervals given
	bool isContained(const GenomicInterval& gi);
	
	// Return Snp Counts
	void SnpCounts(const int& refid, const int& pos, const char& refAllele, const char& altAllele,
				   int& refCounts, int& altCounts, int& refPlus, int& refMinus, int& altPlus, int& altMinus);
	
private:
	void splitAlignment( vector<Interval>& i );
	
};

#endif
