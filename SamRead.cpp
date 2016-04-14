//
//  SamRead.cpp
//  
//
//  Created by Marcus Alvarez on 12/10/15.
//
//

#include "SamRead.h"
#include <string>
#include <cstdlib>
#include <stdio.h>

using namespace std;

SamRead::SamRead(BamAlignment b)
{
	ba = b;
}

bool SamRead::isContained(const GenomicInterval& gi)
{
	// Split alignment into pieces if there are any splices or Ns in CIGAR
	
	// Go through each split alignment and look for overlap of genomic intervals.
	// If there are none at any time, return false
	
	// If we go through all split alignments and an interval is always returned, then return true
	
}

void SamRead::SnpCounts(const int& refid, const int& pos, const char& refAllele, const char& altAllele,
						int& refCounts, int& altCounts, int& refPlus, int& refMinus, int& altPlus, int& altMinus)
{
	refCounts = 0;
	altCounts = 0;
	refPlus = 0;
	refMinus = 0;
	altPlus = 0;
	altMinus = 0;

	// If position is not contained within read, return without modifying variables
	if ( refid != ba.RefID )
		return;
	
	unsigned int zero_pos = pos - 1;
	
	vector<Interval> alignments;
	splitAlignment(alignments);
	
	bool snpInRead = false;
	for (unsigned int i = 0; i < alignments.size(); i++)
	{
		unsigned int start = alignments.at(i).Start();
		unsigned int end = alignments.at(i).End();
		
		if ( zero_pos >= start && zero_pos < end )
			snpInRead = true;
	}
	
	if ( !snpInRead )
		return;
	
	

	
	
}

void SamRead::splitAlignment( vector<Interval>& i )
{
	i.erase();
	
	unsigned int indexPosition = ba.Position;
	
	for (unsigned int i = 0; i < ba.CigarDatasize(); i++)
	{
		// Interval stats
		string c = to_string(ba.RefID);
		unsigned int start = 0;
		unsigned int stop = 0;
		
		char type = ba.cigardata.at(i).Type;
		unsigned int length = ba.cigardata.at(i).Length;
		
		if ( type == 'H' || type == 'S' )
			continue;
		else if ( type == 'M' )
		{
			start = indexPosition;
			stop = indexPosition + length;
			indexPosition += length;
			
			i.push_back( Interval(c, start, top) );
		}
		else
		{
			indexPosition += length;
		}
		
	}
}