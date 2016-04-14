//
//  ASQuantifier.cpp
//  
//
//  Created by Marcus Alvarez on 12/10/15.
//  Pajukanta Lab
//

#include "ASQuantifier.h"
#include <cstdlib>
#include <cctype>
#include <list>
#include <cmath>
#include <climits>

using namespace std;

const int WORD_SIZE = 65536;

ASQuantifier::ASQuantifier::ASQuantifier()
{
	_5primecutoff = 2;
	_3primecutoff = 2;
	_refMaxMismatch = 10;
	_altMaxMisMatch = 11;
	_useExonic = false;
	_exonMinReadFraction = -1;
}

bool ASQuantifier::ASQuantifier::openVcf(const std::string& ifs)
{
	return _vcfifs.openVcf(ifs);
}

void ASQuantifier::ASQuantifier::closeVcf()
{
	_vcfifs.close();
}

bool ASQuantifier::ASQuantifier::openBam(const std::string& bamfn)
{
	_bamfs.Open(bamfn);
	
	if ( !_bamfs.IsOpen() )
	{
		cout << "Could not open BAM file " << bamfn << endl;
		return false;
	}
	
	string bami = bamfn + ".bai";
	if ( ! _bamfs.OpenIndex(bami) )
	{
		cout << "Creating index for " << bamfn << endl;
		_bamfs.CreateIndex();
	}
	
	if ( ! _bamfs.HasIndex() )
	{
		cout << "No access to BAM index" << endl;
		return false;
	}
	return true;
}

void ASQuantifier::ASQuantifier::closeBam()
{
	_bamfs.Close();
}

bool ASQuantifier::ASQuantifier::readGenes(const std::string& ifs)
{
	ifstream fs(ifs);
	
	if (fs.fail())
		return false;
	
	// Read in BED file, skip lines that start with #
	while (fs.good())
	{
		char bedline[WORD_SIZE];
		bedline[0] = '\0';
		fs.getline(bedline, WORD_SIZE, '\n');
		
		if ( bedline[0] == '\0' )
			break;
		
		if ( bedline[0] == '#' )
			continue;
		
		// Create interval
		string c;
		int start, end;
		
		// Store chromosome
		unsigned int i = 0;
		while ( bedline[i] != ' ' && bedline[i] != '\t' )
		{
			c += bedline[i];
			i++;
		}
		
		while ( bedline[i] == ' ' || bedline[i] == '\t' )
			i++;
		
		// Store start
		unsigned int j = 0;
		char number[WORD_SIZE];
		while ( bedline[i] != ' ' && bedline[i] != '\t' )
		{
			if (!isdigit(bedline[i]))
			{
				cout << "BED file of genes does not have a numeric value: " << bedline[i] << " in the 2nd column for start position" << endl;
				exit(1);
			}
			number[j] = bedline[i];
			i++;
			j++;
		}
		number[j] = '\0';
		start = atoi(number);
		
		while ( bedline[i] == ' ' || bedline[i] == '\t' )
			i++;
		
		
		// Store end
		j = 0;
		while ( bedline[i] != ' ' && bedline[i] != '\t' )
		{
			if (!isdigit(bedline[i]))
			{
				cout << "BED file of genes does not have a numeric value: " << bedline[i] << "in the 3nd column for end position" << endl;
				exit(1);
			}
			number[j] = bedline[i];
			i++;
			j++;
		}
		number[j] = '\0';
		end = atoi(number);
		
		Interval interval(c, start, end);
		_gi.addInterval(interval);
	}
	fs.close();
	return true;
}

bool ASQuantifier::ASQuantifier::openOut(const std::string& ofs)
{
	_ofs.open(ofs);
	
	if (_ofs.fail())
		return false;
	return true;
}

void ASQuantifier::ASQuantifier::closeOut()
{
	_ofs.close();
}

void ASQuantifier::ASQuantifier::set5primecutoff(const int& c)
{
	if (c < 0)
	{
		cout << "5 prime cutoff must be value between 0 and read length, not " << c << endl;
		exit(1);
	}
	_5primecutoff = c;
}

void ASQuantifier::ASQuantifier::set3primecutoff(const int& c)
{
	if (c < 0)
	{
		cout << "3 prime cutoff must be value between 0 and read length, not " << c << endl;
		exit(1);
	}
	_3primecutoff = c;
}

void ASQuantifier::ASQuantifier::setRefMaxMismatch(const int& c)
{
	if (c < 0)
	{
		cout << "Max number of reference allele mismatches must be greater than 0, not " << c << endl;
		exit(1);
	}
	_refMaxMismatch = c;
}

void ASQuantifier::ASQuantifier::setAltMaxMismatch(const int& c)
{
	if (c < 1)
	{
		cout << "Max number of alternate allele mismatches must be greater than 1, not " << c << endl;
		exit(1);
	}
	_altMaxMisMatch = c;
}

void ASQuantifier::ASQuantifier::setExonic(bool c)
{
	_useExonic = c;
}

void ASQuantifier::ASQuantifier::setExonMinReadFrac(const double& c)
{
	if ( c < 0 || c > 1)
	{
		cout << "Read fraction covering exon must be between 0 and 1, not " << c << endl;
		exit(1);
	}
	_exonMinReadFraction = c;
}

bool ASQuantifier::ASQuantifier::countReads()
{
	
	////////////////////////////////////////////////////////////////////////
	// Make sure all input and output files are open
	////////////////////////////////////////////////////////////////////////
	if ( !_vcfifs.good() || !_bamfs.IsOpen() || !_ofs.good() )
	{
		cout << "Files are not ready" << endl;
		return false;
	}
	
	
	
	////////////////////////////////////////////////////////////////////////
	// Loop through VCF file
	////////////////////////////////////////////////////////////////////////
	vcfHeader vh = _vcfifs.readHeader();
	while (_vcfifs.good())
	{
		vcfRecord vr = _vcfifs.readRecord();
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		// Skip if not SNP or doesn't fall in exon (when gene bed file is given)
		////////////////////////////////////////////////////////////////////////
		if (vr.ref.size() > 1 || vr.alt.size() > 1)
			continue;
		if ( _gi.size() != 0 )
		{
			Interval snpi(vr.chr, vr.pos-1, vr.pos);
			list<Interval> l = _gi.intersect(snpi);
			
			if ( l.size() == 0 )
				continue;
		}
		
		
		
		////////////////////////////////////////////////////////////////////////
		// Set SNP variables
		////////////////////////////////////////////////////////////////////////
		bool phased = false;
		string f = vr.samples.at(0).value.at(0);
		size_t split;
		int a1;
		int a2;
		bool het = false;
		int refID = _bamfs.GetReferenceID(vr.chr);
		int alleleCounts1Plus = 0;
		int alleleCounts1Minus = 0;
		int alleleCounts2Plus = 0;
		int alleleCounts2Minus = 0;
		int alleleCountsOtherPlus = 0;
		int alleleCountsOtherMinus = 0;
		
		if ( f.find("|") != string::npos)
		{
			split = f.find_first_of('|');
			phased = true;
		}
		else
			split = f.find_first_of('/');
		
		if (split == string::npos)
		{
			cout << "Incorrect format for sample field, no | or / delimiter for genotype." << endl;
			return false;
		}
		
		// Missing genotype. Output 0.
		if ( f.substr(0,split) == "." || f.substr(split+1,  f.find_first_of(':') - split - 1) == ".")
		{
			_ofs << vr.chr << '\t' << vr.pos << '\t' << vr.id << '\t' << vr.ref << '\t' << vr.alt << '\t' << f << '\t' <<
			alleleCounts1Plus << '\t' << alleleCounts1Minus << '\t' << alleleCounts2Plus << '\t' << alleleCounts2Minus << '\t' <<
			alleleCountsOtherPlus << '\t' << alleleCountsOtherMinus << endl;
			
			continue;
		}
		
		a1 = round(atof(f.substr(0,split).c_str()));
		a2 = round(atof(f.substr(split+1,  f.find_first_of(':') - split - 1).c_str()));
		
		if (a1 != a2)
			het = true;
		
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		// Loop through reads in snp location
		////////////////////////////////////////////////////////////////////////
		if ( ! _bamfs.SetRegion(refID, vr.pos-1, refID, vr.pos) ) // 0 based, open ended
		{
			cout << "Could not set region in bam file for " << vr.chr << ":" << vr.pos << endl;
			return false;
		}
		BamTools::BamAlignment al;
		while (_bamfs.GetNextAlignment(al))
		{
			////////////////////////////////////////////////////////////////////////
			// See if read passes parameter filtering
			////////////////////////////////////////////////////////////////////////
			if ( _useExonic && (!_readOverlap(al)) )
				continue;
			

			
			////////////////////////////////////////////////////////////////////////
			// Get SNP index and nucleotide in read
			////////////////////////////////////////////////////////////////////////
			int snp_index = _getLocIndex(al, vr.chr, vr.pos-1); // NEED TO PASS 0-BASED LOCATION TO FUNCTION
			if (snp_index < 0 || snp_index >= al.AlignedBases.length() )
			{
				cout << "Could not find SNP location in aligned bases from BamAlignment" << endl;
				return false;
			}
			char alignedSnpBase = al.AlignedBases.at(snp_index);
			
			////////////////////////////////////////////////////////////////////////
			// Check 5' and 3' cutoff
			////////////////////////////////////////////////////////////////////////
			if ( al.IsReverseStrand() )
			{
				int fivePrime = al.AlignedBases.length() - snp_index;
				int threePrime = snp_index;
				if ( fivePrime <= _5primecutoff )
					continue;
				if (threePrime <= _3primecutoff )
					continue;
			}
			else
			{
				int fivePrime = snp_index;
				int threePrime = al.AlignedBases.length() - snp_index;
				if ( fivePrime <= _5primecutoff )
					continue;
				if (threePrime <= _3primecutoff )
					continue;
			}
			
			
			// Alignment is spliced or deleted at the SNP location.
			if ( alignedSnpBase == 'N' || alignedSnpBase == '-' )
				continue;
			// If not A G C or T, there is something wrong
			if ( alignedSnpBase != 'A' && alignedSnpBase != 'C' && alignedSnpBase != 'G' && alignedSnpBase != 'T' )
			{
				cout << "Alignment doesn't have A, C, G, or T: " << alignedSnpBase << " in read " << al.Name << endl;
				return false;
			}
			
			
			////////////////////////////////////////////////////////////////////////
			// Determine which allele and strand this read is from
			////////////////////////////////////////////////////////////////////////
			bool refBase;
			bool altBase;
			bool otherBase;
			bool reverseStrand = al.IsReverseStrand();
			
			if ( alignedSnpBase == vr.ref.at(0) )
			{
				refBase = true;
				altBase = false;
				otherBase = false;
			}
			else if ( alignedSnpBase == vr.alt.at(0) )
			{
				refBase = false;
				altBase = true;
				otherBase = false;
			}
			else
			{
				refBase = false;
				altBase = false;
				otherBase = true;
			}
			
			
			////////////////////////////////////////////////////////////////////////
			// Read must have less than or equal to max mismatches
			////////////////////////////////////////////////////////////////////////
			unsigned int nMisMatches = UINT_MAX; // NOT SURE WHETHER TO MAKE THIS INT OR STRING??
			if ( al.GetTag("nM", nMisMatches) && refBase && nMisMatches > _refMaxMismatch )
				continue;
			else if ( al.GetTag("NM", nMisMatches) && refBase && nMisMatches > _altMaxMisMatch )
				continue;
			else if ( al.GetTag("nM", nMisMatches) && altBase && nMisMatches > _altMaxMisMatch )
				continue;
			else if ( al.GetTag("NM", nMisMatches) && altBase && nMisMatches > _altMaxMisMatch )
				continue;
			
			if (nMisMatches == UINT_MAX)
			{
				cout << "Could not determine number of mismatches in BamAlignment." << endl;
				return false;
			}
			
			
			
			
			
			////////////////////////////////////////////////////////////////////////
			// Assign read to appropriate allele/haplotype, and strand
			////////////////////////////////////////////////////////////////////////
			// If phased and heterozygous, assign read to phased haplotype
			if (het && phased && refBase && a1 == 0 && !reverseStrand)
				alleleCounts1Plus++;
			else if (het && phased && refBase && a2 == 0 && !reverseStrand)
				alleleCounts2Plus++;
			else if (het && phased && refBase && a1 == 0 && reverseStrand)
				alleleCounts1Minus++;
			else if (het && phased && refBase && a2 == 0 && reverseStrand)
				alleleCounts2Minus++;
			else if (het && phased && altBase && a1 == 1 && !reverseStrand)
				alleleCounts1Plus++;
			else if (het && phased && altBase && a2 == 1 && !reverseStrand)
				alleleCounts2Plus++;
			else if (het && phased && altBase && a1 == 1 && reverseStrand)
				alleleCounts1Minus++;
			else if (het && phased && altBase && a2 == 1 && reverseStrand)
				alleleCounts2Minus++;
			// Unphased genotype, or homozygous at this allele. In this case, allele 1 is ref. allele.
			else if ( refBase && !reverseStrand)
				alleleCounts1Plus++;
			else if ( refBase && reverseStrand)
				alleleCounts1Minus++;
			else if ( altBase && !reverseStrand)
				alleleCounts2Plus++;
			else if ( altBase && reverseStrand)
				alleleCounts2Minus++;
			// Another base besides reference or alternate. Obviously can't phase this so no allele1 or allele2
			else if (otherBase && !reverseStrand)
				alleleCountsOtherPlus++;
			else if (otherBase && reverseStrand)
				alleleCountsOtherMinus++;
			else
			{
				cout << "Something is missing in assigning this read to an allele: " << al.Name << endl;
				return false;
			}
		} // End for each read
		
		
		
		////////////////////////////////////////////////////////////////////////
		// Write information to file
		////////////////////////////////////////////////////////////////////////
		_ofs << vr.chr << '\t' << vr.pos << '\t' << vr.id << '\t' << vr.ref << '\t' << vr.alt << '\t' << f << '\t' <<
		alleleCounts1Plus << '\t' << alleleCounts1Minus << '\t' << alleleCounts2Plus << '\t' << alleleCounts2Minus << '\t' <<
		alleleCountsOtherPlus << '\t' << alleleCountsOtherMinus << endl;
		
	} // End for each SNP
	return true;
}



// TODO: Update so readOverlap can use P, =, and X, in CIGAR string
bool ASQuantifier::ASQuantifier::_readOverlap(const BamTools::BamAlignment& al)
{
	if (_gi.size() == 0)
		return false;
	
	string chr = _bamfs.GetReferenceData().at(al.RefID).RefName;
	int alignStart = al.Position;
	
	for (unsigned int i = 0; i < al.CigarData.size(); i++)
	{
		char t = al.CigarData.at(i).Type;
		unsigned int l = al.CigarData.at(i).Length;
		
		// hard or soft clip is not part of alignment so pass
		if (t == 'H' || t == 'S' || t == 'I')
			continue;
		
		if (t == 'D' || t == 'N')
		{
			alignStart += l;
			continue;
		}
		
		if (t == 'M')
		{
			int alignEnd = alignStart + l;
			Interval subalign(chr, alignStart, alignEnd);
			list<Interval> lint = _gi.intersect(subalign);
			if ( lint.size() == 0 )
				return false;
			alignStart = alignEnd;
		}
	}
	return true;
}



// Expects 0-based location loc
int ASQuantifier::ASQuantifier::_getLocIndex(const BamTools::BamAlignment& al, string chr, int loc)
{
	string alignChr = _bamfs.GetReferenceData().at(al.RefID).RefName;
	int alignStart = al.Position;
	
	if (alignChr != chr)
		return -1;
	
	if ( loc < alignStart )
		return -1;
	
	int alignedIndex = 0;
	
	for (unsigned int i = 0; i < al.CigarData.size(); i++)
	{
		char t = al.CigarData.at(i).Type;
		unsigned int l = al.CigarData.at(i).Length;
		
		// hard or soft clip, or insertion is not part of alignment so pass
		if (t == 'H' || t == 'S')
			continue;
		else if ( t == 'I' )
			alignedIndex += l;
		else if (t == 'D' || t == 'N' || t == 'M')
		{
			int alignEnd = alignStart + l;
			if ( loc >= alignStart && loc < alignEnd )
				return ( alignedIndex += ( loc - alignStart ) );
			
			alignStart += l;
			alignedIndex += l;
		}
		else
		{
			cout << "Unrecognized CIGAR string: " << t << endl;
			exit(1);
		}
	}
	
	return -1;
	
}