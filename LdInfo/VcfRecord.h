// VcfRecord.h
#ifndef __VCFRECORD__
#define __VCFRECORD__

#include <string>
#include <vector>

struct infoField{
	// Each component is separated by a semicolon in the VCF format
	std::vector<std::string> key;
	std::vector<std::string> value;
	std::vector<std::string> flags;
};

struct formatField{
	// Each component is separated by a colon in the VCF format
	std::vector<std::string> value;
};

struct sampleField{
	// Each component is separated by a colon in the VCF format
	// Must have the same length as formatField.
	std::vector<std::string> value;
};

struct vcfHeader
{
	// Meta lines that start with ## 
	std::vector<std::string> meta;
	
	// Header lne beginning with #CHROM
	std::string head;
};

struct vcfRecord{
	std::string chr;
	int pos;
	std::string id;
	std::string ref;
	std::string alt;
	std::string qual;
	std::string filter;
	infoField info;
	formatField format;
	std::vector<sampleField> samples;
	void parseLine(std::string s); // Store VCF line from a string
	void clear();
};

#endif // __VCFRECORD__
