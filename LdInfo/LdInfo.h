// LdInfo.h
#ifndef __LD__
#define __LD__

#include <string>
#include <tabix.hpp>
#include <cstdlib>
#include <list>
#include <unordered_set>
#include "SnpStruct.h"
#include "VcfRecord.h"

struct LdInfo
{
public:
	double r2;
	double Dprime;
};

class LinkageDisequilibrium
{
public:
	// Constructor. Note the file name for the tabix indexed VCF is needed
	LinkageDisequilibrium(std::string vcfFileName);
	
	// Given a SNP, return a list of SNPs and their LD information.
	void GetProxies(const SnpBase& s,
		std::list<SnpLd>* ld,
		const float& minR2,
		const float& minDp,
		const int& windowSize);

	// Given 2 SNPs, return a list of SNPs and their LD information.
	void GetProxies(const SnpBase& s1,
		const SnpBase& s2,
		std::list<SnpLd>* ld,
		const float& minR2,
		const float& minDp,
		const int& windowSize);
	
	// Given a SNP, return a list of SNPs and LD information, restricted to snp list restrL
	void GetProxies(const SnpBase& s, 
		std::list<SnpLd>* ld,
		const std::unordered_set<std::string>& restrL,
		const float& minR2,
		const float& minDp,
		const int& windowSize);
	
	// Check if a chromosome is present in tabix file
	bool hasChr(std::string c);
	
private:
	Tabix _tVcf;
};

#endif // __LD__
