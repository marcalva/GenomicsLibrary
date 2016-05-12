// LdInfo.cpp
#include <string>
#include <cassert>
#include <list>
#include <vector>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <algorithm>
#include "LdInfo.h"
#include "VcfRecord.h"
#include "ParseText.h"

using namespace std;

// Calculate LD function
LdInfo calcLd(const vcfRecord& r1, const vcfRecord& r2)
{
	string gfn = "GT";
	LdInfo l;
	// Need to incorporate unphased LD calculation
	// bool phased;
	size_t nSamples = r1.samples.size();
	size_t nHaps = nSamples * 2;
	float p1, p2, q1, q2; // allele frequency for r1 and r2 snp respectively
	float a1b1 = 0, a1b2 = 0, a2b1 = 0, a2b2 = 0; // Haplotype frequencies
	
	if ( r1.samples.size() != r2.samples.size() )
	{
		cerr << "These two SNPs hve differing numbers of samples:" << endl;
		cerr << r1.id << " " << r2.id << endl;
		exit (EXIT_FAILURE);
	}
	
	size_t r1index = 0, r2index = 0;
	while ( r1index < r1.format.value.size() )
	{
		if (r1.format.value.at(r1index) == gfn)
			break;
		else
			r1index++;
	}
	while ( r2index < r2.format.value.size() )
	{
		if (r2.format.value.at(r2index) == gfn)
			break;
		r2index++;
	}
	if (r1index == r1.format.value.size() || r2index == r2.format.value.size())
	{
		cerr << "Could not find genotype column " << gfn << " in the vcf file for SNPs " << r1.id << " " << r2.id << endl;
		exit (EXIT_FAILURE);
	}
	
	// TODO : Incorporate ability to get LD for unphased genotypes.
	for (size_t i = 0; i < nSamples; i++)
	{
		vector<string> s1g, s2g;
		parse(r1.samples.at(i).value.at(r1index), '|', s1g);
		parse(r2.samples.at(i).value.at(r2index), '|', s2g);
		if (s1g.size() != 2 || s2g.size() != 2){
			cerr << "Genotypes not recognized at SNPs "<< r1.id << " " << r2.id << endl;
			exit (EXIT_FAILURE);
		}
		string h1 = s1g.at(0) + s2g.at(0);
		string h2 = s1g.at(1) + s2g.at(1);
		vector<string> h;
		h.push_back(h1);
		h.push_back(h2);
		
		for (size_t j = 0; j < h.size(); j++)
		{
			if (h.at(j) == "00")
				a1b1++;
			else if (h.at(j) == "01")
				a1b2++;
			else if (h.at(j) == "10")
				a2b1++;
			else if (h.at(j) == "11")
				a2b2++;
			else{
				cerr << "Genotypes not recognized at SNPs "<< r1.id << " " << r2.id << endl;
				exit (EXIT_FAILURE);
			}
		}
	} // end for
	
	a1b1 /= nHaps;
	a1b2 /= nHaps;
	a2b1 /= nHaps;
	a2b2 /= nHaps;
	
	p1 = a1b1 + a1b2;
	p2 = a2b1 + a2b2;
	q1 = a1b1 + a2b1;
	q2 = a1b2 + a2b2;
	
	float d = (a1b1 * a2b2) - (a1b2 * a2b1);
	float dMax;
	if (d < 0){
		if ( (p1*q1) < (p2*q2))
			dMax = p1 * q1 * -1;
		else
			dMax = p2 * q2 * -1;
	}
	else{
		if ( (p1*q2) < (p2*q1))
			dMax = p1 * q2;
		else
			dMax = p2 * q1;
	}
	
	l.Dprime = d / dMax;
	l.r2 = (d * d) / (p1 * p2 * q1 * q2);
	
	return l;
}

// Constructor
LinkageDisequilibrium::LinkageDisequilibrium(std::string vcfFileName) : _tVcf(vcfFileName) {}

// TODO : Put error controls if SNPs or regions could not be set.
void LinkageDisequilibrium::GetProxies(const SnpBase& s, std::list<SnpLd>* ld,  const float& minR2, const int& windowSize)
{
	ld->clear();
	
	int beginRegion = (s.pos - 1) - windowSize;
	int endRegion = s.pos + windowSize;
	string snpRegion = s.chr + ":" + to_string(s.pos - 1) + "-" + to_string(s.pos);
	string windowRegion = s.chr + ":" + to_string(beginRegion) + "-" + to_string(endRegion);
	string tLine;
	vcfRecord refSnp;
	vcfRecord proxySnp;
	vector<vcfRecord> tVector;
	
	if (_tVcf.setRegion(snpRegion))
	{
		while (_tVcf.getNextLine(tLine))
		{
			refSnp.parseLine(tLine);
			if (refSnp.id == s.name)
				break;
		}
		if (refSnp.id != s.name){
			cerr << "Could not find input SNP " << s.name << "." << endl;
			exit (EXIT_FAILURE);
		}
	}
	else
	{
		cerr << "Could not set region from tabix-indexed VCF file. Check chromosomes and positions." << endl;
		exit (EXIT_FAILURE);
	}
	
	if (_tVcf.setRegion(windowRegion))
	{
		while (_tVcf.getNextLine(tLine))
		{
			proxySnp.parseLine(tLine);
			LdInfo l;
			l = calcLd(refSnp, proxySnp);
			if (l.r2 < minR2)
				continue;
			SnpLd proxyLd;
			proxyLd.chr = proxySnp.chr;
			proxyLd.pos = proxySnp.pos;
			proxyLd.name = proxySnp.id;
			proxyLd.r2 = l.r2;
			proxyLd.Dprime = l.Dprime;
			ld->push_back(proxyLd);
		}
	}
}

void LinkageDisequilibrium::GetProxies(const SnpBase& s, 
	std::list<SnpLd>* ld,
	const std::unordered_set<std::string>& restrL,
	const float& minR2,
	const int& windowSize)
{
	ld->clear();
	
	int beginRegion = (s.pos - 1) - windowSize;
	int endRegion = s.pos + windowSize;
	string snpRegion = s.chr + ":" + to_string(s.pos - 1) + "-" + to_string(s.pos);
	string windowRegion = s.chr + ":" + to_string(beginRegion) + "-" + to_string(endRegion);
	string tLine;
	vcfRecord refSnp;
	vcfRecord proxySnp;
	vector<vcfRecord> tVector;
	
	if (_tVcf.setRegion(snpRegion))
	{
		while (_tVcf.getNextLine(tLine))
		{
			refSnp.parseLine(tLine);
			if (refSnp.id == s.name)
				break;
		}
		if (refSnp.id != s.name){
			cerr << "Could not find SNP " << s.name << "." << endl;
			exit (EXIT_FAILURE);
		}
	}
	else
	{
		cerr << "Could not set region from tabix-indexed VCF file." << endl;
		exit (EXIT_FAILURE);
	}
	
	if (_tVcf.setRegion(windowRegion))
	{
		while (_tVcf.getNextLine(tLine))
		{
			proxySnp.parseLine(tLine);
			std::unordered_set<std::string>::const_iterator got = restrL.find (proxySnp.id);
			if (got != restrL.end())
			{
				LdInfo l;
				l = calcLd(refSnp, proxySnp);
				SnpLd proxyLd;
				proxyLd.chr = proxySnp.chr;
				proxyLd.pos = proxySnp.pos;
				proxyLd.name = proxySnp.id;
				proxyLd.r2 = l.r2;
				proxyLd.Dprime = l.Dprime;
				if (l.r2 >= minR2)
					ld->push_back(proxyLd);
				
			}
		}
	}
}

bool LinkageDisequilibrium::hasChr(std::string c)
{
	vector<string>::iterator sp;
	sp = find(_tVcf.chroms.begin(), _tVcf.chroms.end(), c);
	if (sp == _tVcf.chroms.end())
		return false;
	else
		return true;
}