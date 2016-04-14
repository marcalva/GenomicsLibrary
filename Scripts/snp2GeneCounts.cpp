// SnpToGene Counts.
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>
#include <cassert>
#include <sstream>
#include "Interval.h"
#include "GenomicInterval.h"

using namespace std;

const int LINESIZE = 65536;

vector<string> parseString(string line)
{
	stringstream ss(line);
	vector<string> outVector;
	string outString;
	while ( ss >> outString )
		outVector.push_back(outString);
	return outVector;
}

int main(int argc, char* argv[])
{
	
	std::string ac_filename = argv[1];
	ifstream ac_fs(ac_filename.c_str(), std::ifstream::in);
	
	// BED file format
	std::string interval_filename = argv[2];
	ifstream interval_fs(interval_filename.c_str(), std::ifstream::in);
	
	// Store allele counts into matrix
	map< string, vector<double> > snp_counts;
	map< string, Interval> locations;
	vector<string> snp_names;
	unsigned int nSamples;
	bool s = false;
	char line[LINESIZE];
	while (ac_fs.getline(line, LINESIZE, '\n'))
	{
		vector<string> line_vector = parseString(string(line));
		string chr = line_vector.at(0);
		int start = atoi(line_vector.at(1).c_str());
		int end = start;
		start--;
		string name = line_vector.at(2);
		
		vector<double> allele_count;
		for (unsigned int i = 3; i < line_vector.size(); i++)
			allele_count.push_back(atof(line_vector.at(i).c_str()));
		
		if (!s)
			nSamples = allele_count.size();
		else
			assert(nSamples == allele_count.size());
		
		Interval i(chr, start, end, name);
		
		snp_counts.insert(pair<string,vector< double > >(name, allele_count));
		locations.insert(pair<string,Interval>(name,i));
		snp_names.push_back(name);
	}
	
	ac_fs.close();
	
	// Store intervals into GenomicInterval and create gene counts
	GenomicInterval intervals;
	map<string, vector<double> > gene_counts;
	while (interval_fs.getline(line, LINESIZE, '\n'))
	{
		vector<string> line_vector = parseString(string(line));
		string chr = line_vector.at(0);
		int start = atoi(line_vector.at(1).c_str());
		int end = atoi(line_vector.at(2).c_str());
		string genename = line_vector.at(3);
		
		Interval i(chr, start, end, genename);
		
		intervals.addInterval(i);
		
		// Create gene counts dictionary
		vector<double> c(nSamples, 0);
		gene_counts.insert( pair<string, vector<double> >(genename, c) );
		
	}
	
	interval_fs.close();
	
	for (unsigned int i = 0; i < snp_names.size(); i++)
	{
		string snp_name = snp_names.at(i);
		Interval i = locations[snp_name];
		list<Interval> intersections = intervals.intersect(i);
		
		if ( intersections.size() > 1 )
		{
			cout << "Size " << intersections.size() << endl;
			for ( list<Interval>::const_iterator lit = intersections.begin() ; lit != intersections.end(); ++lit )
			{
				cout << lit->Chr() << " " << lit->Start() << " " << lit->End() << " " << lit->Name() << endl;
			}
		}
		/**
		if (intersections.size() > 1)
		{
			cout << "BED file has overlapping intervals" << endl;
			for ()
				cout << intersections.
			exit(1);
		}
		*/
		
		if (intersections.size() == 0)
		{
			continue;
		}
		
		
		Interval overlap_interval = intersections.front();
		
		string overlap_gene = overlap_interval.Name();
		
		for (unsigned int j = 0; j < gene_counts[overlap_gene].size(); j++)
		{
			gene_counts[overlap_gene].at(j) += snp_counts[snp_name].at(j);
		}
	}
	

	// Write gene level counts to file
	
	string ofsfilename  = argv[3];
	ofstream ofs(ofsfilename.c_str(), std::ofstream::out);
	for ( map<string, vector<double> >::const_iterator gcit = gene_counts.begin(); gcit != gene_counts.end(); gcit++ )
	{
		string g = gcit->first;
		ofs << g;
		
		for (unsigned int i = 0; i < nSamples; i++)
			ofs << ',' << gcit->second.at(i);
		ofs << endl;
	}
	
	ofs.close();
	
}