// VcfRecord.cpp
#include "VcfRecord.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include "ParseText.h"

using namespace std;

void vcfRecord::parseLine(std::string s)
{
	clear();
	std::vector<std::string> line_vector;
	char delimiter = '\t';
	parse(s, delimiter, line_vector);
	if ( line_vector.size() < 4 )
	{
		cerr << "VCF file must have at least 3 columns for chr, position, and id." << endl;
		cerr << "The file must be tab-delimited." << endl;
		exit (EXIT_FAILURE);
	}
	
	chr = line_vector.at(0);
	pos = std::stoi(line_vector.at(1));
	id = line_vector.at(2);
	ref = line_vector.at(3);
	alt = line_vector.at(4);
	qual = line_vector.at(5);
	filter = line_vector.at(6);
	
	vector<string> iField;
	parse(line_vector.at(7), ';', iField);
	
	for ( unsigned long i = 0; i < iField.size(); i++ )
	{
		size_t keyValue = iField.at(i).find('=');
		if ( keyValue != string::npos)
		{
			vector<string> v;
			parse(iField.at(i), '=', v);
			info.key.push_back(v.at(0));
			info.value.push_back(v.at(1));
		}
		else
			info.flags.push_back(iField.at(i));
	}
	
	parse(line_vector.at(8), ':', format.value);
	
	for (size_t i = 9; i < line_vector.size(); i++)
	{
		sampleField sf;
		parse(line_vector.at(i), ':', sf.value);
		samples.push_back(sf);
	}
}

void vcfRecord::clear()
{
	chr.clear();
	pos = -1;
	id.clear();
	ref.clear();
	alt.clear();
	qual.clear();
	filter.clear();
	info.key.clear();
	info.value.clear();
	info.flags.clear();
	format.value.clear();
	samples.clear();
}
