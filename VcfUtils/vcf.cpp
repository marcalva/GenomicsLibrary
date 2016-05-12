#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <sstream>
#include "vcf.h"

using namespace std;
const int WORD_SIZE = 65536;

vcfFileInStream::vcfFileInStream()
{
}


vcfFileInStream::vcfFileInStream(std::string fn)
{
	ifs.open(fn.c_str()); 
}

bool vcfFileInStream::openVcf(std::string fn)
{
	ifs.open(fn.c_str());
	
	if ( ifs.fail() )
		return false;
	else
		return true;
}

vcfRecord vcfFileInStream::readRecord()
{
	vcfRecord r;
	ifs.peek();
	if (!ifs.good())
		return r;
	
	char line_array[100000];
	ifs.getline(line_array, 100000, '\n');
	
	std::string line = line_array;
	
	std::vector<std::string> line_vector;
	
	// MAKE SURE THAT VCF RECORD IS TAB-DELIMITED
	_split(line, '\t', line_vector);
	
	if (line_vector.size() < 10)
	{
		std::cout << "There are less than 10 columns in this VCF record. Please correct file format." << std::endl;
		std::cout << "Number of columns: " << line_vector.size() << endl;
		for (unsigned int i = 0 ; i < line_vector.size(); i++)
			cout << line_vector.at(i) << "\t";
		cout << endl;
		std::exit(1);
	}
	
	r.chr = line_vector.at(0);
	r.pos = std::stoi(line_vector.at(1));
	r.id = line_vector.at(2);
	r.ref = line_vector.at(3);
	r.alt = line_vector.at(4);
	r.qual = line_vector.at(5).at(0);
	r.filter = line_vector.at(6);
	
	std::vector<std::string> iField;
	_split(line_vector.at(7), ';', iField);
	
	if (! ( line_vector.at(7) == "" || line_vector.at(7) == "." ) ) // If the info field is not empty
	{
		for (unsigned long i = 0; i < iField.size(); i++)
		{
			std::vector<std::string> v;
			_split(iField.at(i), '=', v);
			
			// If this info field is a flag
			if ( v.size() == 1 )
			{
				r.info.key.push_back(v.at(0));
				r.info.value.push_back("");
			}
			else // If the info field is a key value pair
			{
				r.info.key.push_back(v.at(0));
				r.info.value.push_back(v.at(1));
			}
		}
	}
	
	_split(line_vector.at(8), ':', r.format.value);
	
	for (unsigned int i = 9; i < line_vector.size(); i++)
	{
		sampleField s;
		_split(line_vector.at(i), ':', s.value);
		r.samples.push_back(s);
	}

	return r;
	
}

vcfHeader vcfFileInStream::readHeader()
{
	vcfHeader h;
	
	if (!ifs.good())
		return h;
	
	// Go to beginning of file
	ifs.seekg(0, ifs.beg);
	
	char line_array[100000];
	std::string line;
	do
	{
		ifs.getline(line_array, 100000, '\n');
		
		line = line_array;
		if (line.substr(0, 2) == "##")
			h.meta.push_back(line);
		
	} while (line.substr(0, 2) == "##");
	
	assert(line.substr(0, 6) == "#CHROM");
	
	h.head = line;
	
	return h;
	
}

bool vcfFileInStream::atEnd()
{
	ifs.peek();
	return ifs.eof();
}

void vcfFileInStream::close()
{
	ifs.close();
}

bool vcfFileInStream::good()
{
	return ifs.good();
};

vcfFileOutStream::vcfFileOutStream(std::string fn)
{
	ofs.open(fn.c_str());
}

void vcfFileOutStream::writeRecord(vcfRecord v)
{
	if (!ofs.good())
	{
		std::cout << "Attempting to write to non-good output file." << std::endl;
		std::exit(1);
	}
	
	long long unsigned int p = v.pos;
	std::string q = string(1, v.qual);
	std::string output_line = v.chr + "\t" + \
		std::to_string(p) + "\t" + \
			v.id + "\t" + v.ref + "\t" + v.alt + "\t" + \
				q + "\t" + v.filter + "\t";

	if (v.info.key.size() > 0)
	{
		for (unsigned int i = 0; i < v.info.key.size(); i++)
		{
			if ( v.info.value.at(i) == "" ) // If key is a flag
				output_line += v.info.key.at(i);
			else // If key-value pair
				output_line += v.info.key.at(i) + "=" + output_line += v.info.value.at(i);
			
			// Don't place ; delimiter if at the last value.
			if (i < v.info.key.size() - 1)
				output_line += ";";
		}
	}
	else
		output_line += ".";
		
	output_line += "\t";
	
	for (unsigned int i = 0; i < v.format.value.size(); i++)
	{
		output_line += v.format.value.at(i);
		
		if (i < v.format.value.size() - 1)
			output_line += ":";
	}
	
	for (unsigned int i = 0; i < v.samples.size(); i++)
	{
		output_line += "\t";
		
		for (unsigned int j = 0; j < v.samples.at(i).value.size(); j++)
		{
			output_line += v.samples.at(i).value.at(j);
			
			if (j < v.samples.at(i).value.size() - 1)
				output_line += ":";
		}
	}
		
	output_line += "\n";
	
	ofs << output_line;
	
}

void vcfFileOutStream::writeHeader(vcfHeader h)
{
	if (!ofs.good())
	{
		std::cout << "Attempting to write to non-good output file." << std::endl;
		std::exit(1);
	}
	
	ofs.seekp(0, ofs.beg);
	
	for (int i = 0; i < h.meta.size(); i++)
	{
		ofs << h.meta.at(i) << "\n";
	}
	
	ofs << h.head << "\n";
}

void vcfFileOutStream::close()
{
	ofs.close();
}

bool vcfFileOutStream::good()
{
	return ofs.good();
}

std::vector<std::string>& vcfFileInStream::_split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}
