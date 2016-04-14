//  GenomicIntervalTabix.cpp

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include "Interval.h"
#include "GenomicIntervalTabix.h"


using namespace std;

vector<string> parse_string(const string in, const char delim, vector<string>& out)
{
	out.clear();
	stringstream ss(in);
	string token;
	
	while (getline(ss, token, delim))
		out.push_back(token);
	
	return out;
}

GenomicIntervalTabix::GenomicIntervalTabix(const std::string& tFileName) : _t(tFileName) {}

std::list<Interval> GenomicIntervalTabix::intersect(const Interval& i) const
{
	std::list<Interval> intersections;
	
	string region = i.Chr() + ":" + to_string(i.Start()) + "-" + to_string(i.End());
	string tLine;
	vector<string> tVector;
	if (_t.setRegion(region))
	{
		while (_t.getNextLine(tLine))
		{
			parse_string(tLine, '\t', tVector);
			chr = tVector.at(0);
			start = atoi(tVector.at(1));
			end = atoi(tVector.at(2));
			name = tVector.at(3);
			Interval thisI(chr, start, end, name);
			if (chr == i.Chr() && i.Start() >= start && i.End() <= end)
				intersections.push_back(i);
		}
	}
	return intersections;
}