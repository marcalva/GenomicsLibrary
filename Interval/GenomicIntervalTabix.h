//  GenomicIntervalTabix.h

#ifndef GenomicIntervalTabix_h
#define GenomicIntervalTabix_h
#include <cstdlib>
#include <string>
#include "Interval.h"
#include <tabix.hpp>

class GenomicIntervalTabix{
private:
	Tabix _t;
public:
	// Create tabix file. Only way is through constructor because we cannot have intervals without initializing the tabix file
	GenomicIntervalTabix(const std::string& tFileName);
	
	// Find intersecting intervals from interval, return a list.
	std::list<Interval> intersect(const Interval& i) const;
	
	// Find intersecting intervals from another set of intervals:
	// TODO Need to figure out how to implement this
	// std::list <Interval> intersect(const GenomicInterval& gi) const;
	
};

#endif /* GenomicIntervalTabix_h */
