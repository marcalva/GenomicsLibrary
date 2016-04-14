// GenomicInterval.h

#include <cstdlib>
#include <string>
#include <iostream>
#include <set>
#include <list>
#include "Interval.h"

class GenomicInterval{
private:
	std::multiset<Interval> _intervals;
	
public:
	// Add interval to _intervals BST
	void addInterval(const Interval& i);
	
	unsigned int size();
	
	// Find intersecting intervals from interval, return a list.
	std::list<Interval> intersect(const Interval& i) const;
	
	// Find intersecting intervals from another set of intervals:
	// TODO Need to figure out how to implement this
	// std::list <Interval> intersect(const GenomicInterval& gi) const;
	
};