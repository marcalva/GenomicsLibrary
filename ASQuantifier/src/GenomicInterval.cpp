// GenomicInterval.cpp
#include "GenomicInterval.h"
#include "Interval.h"

using namespace std;

void GenomicInterval::addInterval(const Interval& i)
{
	_intervals.insert(i);
}

std::list<Interval> GenomicInterval::intersect(const Interval& i) const
{
	std::list<Interval> intersections;
	
	std::multiset<Interval>::iterator it;
	
	Interval set_int(i.Chr(), 0, 1); // Go through all intervals in chromosome
	
	it = _intervals.lower_bound(set_int);
	
	while ( it != _intervals.end() )
	{
		if ( i.Chr() == it->Chr() && i.Start() >= it->Start() && i.End() <= it->End() )
			intersections.push_back(*it);
		else if ( it->Start() > i.End() )
			break; // _intervals are sorted
		++it;
	}
		
	return intersections;
}

unsigned int GenomicInterval::size()
{
	return _intervals.size();
}