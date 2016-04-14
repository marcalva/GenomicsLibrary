/** Interval.h */

#ifndef _INTERVAL_
#define _INTERVAL_
#include <string>
#include <iostream>
#include <cstdlib>

class Interval
{
private:
	std::string _chr;
	int _start;
	int _end;
	std::string _name;
	std::string _score;
	std::string _strand;
public:
	Interval()
	{}
	
	// Construct interval with known coordinates
	Interval(const std::string& c, int st, int en)
	{
		if (st < 0 || en < 0)
		{
			std::cout << "Intervals can only have 0 or positive integers" << std::endl;
			std::exit(1);
		}
		_chr = c;
		_start = st;
		_end = en;
	}
	
	Interval(const std::string& c, int st, int en, std::string name)
	{
		if (st < 0 || en < 0)
		{
			std::cout << "Intervals can only have 0 or positive integers" << std::endl;
			std::exit(1);
		}
		_chr = c;
		_start = st;
		_end = en;
		_name = name;
	}
	
	Interval(const std::string& c, int st, int en, std::string name, std::string score, std::string strand)
	{
		if (st < 0 || en < 0)
		{
			std::cout << "Intervals can only have 0 or positive integers" << std::endl;
			std::exit(1);
		}
		_chr = c;
		_start = st;
		_end = en;
		_name = name;
		_strand = strand;
	}
	
	void setName(const std::string& name) {_name = name;}
	void setScore(const std::string& score) {_score = score;}
	void setStrand(const std::string& strand) {_strand = strand;}
	
	// Accessor functions
	std::string Chr() const {return _chr;}
	int Start() const {return _start;}
	int End() const {return _end;}
	std::string Name() const {return _name;}
	std::string Score() const {return _score;}
	std::string Strand() const {return _strand;}
	
	/**
	Important:
	Intervals are compared as follows:
		Chromosomes are compared with string comparison, not numeric comparison.
		Only the start coordinate is used. Intervals with the same chromosome and start, but different end will be equal.
		Start coordinate is compared as an integer.
	*/
	
	bool operator<(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		if (leftChr < rightChr)
			return true;
		else if (leftChr > rightChr)
			return false;
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart < rightStart)
				return true;
			else // leftStart >= rightStart
				return false;
		}
	}
	
	bool operator<=(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		if (leftChr < rightChr)
			return true;
		else if (leftChr > rightChr)
			return false;
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart <= rightStart)
				return true;
			else // leftStart > rightStart
				return false;
		}
	}
	
	bool operator==(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		if (leftChr != rightChr)
			return false;
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart == rightStart)
				return true;
			else // leftStart != rightStart
				return false;
		}
	}
	
	bool operator!=(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		
		if (leftChr != rightChr)
			return true;
		
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart != rightStart)
				return true;
			else // leftStart == rightStart
				return false;
		}
	}
	
	bool operator>(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		
		if (leftChr > rightChr)
			return true;
		else if (leftChr < rightChr)
			return false;
		
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart > rightStart)
				return true;
			else // leftStart <= rightStart
				return false;
		}
	}
	
	bool operator>=(const Interval& rhs) const
	{
		std::string leftChr = Chr();
		std::string rightChr = rhs.Chr();
		
		if (leftChr > rightChr)
			return true;
		else if (leftChr < rightChr)
			return false;
		
		else // leftChr == rightChr
		{
			int leftStart = Start();
			int rightStart = rhs.Start();
			
			if (leftStart >= rightStart)
				return true;
			else // leftStart < rightStart
				return false;
		}
	}
};

#endif // Interval.h