// ParseText.h
#ifndef __PARSETEXT__
#define __PARSETEXT__
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

inline void parse(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}
#endif // __PARSETEXT__