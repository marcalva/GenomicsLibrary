// SnpStruct.h

#ifndef __SNPSTRUCT__
#define __SNPSTRUCT__

#include <string>

struct SnpBase {
public:
	std::string chr;
	unsigned int pos;
	std::string name;
};

struct SnpLd : public SnpBase {
public:
	double r2;
	double Dprime;
};

#endif // __SNPSTRUCT__