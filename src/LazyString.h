#ifndef LazyString_h
#define LazyString_h

#include <string>
#include <string_view>
#include "TwobitString.h"

class LazyString
{
public:
	LazyString(TwobitView first, TwobitView second, size_t overlap);
	char operator[](size_t index) const;
	std::string_view view(size_t start, size_t kmerSize);
	size_t size() const;
private:
	TwobitView first;
	TwobitView second;
	const size_t overlap;
	std::string str;
};

#endif