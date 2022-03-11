#ifndef DumbSelect_h
#define DumbSelect_h

#include "RankBitvector.h"

class DumbSelect
{
public:
	DumbSelect(size_t bitvectorLength);
	size_t selectOne(size_t count) const;
	size_t countOnes() const;
	void build();
	RankBitvector bitvector;
private:
	size_t numOnes;
};

#endif
