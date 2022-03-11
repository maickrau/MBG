#include <cassert>
#include "DumbSelect.h"

DumbSelect::DumbSelect(size_t bitvectorLength) :
bitvector(bitvectorLength)
{
}

size_t DumbSelect::selectOne(size_t count) const
{
	assert(count < countOnes());
	if (count == 0) return 0;
	size_t start = 0;
	size_t end = bitvector.size();
	while (end > start+1)
	{
		size_t mid = (end - start) / 2 + start;
		size_t midCount = bitvector.getRank(mid) + (bitvector.get(mid) ? 1 : 0);
		if (midCount < count)
		{
			start = mid;
		}
		else
		{
			end = mid;
		}
	}
	assert(end == start+1);
	assert(end < bitvector.size());
	assert(bitvector.get(end));
	assert(bitvector.getRank(end) + 1 == count);
	return end;
}

size_t DumbSelect::countOnes() const
{
	return numOnes + 1;
}

void DumbSelect::build()
{
	bitvector.buildRanks();
	numOnes = bitvector.getRank(bitvector.size()-1) + (bitvector.get(bitvector.size()-1) ? 1 : 0);
}
