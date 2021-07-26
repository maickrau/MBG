#include <cassert>
#include "RankBitvector.h"

int popcount(uint64_t x)
{
	//https://gcc.gnu.org/onlinedocs/gcc-4.8.4/gcc/X86-Built-in-Functions.html
	// return __builtin_popcountll(x);
	//for some reason __builtin_popcount takes 21 instructions so call assembly directly
	__asm__("popcnt %0, %0" : "+r" (x));
	return x;
}

RankBitvector::RankBitvector(size_t size) :
	ranksBuilt(false)
{
	bits.resize((size + BitsPerChunk - 1) / BitsPerChunk, 0);
	realSize = size;
}

size_t RankBitvector::size() const
{
	return realSize;
}

void RankBitvector::set(size_t index, bool value)
{
	assert(!ranksBuilt);
	size_t chunk = index / BitsPerChunk;
	size_t offset = index % BitsPerChunk;
	if (value)
	{
		bits[chunk] |= (uint64_t)1 << (uint64_t)offset;
	}
	else
	{
		bits[chunk] &= ~((uint64_t)1 << (uint64_t)offset);
	}
}

bool RankBitvector::get(size_t index) const
{
	size_t chunk = index / BitsPerChunk;
	size_t offset = index % BitsPerChunk;
	return ((bits[chunk] >> offset) & 1) == 1;
}

void RankBitvector::buildRanks()
{
	assert(!ranksBuilt);
	smallRanks.resize(bits.size());
	bigRanks.reserve((bits.size() + SmallRanksPerBig - 1) / SmallRanksPerBig);
	size_t runningSum = 0;
	size_t runningBigSum = 0;
	for (size_t i = 0; i < bits.size(); i++)
	{
		if (i % SmallRanksPerBig == 0)
		{
			runningBigSum += runningSum;
			bigRanks.push_back(runningBigSum);
			runningSum = 0;
		}
		smallRanks[i] = runningSum;
		runningSum += popcount(bits[i]);
	}
	assert(bigRanks.size() == (bits.size() + SmallRanksPerBig - 1) / SmallRanksPerBig);
	ranksBuilt = true;
	size_t runningRank = 0;
	for (size_t i = 0; i < size(); i++)
	{
		assert(getRank(i) == runningRank);
		if (get(i)) runningRank += 1;
	}
}

size_t RankBitvector::getRank(size_t index) const
{
	assert(ranksBuilt);
	size_t chunk = index / BitsPerChunk;
	size_t offset = index % BitsPerChunk;
	uint64_t mask = ((uint64_t)1 << (uint64_t)offset) - 1;
	size_t chunkSum = popcount(bits[chunk] & mask);
	size_t bigChunk = chunk / SmallRanksPerBig;
	return chunkSum + smallRanks[chunk] + bigRanks[bigChunk];
}
