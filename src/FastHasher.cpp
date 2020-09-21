#include "FastHasher.h"

#include <vector>

FastHasher::FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash) :
fwHash(fwHash),
bwHash(bwHash),
kmerSize(kmerSize % 64)
{
	precalcRots();
}

FastHasher::FastHasher(size_t kmerSize) :
fwHash(0),
bwHash(0),
kmerSize(kmerSize % 64)
{
	precalcRots();
}

char FastHasher::complement(char c) const
{
	static std::vector<char> comp { 0, 4, 3, 2, 1 };
	return comp[c];
}

void FastHasher::precalcRots()
{
	for (int i = 0; i < 5; i++)
	{
		fwAdd[i] = charHashes[i];
	}
	for (int i = 0; i < 5; i++)
	{
		fwRemove[i] = rotlk(charHashes[i]);
	}
	for (int i = 0; i < 5; i++)
	{
		bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)]);
	}
	for (int i = 0; i < 5; i++)
	{
		bwRemove[i] = rotrone(charHashes[(int)complement(i)]);
	}
}
