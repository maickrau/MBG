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
	static char mapping[29] { 0, 4, 3, 2, 1, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 25, 22, 19, 28, 21, 18, 27, 24, 17, 26, 23, 20 };
	return mapping[(int)c];
}

void FastHasher::precalcRots()
{
	for (int i = 0; i < 29; i++)
	{
		fwAdd[i] = charHashes[i];
	}
	for (int i = 0; i < 29; i++)
	{
		fwRemove[i] = rotlk(charHashes[i]);
	}
	for (int i = 0; i < 29; i++)
	{
		bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)]);
	}
	for (int i = 0; i < 29; i++)
	{
		bwRemove[i] = rotrone(charHashes[(int)complement(i)]);
	}
}
