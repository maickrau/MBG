#include <vector>
#include "FastHasher.h"
#include "MBGCommon.h"

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
