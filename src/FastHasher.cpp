#include <vector>
#include <cassert>
#include "FastHasher.h"
#include "MBGCommon.h"
#include "ErrorMaskHelper.h"

std::mutex FastHasher::precalcMutex;
size_t FastHasher::precalcedK = 0;
std::vector<uint64_t> FastHasher::charHashes;
std::vector<uint64_t> FastHasher::fwAdd;
std::vector<uint64_t> FastHasher::fwRemove;
std::vector<uint64_t> FastHasher::bwAdd;
std::vector<uint64_t> FastHasher::bwRemove;

FastHasher::FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash) :
fwHash(fwHash),
bwHash(bwHash),
kmerSize(kmerSize % 64)
{
	std::lock_guard<std::mutex> guard { precalcMutex };
	if (charHashes.size() == 0 || kmerSize != precalcedK) precalcRots();
}

FastHasher::FastHasher(size_t kmerSize) :
fwHash(0),
bwHash(0),
kmerSize(kmerSize % 64)
{
	std::lock_guard<std::mutex> guard { precalcMutex };
	if (charHashes.size() == 0 || kmerSize != precalcedK) precalcRots();
}

// https://naml.us/post/inverse-of-a-hash-function/
uint64_t getHash(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

void FastHasher::precalcRots()
{
	assert((size_t)maxCode() < (size_t)std::numeric_limits<CharType>::max());
	charHashes.resize(maxCode());
	fwAdd.resize(maxCode());
	fwRemove.resize(maxCode());
	bwAdd.resize(maxCode());
	bwRemove.resize(maxCode());
	for (size_t i = 0; i < maxCode(); i++)
	{
		charHashes[i] = getHash(i);
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		fwAdd[i] = charHashes[i];
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		fwRemove[i] = rotlk(charHashes[i]);
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)]);
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		bwRemove[i] = rotrone(charHashes[(int)complement(i)]);
	}
	precalcedK = kmerSize;
}
