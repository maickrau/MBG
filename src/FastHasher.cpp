#include "FastHasher.h"

#include <vector>

extern std::vector<uint16_t> multiRLEReverseComplements;

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

uint16_t FastHasher::complement(uint16_t c) const
{
	return multiRLEReverseComplements[c];
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
	charHashes.resize(multiRLEReverseComplements.size());
	fwAdd.resize(multiRLEReverseComplements.size());
	fwRemove.resize(multiRLEReverseComplements.size());
	bwAdd.resize(multiRLEReverseComplements.size());
	bwRemove.resize(multiRLEReverseComplements.size());
	for (size_t i = 0; i < multiRLEReverseComplements.size(); i++)
	{
		charHashes[i] = getHash(i);
	}
	for (size_t i = 0; i < multiRLEReverseComplements.size(); i++)
	{
		fwAdd[i] = charHashes[i];
	}
	for (size_t i = 0; i < multiRLEReverseComplements.size(); i++)
	{
		fwRemove[i] = rotlk(charHashes[i]);
	}
	for (size_t i = 0; i < multiRLEReverseComplements.size(); i++)
	{
		bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)]);
	}
	for (size_t i = 0; i < multiRLEReverseComplements.size(); i++)
	{
		bwRemove[i] = rotrone(charHashes[(int)complement(i)]);
	}
	precalcedK = kmerSize;
}
