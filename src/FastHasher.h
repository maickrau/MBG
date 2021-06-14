#ifndef FastHasher_h
#define FastHasher_h

#include <cstdint>
#include <cstring>
#include <algorithm>

constexpr uint64_t combineHash(uint64_t first, uint64_t second)
{
	return ((first << 32) | (first >> 32)) ^ second;
}

class FastHasher
{
public:
	FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash);
	FastHasher(size_t kmerSize);
	__attribute__((always_inline))
	inline void addChar(char c)
	{
		fwHash = rotlone(fwHash) ^ fwAddHash(c);
		bwHash = rotrone(bwHash) ^ bwAddHash(c);
	}
	__attribute__((always_inline))
	inline void removeChar(char c)
	{
		fwHash ^= fwRemoveHash(c);
		bwHash ^= bwRemoveHash(c);
	}
	__attribute__((always_inline))
	inline uint64_t hash() const
	{
		return std::min(fwHash, bwHash);
	}
	__attribute__((always_inline))
	inline uint64_t getFwHash() const
	{
		return fwHash;
	}
	__attribute__((always_inline))
	inline uint64_t getBwHash() const
	{
		return bwHash;
	}
private:
	__attribute__((always_inline))
	inline uint64_t rotlone(uint64_t val) const
	{
		return (val << 1) | (val >> (64-1));
	};
	__attribute__((always_inline))
	inline uint64_t rotrone(uint64_t val) const
	{
		return (val >> 1) | (val << (64-1));
	};
	__attribute__((always_inline))
	inline uint64_t rotlk(uint64_t val) const
	{
		return (val << kmerSize) | (val >> (64-kmerSize));
	};
	__attribute__((always_inline))
	inline uint64_t rotlkmin1(uint64_t val) const
	{
		return (val << (kmerSize-1)) | (val >> (64-(kmerSize-1)));
	};
	void precalcRots();
	__attribute__((always_inline))
	inline uint64_t fwAddHash(char c)
	{
		return fwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t fwRemoveHash(char c)
	{
		return fwRemove[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwAddHash(char c)
	{
		return bwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwRemoveHash(char c)
	{
		return bwRemove[(int)c];
	}
	uint64_t fwAdd[29];
	uint64_t fwRemove[29];
	uint64_t bwAdd[29];
	uint64_t bwRemove[29];
	// https://bioinformatics.stackexchange.com/questions/19/are-there-any-rolling-hash-functions-that-can-hash-a-dna-sequence-and-its-revers
	static constexpr uint64_t hashA = 0x3c8bfbb395c60474;
	static constexpr uint64_t hashC = 0x3193c18562a02b4c;
	static constexpr uint64_t hashG = 0x20323ed082572324;
	static constexpr uint64_t hashT = 0x295549f54be24456;
	uint64_t charHashes[29] { 0, 0x3c8bfbb395c60474, 0x3193c18562a02b4c, 0x20323ed082572324, 0x295549f54be24456, combineHash(hashA, hashC), combineHash(hashA, hashG), combineHash(hashA, hashT), combineHash(hashC, hashA), combineHash(hashC, hashG), combineHash(hashC, hashT), combineHash(hashG, hashA), combineHash(hashG, hashC), combineHash(hashG, hashT), combineHash(hashT, hashA), combineHash(hashT, hashC), combineHash(hashT, hashG), combineHash(hashA, hashC), combineHash(hashA, hashG), combineHash(hashA, hashT), combineHash(hashC, hashA), combineHash(hashC, hashG), combineHash(hashC, hashT), combineHash(hashG, hashA), combineHash(hashG, hashC), combineHash(hashG, hashT), combineHash(hashT, hashA), combineHash(hashT, hashC), combineHash(hashT, hashG) };
	uint64_t fwHash;
	uint64_t bwHash;
	size_t kmerSize;
};

#endif
