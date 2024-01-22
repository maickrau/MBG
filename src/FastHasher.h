#ifndef FastHasher_h
#define FastHasher_h

#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>
#include <mutex>
#include "MBGCommon.h"
#include "ErrorMaskHelper.h"

template <int Instance> // template so multiple hashers can be active with different k-mer sizes at the same time
class FastHasherInstance
{
public:
	FastHasherInstance(size_t kmerSize, uint64_t fwHash, uint64_t bwHash) :
		fwHash(fwHash),
		bwHash(bwHash),
		kmerSize(kmerSize % 64)
	{
		std::lock_guard<std::mutex> guard { precalcMutex };
		if (charHashes.size() == 0 || this->kmerSize != precalcedK) precalcRots();
	}
	FastHasherInstance(size_t kmerSize) :
		fwHash(0),
		bwHash(0),
		kmerSize(kmerSize % 64)
	{
		std::lock_guard<std::mutex> guard { precalcMutex };
		if (charHashes.size() == 0 || this->kmerSize != precalcedK) precalcRots();
	}
	__attribute__((always_inline))
	inline void addChar(CharType c)
	{
		fwHash = rotlone(fwHash) ^ fwAddHash(c);
		bwHash = rotrone(bwHash) ^ bwAddHash(c);
	}
	__attribute__((always_inline))
	inline void removeChar(CharType c)
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
	__attribute__((always_inline))
	inline uint64_t fwAddHash(CharType c)
	{
		return fwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t fwRemoveHash(CharType c)
	{
		return fwRemove[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwAddHash(CharType c)
	{
		return bwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwRemoveHash(CharType c)
	{
		return bwRemove[(int)c];
	}
	void precalcRots()
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
	uint64_t fwHash;
	uint64_t bwHash;
	size_t kmerSize;
	static std::mutex precalcMutex;
	static size_t precalcedK;
	static std::vector<uint64_t> charHashes;
	static std::vector<uint64_t> fwAdd;
	static std::vector<uint64_t> fwRemove;
	static std::vector<uint64_t> bwAdd;
	static std::vector<uint64_t> bwRemove;
};

template <int Instance>
size_t FastHasherInstance<Instance>::precalcedK = 0;
template <int Instance>
std::vector<uint64_t> FastHasherInstance<Instance>::charHashes;
template <int Instance>
std::vector<uint64_t> FastHasherInstance<Instance>::fwAdd;
template <int Instance>
std::vector<uint64_t> FastHasherInstance<Instance>::fwRemove;
template <int Instance>
std::vector<uint64_t> FastHasherInstance<Instance>::bwAdd;
template <int Instance>
std::vector<uint64_t> FastHasherInstance<Instance>::bwRemove;
template <int Instance>
std::mutex FastHasherInstance<Instance>::precalcMutex;

using FastHasher = FastHasherInstance<0>;

#endif
