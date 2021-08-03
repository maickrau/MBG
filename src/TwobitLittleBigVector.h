#ifndef TwobitLittleBigVector_h
#define TwobitLittleBigVector_h

#include <cassert>
#include <limits>
#include <vector>
#include <phmap.h>

template <typename BigType>
class TwobitLittleBigVector
{
public:
	void resize(size_t size)
	{
		realSize = size;
		littles.resize((realSize+3)/4);
	}
	void resize(size_t size, BigType v)
	{
		assert(v >= 0 && v <= 3);
		realSize = size;
		uint8_t pattern = (v << 6) + (v << 4) + (v << 2) + v;
		littles.resize((realSize+3)/4, pattern);
	}
	BigType get(size_t i) const
	{
		auto found = bigs.find(i);
		if (found != bigs.end()) return found->second;
		size_t index = i / 4;
		size_t offset = (i % 4) * 2;
		return (littles[index] >> offset) & 3;
	}
	void set(size_t i, BigType v)
	{
		if (v >= 0 && v <= 3)
		{
			auto found = bigs.find(i);
			if (found != bigs.end()) bigs.erase(found);
			size_t index = i / 4;
			size_t offset = (i % 4) * 2;
			uint8_t removeMask = ~(3 << offset);
			littles[index] &= removeMask;
			uint8_t addMask = (uint8_t)v << offset;
			littles[index] |= addMask;
			return;
		}
		bigs[i] = v;
	}
	void emplace_back(BigType v)
	{
		realSize += 1;
		if (realSize % 4 == 0) littles.emplace_back(0);
		set(realSize-1, v);
	}
	size_t size() const
	{
		return realSize;
	}
private:
	size_t realSize;
	std::vector<uint8_t> littles;
	phmap::flat_hash_map<uint32_t, BigType> bigs;
};

#endif
