#ifndef LittleBigVector_h
#define LittleBigVector_h

#include <cassert>
#include <limits>
#include <vector>
#include <phmap.h>

template <typename LittleType, typename BigType>
class LittleBigVector
{
	// static_assert(sizeof(BigType) > sizeof(LittleType), "Big type must be bigger than little type");
public:
	void resize(size_t size)
	{
		littles.resize(size);
	}
	void resize(size_t size, BigType v)
	{
		assert(v <= (BigType)std::numeric_limits<LittleType>::max());
		littles.resize(size, (LittleType)v);
	}
	BigType get(size_t i) const
	{
		auto found = bigs.find(i);
		if (found != bigs.end()) return found->second;
		return littles[i];
	}
	void set(size_t i, BigType v)
	{
		if (v <= (BigType)std::numeric_limits<LittleType>::max())
		{
			auto found = bigs.find(i);
			if (found != bigs.end()) bigs.erase(found);
			littles[i] = (LittleType)v;
			return;
		}
		bigs[i] = v;
	}
	void emplace_back(BigType v)
	{
		if (v <= (BigType)std::numeric_limits<LittleType>::max())
		{
			littles.emplace_back((LittleType)v);
			return;
		}
		littles.emplace_back(0);
		set(littles.size()-1, v);
	}
	size_t size() const
	{
		return littles.size();
	}
private:
	std::vector<LittleType> littles;
	phmap::flat_hash_map<uint32_t, BigType> bigs;
};

#endif
