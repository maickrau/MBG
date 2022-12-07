#ifndef MostlySparse2DHashmap_h
#define MostlySparse2DHashmap_h

#include <algorithm>
#include <vector>
#include <limits>
#include <tuple>
#include <phmap.h>
#include "VectorWithDirection.h"

template <typename SmallType, typename BigType>
class MostlySparse2DHashmap
{
public:
	void resize(size_t newSize)
	{
		size_t oldSize = size();
		firstKey.resize(newSize);
		firstValue.resize(newSize);
		for (size_t i = oldSize; i < newSize; i++)
		{
			firstKey[std::make_pair(i, true)] = std::numeric_limits<uint32_t>::max();
			firstKey[std::make_pair(i, false)] = std::numeric_limits<uint32_t>::max();
		}
	}
	bool hasValue(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (pairToInt(to) < std::numeric_limits<uint32_t>::max())
		{
			if (firstKey[from] == pairToInt(to)) return true;
		}
		auto found = additionalKeyValues.find(from);
		if (found == additionalKeyValues.end()) return false;
		auto found2 = found->second.find(to);
		return found2 != found->second.end();
	}
	size_t size() const
	{
		return firstKey.size();
	}
	void emplace_back()
	{
		firstKey.emplace_back();
		firstValue.emplace_back();
		firstKey[std::make_pair(firstKey.size()-1, true)] = std::numeric_limits<uint32_t>::max();
		firstKey[std::make_pair(firstKey.size()-1, false)] = std::numeric_limits<uint32_t>::max();
	}
	BigType get(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (pairToInt(to) < std::numeric_limits<uint32_t>::max() && firstKey[from] == pairToInt(to)) return (BigType)firstValue[from];
		auto found = additionalKeyValues.find(from);
		assert(found != additionalKeyValues.end());
		auto found2 = found->second.find(to);
		assert(found2 != found->second.end());
		return found2->second;
	}
	void set(std::pair<size_t, bool> from, std::pair<size_t, bool> to, BigType value)
	{
		if (pairToInt(to) < std::numeric_limits<uint32_t>::max() && value >= (BigType)std::numeric_limits<SmallType>::min() && value <= (BigType)std::numeric_limits<SmallType>::max())
		{
			bool setted = false;
			if (firstKey[from] == std::numeric_limits<uint32_t>::max())
			{
				firstKey[from] = pairToInt(to);
				firstValue[from] = (SmallType)value;
				setted = true;
			}
			if (firstKey[from] == pairToInt(to))
			{
				firstValue[from] = (SmallType)value;
				setted = true;
			}
			if (setted)
			{
				auto found = additionalKeyValues.find(from);
				if (found != additionalKeyValues.end())
				{
					auto found2 = found->second.find(to);
					if (found2 != found->second.end())
					{
						found->second.erase(found2);
					}
				}
				return;
			}
		}
		else
		{
			if (pairToInt(to) < std::numeric_limits<uint32_t>::max() && firstKey[from] == pairToInt(to)) firstKey[from] = std::numeric_limits<uint32_t>::max();
		}
		additionalKeyValues[from][to] = value;
	}
	std::vector<std::pair<std::pair<size_t, bool>, BigType>> getValues(std::pair<size_t, bool> from) const
	{
		std::vector<std::pair<std::pair<size_t, bool>, BigType>> result;
		if (firstKey[from] != std::numeric_limits<uint32_t>::max()) result.emplace_back(intToPair(firstKey[from]), firstValue[from]);
		if (additionalKeyValues.count(from) == 1)
		{
			result.insert(result.end(), additionalKeyValues.at(from).begin(), additionalKeyValues.at(from).end());
		}
		std::sort(result.begin(), result.end(), [](std::pair<std::pair<size_t, bool>, BigType> left, std::pair<std::pair<size_t, bool>, BigType> right) { return left.first.first < right.first.first || (left.first.first == right.first.first && left.first.second < right.first.second); });
		return result;
	}
private:
	uint32_t pairToInt(std::pair<size_t, bool> value) const
	{
		if (value.first >= (size_t)std::numeric_limits<uint32_t>::max() / 2) return std::numeric_limits<uint32_t>::max();
		return (uint32_t)value.first * 2 + (value.second ? 1 : 0);
	}
	std::pair<size_t, bool> intToPair(uint32_t value) const
	{
		return std::make_pair(value / 2, (value % 2) == 1);
	}
	VectorWithDirection<uint32_t> firstKey;
	VectorWithDirection<SmallType> firstValue;
	phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_map<std::pair<size_t, bool>, BigType>> additionalKeyValues;
};

#endif
