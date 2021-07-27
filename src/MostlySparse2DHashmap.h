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
		firstKeyValue.resize(newSize);
		for (size_t i = oldSize; i < newSize; i++)
		{
			firstKeyValue[std::make_pair(i, true)].first = std::numeric_limits<uint32_t>::max();
			firstKeyValue[std::make_pair(i, false)].first = std::numeric_limits<uint32_t>::max();
		}
	}
	bool hasValue(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (pairToInt(to) < std::numeric_limits<uint32_t>::max())
		{
			if (firstKeyValue[from].first == pairToInt(to)) return true;
		}
		auto found = additionalKeyValues.find(from);
		if (found == additionalKeyValues.end()) return false;
		auto found2 = found->second.find(to);
		return found2 != found->second.end();
	}
	size_t size() const
	{
		return firstKeyValue.size();
	}
	void emplace_back()
	{
		firstKeyValue.emplace_back();
		firstKeyValue[std::make_pair(firstKeyValue.size()-1, true)].first = std::numeric_limits<uint32_t>::max();
		firstKeyValue[std::make_pair(firstKeyValue.size()-1, false)].first = std::numeric_limits<uint32_t>::max();
	}
	BigType get(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (pairToInt(to) < std::numeric_limits<uint32_t>::max() && firstKeyValue[from].first == pairToInt(to)) return (BigType)firstKeyValue[from].second;
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
			if (firstKeyValue[from].first == std::numeric_limits<uint32_t>::max())
			{
				firstKeyValue[from].first = pairToInt(to);
				firstKeyValue[from].second = (SmallType)value;
				setted = true;
			}
			if (firstKeyValue[from].first == pairToInt(to))
			{
				firstKeyValue[from].second = (SmallType)value;
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
			if (pairToInt(to) < std::numeric_limits<uint32_t>::max() && firstKeyValue[from].first == pairToInt(to)) firstKeyValue[from].first = std::numeric_limits<uint32_t>::max();
		}
		additionalKeyValues[from][to] = value;
	}
	std::vector<std::pair<std::pair<size_t, bool>, BigType>> getValues(std::pair<size_t, bool> from) const
	{
		std::vector<std::pair<std::pair<size_t, bool>, BigType>> result;
		if (firstKeyValue[from].first != std::numeric_limits<uint32_t>::max()) result.emplace_back(intToPair(firstKeyValue[from].first), firstKeyValue[from].second);
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
	VectorWithDirection<std::pair<uint32_t, SmallType>> firstKeyValue;
	phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_map<std::pair<size_t, bool>, BigType>> additionalKeyValues;
};

#endif
