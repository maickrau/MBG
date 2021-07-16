#ifndef MostlySparse2DHashmap_h
#define MostlySparse2DHashmap_h

#include <algorithm>
#include <vector>
#include <limits>
#include <tuple>
#include <phmap.h>
#include "VectorWithDirection.h"

template <typename ValueType>
class MostlySparse2DHashmap
{
public:
	void resize(size_t newSize)
	{
		size_t oldSize = size();
		firstKeyValue.resize(newSize);
		for (size_t i = oldSize; i < newSize; i++)
		{
			firstKeyValue[std::make_pair(i, true)].first.first = std::numeric_limits<size_t>::max();
			firstKeyValue[std::make_pair(i, false)].first.first = std::numeric_limits<size_t>::max();
		}
	}
	bool hasValue(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (firstKeyValue[from].first == to) return true;
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
		firstKeyValue[std::make_pair(firstKeyValue.size()-1, true)].first.first = std::numeric_limits<size_t>::max();
		firstKeyValue[std::make_pair(firstKeyValue.size()-1, false)].first.first = std::numeric_limits<size_t>::max();
	}
	ValueType get(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		if (firstKeyValue[from].first == to) return firstKeyValue[from].second;
		auto found = additionalKeyValues.find(from);
		assert(found != additionalKeyValues.end());
		auto found2 = found->second.find(to);
		assert(found2 != found->second.end());
		return found2->second;
	}
	void set(std::pair<size_t, bool> from, std::pair<size_t, bool> to, ValueType value)
	{
		if (firstKeyValue[from].first.first == std::numeric_limits<size_t>::max())
		{
			firstKeyValue[from].first = to;
			firstKeyValue[from].second = value;
			return;
		}
		if (firstKeyValue[from].first == to)
		{
			firstKeyValue[from].second = value;
			return;
		}
		additionalKeyValues[from][to] = value;
	}
	std::vector<std::pair<std::pair<size_t, bool>, ValueType>> getValues(std::pair<size_t, bool> from) const
	{
		std::vector<std::pair<std::pair<size_t, bool>, ValueType>> result;
		if (firstKeyValue[from].first.first != std::numeric_limits<size_t>::max()) result.emplace_back(firstKeyValue[from]);
		if (additionalKeyValues.count(from) == 1)
		{
			result.insert(result.end(), additionalKeyValues.at(from).begin(), additionalKeyValues.at(from).end());
		}
		std::sort(result.begin(), result.end(), [](std::pair<std::pair<size_t, bool>, ValueType> left, std::pair<std::pair<size_t, bool>, ValueType> right) { return left.first.first < right.first.first || (left.first.first == right.first.first && left.first.second < right.first.second); });
		return result;
	}
private:
	VectorWithDirection<std::pair<std::pair<size_t, bool>, ValueType>> firstKeyValue;
	phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_map<std::pair<size_t, bool>, ValueType>> additionalKeyValues;
};

#endif
