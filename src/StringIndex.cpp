#include "StringIndex.h"
#include "ErrorMaskHelper.h"

void StringIndex::init(size_t maxCode)
{
	index.resize(maxCode);
}

void StringIndex::buildReverseIndex()
{
	reverseIndex.resize(index.size());
	for (size_t i = 0; i < index.size(); i++)
	{
		reverseIndex[i].resize(index[i].size(), "");
		for (auto pair : index[i])
		{
			assert(reverseIndex[i][pair.second] == "");
			assert(pair.second < reverseIndex[i].size());
			reverseIndex[i][pair.second] = pair.first;
		}
		for (size_t j = 0; j < reverseIndex[i].size(); j++)
		{
			assert(reverseIndex[i][j] != "");
		}
	}
}

std::string StringIndex::getString(uint16_t compressed, uint32_t index) const
{
	if (compressed >= 0 && compressed <= 3)
	{
		std::string result;
		for (size_t i = 0; i < index; i++)
		{
			result += "ACGT"[compressed];
		}
		return result;
	}
	if (complement(compressed) < compressed)
	{
		compressed = complement(compressed);
		assert(compressed < reverseIndex.size());
		assert(index < reverseIndex[compressed].size());
		return revCompRaw(reverseIndex[compressed][index]);
	}
	else if (complement(compressed) == compressed)
	{
		assert(compressed < reverseIndex.size());
		assert(index < reverseIndex[compressed].size());
		std::string result = reverseIndex[compressed][index];
		if (result.back() == '_') result.pop_back();
		return result;
	}
	else
	{
		assert(compressed < reverseIndex.size());
		assert(index < reverseIndex[compressed].size());
		return reverseIndex[compressed][index];
	}
}

uint32_t StringIndex::getReverseIndex(uint16_t compressed, uint32_t index) const
{
	if (complement(compressed) != compressed) return index;
	return (index / 2) * 2 + (1 - (index % 2));
}

uint32_t StringIndex::getIndex(uint16_t compressed, std::string expanded)
{
	if (compressed >= 0 && compressed <= 3)
	{
		return expanded.size();
	}
	if (complement(compressed) < compressed)
	{
		compressed = complement(compressed);
		expanded = revCompRaw(expanded);
	}
	uint32_t maybeResult = index[compressed].size();
	auto found = index[compressed].find(expanded);
	if (found != index[compressed].end())
	{
		return found->second;
	}
	index[compressed][expanded] = maybeResult;
	// need to handle compressed sequences which are their own reverse complement
	if (complement(compressed) == compressed)
	{
		std::string revExpanded = revCompRaw(expanded);
		// ...and need to make sure that if the expanded sequence is its own reverse complement it doesn't break anything
		if (revExpanded == expanded) revExpanded = revExpanded + "_";
		assert(index[compressed].count(revExpanded) == 0);
		index[compressed][revExpanded] = maybeResult+1;
	}
	return maybeResult;
}
