#include "ConsensusMaker.h"
#include "ErrorMaskHelper.h"

ConsensusMaker::~ConsensusMaker()
{
	for (size_t i = 0; i < seqMutexes.size(); i++)
	{
		for (size_t j = 0; j < seqMutexes[i].size(); j++)
		{
			delete seqMutexes[i][j];
		}
	}
}

void ConsensusMaker::init(const std::vector<size_t>& unitigLengths)
{
	compressedSequences.resize(unitigLengths.size());
	simpleCounts.resize(unitigLengths.size());
	complexCounts.resize(unitigLengths.size());
	seqMutexes.resize(unitigLengths.size());
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		assert(unitigLengths[i] >= 1);
		compressedSequences[i].resize(unitigLengths[i], 0);
		simpleCounts[i].resize(unitigLengths[i]);
		complexCountMutexes.emplace_back(new std::mutex);
		for (size_t j = 0; j < unitigLengths[i]; j += MutexLength)
		{
			seqMutexes[i].emplace_back(new std::mutex);
		}
	}
	stringIndex.init(maxCode());
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> ConsensusMaker::getSequences()
{
	stringIndex.buildReverseIndex();
	std::vector<CompressedSequenceType> result;
	result.reserve(compressedSequences.size());
	for (size_t i = 0; i < compressedSequences.size(); i++)
	{
		std::vector<uint8_t> simpleExpanded;
		simpleExpanded.resize(compressedSequences[i].size(), 0);
		std::vector<std::pair<uint32_t, uint32_t>> complexExpanded;
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> complexCountsHere;
		for (auto pair : complexCounts[i])
		{
			complexCountsHere.emplace_back(pair.first.first, pair.first.second, pair.second);
		}
		std::sort(complexCountsHere.begin(), complexCountsHere.end(), [](auto& left, auto& right) { return std::get<0>(left) > std::get<0>(right); });
		for (size_t j = 0; j < compressedSequences[i].size(); j++)
		{
			size_t maxCount = simpleCounts[i][j].second;
			uint32_t maxIndex = simpleCounts[i][j].first;
			assert(complexCountsHere.size() == 0 || std::get<0>(complexCountsHere.back()) >= j);
			while (complexCountsHere.size() > 0 && std::get<0>(complexCountsHere.back()) == j)
			{
				uint32_t index = std::get<1>(complexCountsHere.back());
				uint32_t count = std::get<2>(complexCountsHere.back());
				complexCountsHere.pop_back();
				if (index == (uint32_t)simpleCounts[i][j].first) count += (uint32_t)simpleCounts[i][j].second;
				if (count <= maxCount) continue;
				maxIndex = index;
				maxCount = count;
			}
			assert(maxCount > 0);
			assert(stringIndex.getString(compressedSequences[i].get(j), maxIndex) != "");
			if (maxIndex <= (uint32_t)std::numeric_limits<uint8_t>::max())
			{
				simpleExpanded[j] = maxIndex;
			}
			else
			{
				complexExpanded.emplace_back(j, maxIndex);
			}
		}
		assert(compressedSequences[i].size() >= 1);
		assert(compressedSequences[i].size() == simpleExpanded.size());
		result.emplace_back(std::move(compressedSequences[i]), std::move(simpleExpanded), std::move(complexExpanded));
	}
	assert(result.size() == compressedSequences.size());
	return std::make_pair(std::move(result), stringIndex);
}
