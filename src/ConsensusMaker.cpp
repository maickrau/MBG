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
	stringIndices.resize(maxCode());
}

std::vector<CompressedSequenceType> ConsensusMaker::getSequences()
{
	std::vector<CompressedSequenceType> result;
	result.reserve(compressedSequences.size());
	for (size_t i = 0; i < compressedSequences.size(); i++)
	{
		std::vector<std::string> expanded;
		std::vector<std::vector<std::pair<uint32_t, uint32_t>>> complexCountsHere;
		complexCountsHere.resize(compressedSequences[i].size());
		for (auto pair : complexCounts[i])
		{
			complexCountsHere[pair.first.first].emplace_back(pair.first.second, pair.second);
		}
		for (size_t j = 0; j < compressedSequences[i].size(); j++)
		{
			size_t maxCount = simpleCounts[i][j].second;
			uint32_t maxIndex = simpleCounts[i][j].first;
			for (auto pair : complexCountsHere[j])
			{
				if (pair.first == (uint32_t)simpleCounts[i][j].first) pair.second += (uint32_t)simpleCounts[i][j].second;
				if (pair.second <= maxCount) continue;
				maxIndex = pair.first;
				maxCount = pair.second;
			}
			assert(maxCount > 0);
			expanded.push_back(getString(compressedSequences[i][j], maxIndex));
		}
		assert(compressedSequences[i].size() >= 1);
		assert(compressedSequences[i].size() == expanded.size());
		result.emplace_back(compressedSequences[i], expanded);
		{
			std::vector<uint16_t> tmp;
			std::swap(compressedSequences[i], tmp);
			phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint32_t> tmp2;
			std::swap(complexCounts[i], tmp2);
			std::vector<std::pair<uint8_t, uint8_t>> tmp3;
			std::swap(simpleCounts[i], tmp3);
		}
	}
	assert(result.size() == compressedSequences.size());
	return result;
}

std::string ConsensusMaker::getString(uint16_t compressed, uint32_t index) const
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
	// todo optimize
	for (auto pair : stringIndices[compressed])
	{
		if (pair.second == index) return pair.first;
	}
	assert(false);
	return "";
}

uint32_t ConsensusMaker::getStringIndex(uint16_t compressed, const std::string& expanded)
{
	if (compressed >= 0 && compressed <= 3)
	{
		return expanded.size();
	}
	auto found = stringIndices[compressed].find(expanded);
	if (found != stringIndices[compressed].end())
	{
		return found->second;
	}
	uint32_t result = stringIndices[compressed].size();
	stringIndices[compressed][expanded] = result;
	return result;
}
