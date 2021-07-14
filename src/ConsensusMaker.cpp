#include "ConsensusMaker.h"

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
	expandedCounts.resize(unitigLengths.size());
	seqMutexes.resize(unitigLengths.size());
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		assert(unitigLengths[i] >= 1);
		compressedSequences[i].resize(unitigLengths[i], 0);
		expandedCounts[i].resize(unitigLengths[i]);
		for (size_t j = 0; j < unitigLengths[i]; j += MutexLength)
		{
			seqMutexes[i].emplace_back(new std::mutex);
		}
	}
}

std::vector<CompressedSequenceType> ConsensusMaker::getSequences()
{
	std::vector<CompressedSequenceType> result;
	result.reserve(compressedSequences.size());
	for (size_t i = 0; i < compressedSequences.size(); i++)
	{
		std::vector<std::string> expanded;
		for (size_t j = 0; j < compressedSequences[i].size(); j++)
		{
			size_t maxCount = 0;
			std::string maxSeq = "";
			for (auto pair : expandedCounts[i][j])
			{
				if (pair.second <= maxCount) continue;
				maxSeq = pair.first;
				maxCount = pair.second;
			}
			assert(maxCount > 0);
			assert(maxSeq != "");
			expanded.push_back(maxSeq);
		}
		assert(compressedSequences[i].size() >= 1);
		assert(compressedSequences[i].size() == expanded.size());
		result.emplace_back(compressedSequences[i], expanded);
		{
			std::vector<uint16_t> tmp;
			std::swap(compressedSequences[i], tmp);
			std::vector<phmap::flat_hash_map<std::string, size_t>> tmp2;
			std::swap(expandedCounts[i], tmp2);
		}
	}
	assert(result.size() == compressedSequences.size());
	return result;
}
