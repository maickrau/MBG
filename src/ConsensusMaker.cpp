#include "ConsensusMaker.h"
#include "ErrorMaskHelper.h"

void ConsensusMaker::init(const std::vector<size_t>& unitigLengths)
{
	compressedSequences.resize(unitigLengths.size());
	simpleCounts.resize(unitigLengths.size());
	parent.resize(unitigLengths.size());
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		assert(unitigLengths[i] >= 1);
		compressedSequences[i].resize(unitigLengths[i]);
		parent[i].resize(unitigLengths[i]);
		for (size_t j = 0; j < unitigLengths[i]; j++)
		{
			parent[i][j] = std::make_tuple(i, j, true);
		}
		simpleCounts[i].resize(unitigLengths[i]);
	}
	stringIndex.init(maxCode());
}

std::tuple<size_t, size_t, bool> ConsensusMaker::find(size_t unitig, size_t index) const
{
	std::tuple<size_t, size_t, bool> result;
	auto found = parent[unitig][index];
	if (std::get<0>(found) == unitig && std::get<1>(found) == index)
	{
		assert(std::get<2>(found) == true);
		return std::make_tuple(unitig, index, true);
	}
	auto finalParent = find(std::get<0>(found), std::get<1>(found));
	parent[unitig][index] = std::make_tuple(std::get<0>(finalParent), std::get<1>(finalParent), std::get<2>(finalParent) == std::get<2>(found));
	return parent[unitig][index];
}

void ConsensusMaker::merge(size_t leftUnitig, size_t leftIndex, size_t rightUnitig, size_t rightIndex, bool fw)
{
	auto leftParent = find(leftUnitig, leftIndex);
	auto rightParent = find(rightUnitig, rightIndex);
	if (std::get<0>(rightParent) < std::get<0>(leftParent) || (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) > std::get<1>(rightParent)))
	{
		std::swap(leftParent, rightParent);
	}
	if (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) == std::get<1>(rightParent))
	{
		assert(fw == (std::get<2>(leftParent) == std::get<2>(rightParent)));
		return;
	}
	assert(std::get<0>(leftParent) < std::get<0>(rightParent) || (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) < std::get<1>(rightParent)));
	parent[std::get<0>(rightParent)][std::get<1>(rightParent)] = std::make_tuple(std::get<0>(leftParent), std::get<1>(leftParent), fw == (std::get<2>(leftParent) == std::get<2>(rightParent)));
}

void ConsensusMaker::addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap)
{
	for (size_t i = 0; i < overlap; i++)
	{
		size_t fromIndex = parent[from.first].size() - overlap + i;
		if (!from.second) fromIndex = overlap - i - 1;
		size_t toIndex = i;
		if (!to.second) toIndex = parent[to.first].size() - 1 - i;
		merge(from.first, fromIndex, to.first, toIndex, from.second == to.second);
	}
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> ConsensusMaker::getSequences()
{
	stringIndex.buildReverseIndex();
	std::vector<CompressedSequenceType> result;
	result.reserve(simpleCounts.size());
	for (size_t i = 0; i < simpleCounts.size(); i++)
	{
		std::vector<uint8_t> simpleExpanded;
		simpleExpanded.resize(simpleCounts[i].size(), 0);
		CompressedSequenceType compressedSequence;
		compressedSequence.resize(simpleExpanded.size());
		for (size_t j = 0; j < simpleCounts[i].size(); j++)
		{
			auto found = find(i, j);
			size_t realI = std::get<0>(found);
			size_t realJ = std::get<1>(found);
			size_t maxCount = simpleCounts[realI][realJ].second;
			uint32_t maxIndex = simpleCounts[realI][realJ].first;
			if (complexCounts.count(realI) == 1)
			{
				if (complexCounts.at(realI).count(realJ) == 1)
				{
					for (auto pair : complexCounts.at(realI).at(realJ))
					{
						uint32_t index = pair.first;
						uint32_t count = pair.second;
						if (index == simpleCounts[realI][realJ].first) count += (uint32_t)(simpleCounts[realI][realJ].second);
						if (count <= maxCount) continue;
						maxIndex = index;
						maxCount = count;
					}
				}
			}
			uint16_t compressed = compressedSequences[realI].get(realJ);
			if (!std::get<2>(found))
			{
				maxIndex = stringIndex.getReverseIndex(compressed, maxIndex);
				compressed = complement(compressed);
			}
			assert(maxCount > 0);
			// // if (maxCount == 0)
			// // {
			// 	compressedSequences[i][j / MutexLength].set(j % MutexLength, 0);
			// 	maxIndex = stringIndex.getIndex(0, "N");
			// // }
			assert(stringIndex.getString(compressedSequences[realI].get(realJ), maxIndex) != "");
			compressedSequence.setCompressed(j, compressed);
			compressedSequence.setExpanded(j, maxIndex);
		}
		result.emplace_back(std::move(compressedSequence));
	}
	assert(result.size() == compressedSequences.size());
	return std::make_pair(std::move(result), stringIndex);
}
