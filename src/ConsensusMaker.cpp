#include "ConsensusMaker.h"
#include "ErrorMaskHelper.h"

ConsensusMaker::~ConsensusMaker()
{
	for (size_t i = 0; i < simpleSequenceMutexes.size(); i++)
	{
		delete simpleSequenceMutexes[i];
	}
}

void ConsensusMaker::init(const std::vector<size_t>& unitigLengths)
{
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		simpleSequenceMutexes.emplace_back(new std::mutex);
	}
	compressedSequences.resize(unitigLengths.size());
	simpleCounts.resize(unitigLengths.size());
	parent.resize(unitigLengths.size());
	longestLeftOverlap.resize(unitigLengths.size(), 0);
	longestRightOverlap.resize(unitigLengths.size(), 0);
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		assert(unitigLengths[i] >= 1);
		compressedSequences[i].resize(unitigLengths[i]);
		simpleCounts[i].resize(unitigLengths[i]);
	}
	stringIndex.init(maxCode());
}

std::tuple<size_t, size_t, bool> ConsensusMaker::find(size_t unitig, size_t index)
{
	if (index >= longestLeftOverlap[unitig] && index < unitigLength(unitig) - longestRightOverlap[unitig])
	{
		return std::make_tuple(unitig, index, true);
	}
	std::tuple<size_t, size_t, bool> result;
	auto found = getParent(unitig, index);
	if (std::get<0>(found) == unitig && std::get<1>(found) == index)
	{
		assert(std::get<2>(found) == true);
		return std::make_tuple(unitig, index, true);
	}
	auto finalParent = find(std::get<0>(found), std::get<1>(found));
	assert(std::get<1>(finalParent) < longestLeftOverlap[std::get<0>(finalParent)] || std::get<1>(finalParent) >= unitigLength(std::get<0>(finalParent)) - longestRightOverlap[std::get<0>(finalParent)]);
	result = std::make_tuple(std::get<0>(finalParent), std::get<1>(finalParent), std::get<2>(finalParent) == std::get<2>(found));
	if (longestLeftOverlap[unitig] + longestRightOverlap[unitig] < unitigLength(unitig))
	{
		assert(index >= longestLeftOverlap[unitig] || index < unitigLength(unitig) - longestRightOverlap[unitig]);
		if (index < longestLeftOverlap[unitig])
		{
			assert(index < parent[unitig].size());
			parent[unitig][index] = result;
		}
		else
		{
			assert(index >= unitigLength(unitig) - longestRightOverlap[unitig]);
			size_t parentIndex = longestLeftOverlap[unitig] + index - (unitigLength(unitig) - longestRightOverlap[unitig]);
			assert(parentIndex < parent[unitig].size());
			parent[unitig][parentIndex] = result;
		}
	}
	else
	{
		assert(index < parent[unitig].size());
		parent[unitig][index] = result;
	}
	return result;
}

void ConsensusMaker::merge(size_t leftUnitig, size_t leftIndex, size_t rightUnitig, size_t rightIndex, bool fw)
{
	assert(leftIndex < longestLeftOverlap[leftUnitig] || leftIndex >= unitigLength(leftUnitig) - longestRightOverlap[leftUnitig]);
	assert(rightIndex < longestLeftOverlap[rightUnitig] || rightIndex >= unitigLength(rightUnitig) - longestRightOverlap[rightUnitig]);
	auto leftParent = find(leftUnitig, leftIndex);
	auto rightParent = find(rightUnitig, rightIndex);
	if (std::get<0>(rightParent) < std::get<0>(leftParent) || (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) > std::get<1>(rightParent)))
	{
		std::swap(leftParent, rightParent);
	}
	if (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) == std::get<1>(rightParent))
	{
		if (fw != (std::get<2>(leftParent) == std::get<2>(rightParent)))
		{
			needsComplementVerification.emplace_back(std::get<0>(leftParent), std::get<1>(leftParent));
		}
		return;
	}
	assert(std::get<0>(leftParent) < std::get<0>(rightParent) || (std::get<0>(leftParent) == std::get<0>(rightParent) && std::get<1>(leftParent) < std::get<1>(rightParent)));
	size_t leftParentUnitig = std::get<0>(leftParent);
	size_t leftParentIndex = std::get<1>(leftParent);
	bool leftParentOrient = std::get<2>(leftParent);
	size_t rightParentUnitig = std::get<0>(rightParent);
	size_t rightParentIndex = std::get<1>(rightParent);
	bool rightParentOrient = std::get<2>(rightParent);
	assert(leftParentIndex < longestLeftOverlap[leftParentUnitig] || leftParentIndex >= unitigLength(leftParentUnitig) - longestRightOverlap[leftParentUnitig]);
	assert(rightParentIndex < longestLeftOverlap[rightParentUnitig] || rightParentIndex >= unitigLength(rightParentUnitig) - longestRightOverlap[rightParentUnitig]);
	if (longestLeftOverlap[rightParentUnitig] + longestRightOverlap[rightParentUnitig] < unitigLength(rightParentUnitig))
	{
		if (rightParentIndex >= longestLeftOverlap[rightParentUnitig])
		{
			assert(rightParentIndex >= unitigLength(rightParentUnitig) - longestRightOverlap[rightParentUnitig]);
			rightParentIndex = longestLeftOverlap[rightParentUnitig] + rightParentIndex - (unitigLength(rightParentUnitig) - longestRightOverlap[rightParentUnitig]);
		}
	}
	assert(rightParentIndex < parent[rightParentUnitig].size());
	parent[rightParentUnitig][rightParentIndex] = std::make_tuple(leftParentUnitig, leftParentIndex, fw == (leftParentOrient == rightParentOrient));
}

void ConsensusMaker::addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap)
{
	for (size_t i = 0; i < overlap; i++)
	{
		size_t fromIndex = unitigLength(from.first) - overlap + i;
		if (!from.second) fromIndex = overlap - i - 1;
		size_t toIndex = i;
		if (!to.second) toIndex = unitigLength(to.first) - 1 - i;
		merge(from.first, fromIndex, to.first, toIndex, from.second == to.second);
	}
}

void ConsensusMaker::prepareEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap)
{
	if (from.second)
	{
		longestRightOverlap[from.first] = std::max(longestRightOverlap[from.first], overlap);
	}
	else
	{
		longestLeftOverlap[from.first] = std::max(longestLeftOverlap[from.first], overlap);
	}
	if (to.second)
	{
		longestLeftOverlap[to.first] = std::max(longestLeftOverlap[to.first], overlap);
	}
	else
	{
		longestRightOverlap[to.first] = std::max(longestRightOverlap[to.first], overlap);
	}
}

void ConsensusMaker::allocateParent()
{
	for (size_t i = 0; i < parent.size(); i++)
	{
		parent[i].resize(std::min(unitigLength(i), longestLeftOverlap[i] + longestRightOverlap[i]));
		if (longestLeftOverlap[i] + longestRightOverlap[i] < unitigLength(i))
		{
			for (size_t j = 0; j < longestLeftOverlap[i]; j++)
			{
				parent[i][j] = std::make_tuple(i, j, true);
			}
			for (size_t j = longestLeftOverlap[i]; j < parent[i].size(); j++)
			{
				parent[i][j] = std::make_tuple(i, (j - longestLeftOverlap[i]) + (unitigLength(i) - longestRightOverlap[i]), true);
			}
		}
		else
		{
			for (size_t j = 0; j < parent[i].size(); j++)
			{
				parent[i][j] = std::make_tuple(i, j, true);
			}
		}
	}
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> ConsensusMaker::getSequences()
{
	for (auto pair : needsComplementVerification)
	{
		auto found = find(pair.first, pair.second);
		assert(complement(compressedSequences[std::get<0>(found)].get(std::get<1>(found))) == compressedSequences[std::get<0>(found)].get(std::get<1>(found)));
	}
	for (size_t i = 0; i < simpleCounts.size(); i++)
	{
		for (size_t j = 0; j < simpleCounts[i].size(); j++)
		{
			auto found = find(i, j);
			size_t realI = std::get<0>(found);
			size_t realJ = std::get<1>(found);
			if (realI != i || realJ != j)
			{
				compressedSequences[i].set(j, compressedSequences[realI].get(realJ));
				if (!std::get<2>(found))
				{
					compressedSequences[i].set(j, complement(compressedSequences[i].get(j)));
				}
			}
			else
			{
				assert(std::get<2>(found));
			}
		}
	}
	stringIndex.buildReverseIndex();
	std::vector<CompressedSequenceType> result;
	result.resize(simpleCounts.size());
	for (size_t i = 0; i < simpleCounts.size(); i++)
	{
		std::vector<uint8_t> simpleExpanded;
		simpleExpanded.resize(simpleCounts[i].size(), 0);
		result[i].setCompressedAndClearInputVectorAndResizeExpanded(compressedSequences[i]);
		for (size_t j = 0; j < simpleCounts[i].size(); j++)
		{
			auto found = find(i, j);
			size_t realI = std::get<0>(found);
			size_t realJ = std::get<1>(found);
			size_t maxCount = simpleCounts[realI][realJ].second;
			uint32_t maxIndex = simpleCounts[realI][realJ].first;
			uint16_t compressed = result[i].getCompressed(j);
			if (complexCounts.count(realI) == 1)
			{
				if (complexCounts.at(realI).count(realJ) == 1)
				{
					for (auto pair : complexCounts.at(realI).at(realJ))
					{
						uint32_t index = pair.first;
						uint32_t count = pair.second;
						if (index == simpleCounts[realI][realJ].first) count += (uint32_t)(simpleCounts[realI][realJ].second);
						if (count < maxCount) continue;
						if (count == maxCount && stringIndex.getString(compressed, index) < stringIndex.getString(compressed, maxIndex)) continue;
						maxIndex = index;
						maxCount = count;
					}
				}
			}
			if (!std::get<2>(found))
			{
				maxIndex = stringIndex.getReverseIndex(compressed, maxIndex);
			}
			assert(maxCount > 0);
			assert(stringIndex.getString(compressed, maxIndex) != "");
			result[i].setCompressed(j, compressed);
			result[i].setExpanded(j, maxIndex);
		}
	}
	assert(result.size() == compressedSequences.size());
	return std::make_pair(std::move(result), stringIndex);
}

size_t ConsensusMaker::unitigLength(size_t unitig) const
{
	return simpleCounts[unitig].size();
}

std::tuple<size_t, size_t, bool> ConsensusMaker::getParent(size_t unitig, size_t index) const
{
	if (index >= longestLeftOverlap[unitig] && index < unitigLength(unitig) - longestRightOverlap[unitig])
	{
		return std::make_tuple(unitig, index, true);
	}
	if (longestLeftOverlap[unitig] + longestRightOverlap[unitig] >= unitigLength(unitig))
	{
		assert(parent[unitig].size() == unitigLength(unitig));
		assert(index < parent[unitig].size());
		return parent[unitig][index];
	}
	assert(parent[unitig].size() < unitigLength(unitig));
	assert(index >= longestLeftOverlap[unitig] || index < unitigLength(unitig) - longestRightOverlap[unitig]);
	if (index < longestLeftOverlap[unitig])
	{
		assert(index < parent[unitig].size());
		assert(std::get<1>(parent[unitig][index]) < longestLeftOverlap[std::get<0>(parent[unitig][index])] || std::get<1>(parent[unitig][index]) >= unitigLength(std::get<0>(parent[unitig][index])) - longestRightOverlap[std::get<0>(parent[unitig][index])]);
		return parent[unitig][index];
	}
	assert(index >= unitigLength(unitig) - longestRightOverlap[unitig]);
	index = longestLeftOverlap[unitig] + index - (unitigLength(unitig) - longestRightOverlap[unitig]);
	assert(index < parent[unitig].size());
	assert(std::get<1>(parent[unitig][index]) < longestLeftOverlap[std::get<0>(parent[unitig][index])] || std::get<1>(parent[unitig][index]) >= unitigLength(std::get<0>(parent[unitig][index])) - longestRightOverlap[std::get<0>(parent[unitig][index])]);
	return parent[unitig][index];
}

void ConsensusMaker::findParentLinks()
{
	for (size_t i = 0; i < parent.size(); i++)
	{
		if (longestLeftOverlap[i] + longestRightOverlap[i] < unitigLength(i))
		{
			assert(parent[i].size() < unitigLength(i));
			for (size_t j = 0; j < longestLeftOverlap[i]; j++)
			{
				find(i, j);
			}
			for (size_t j = 0; j < longestRightOverlap[i]; j++)
			{
				find(i, unitigLength(i) - longestRightOverlap[i] + j);
			}
		}
		else
		{
			assert(parent[i].size() == unitigLength(i));
			for (size_t j = 0; j < parent[i].size(); j++)
			{
				find(i, j);
			}
		}
	}
}

void ConsensusMaker::prepareHpcVariants()
{
	for (auto pair : needsComplementVerification)
	{
		auto found = find(pair.first, pair.second);
		assert(complement(compressedSequences[std::get<0>(found)].get(std::get<1>(found))) == compressedSequences[std::get<0>(found)].get(std::get<1>(found)));
	}
	for (size_t i = 0; i < simpleCounts.size(); i++)
	{
		for (size_t j = 0; j < simpleCounts[i].size(); j++)
		{
			auto found = find(i, j);
			size_t realI = std::get<0>(found);
			size_t realJ = std::get<1>(found);
			if (realI != i || realJ != j)
			{
				compressedSequences[i].set(j, compressedSequences[realI].get(realJ));
				if (!std::get<2>(found))
				{
					compressedSequences[i].set(j, complement(compressedSequences[i].get(j)));
				}
			}
			else
			{
				assert(std::get<2>(found));
			}
		}
	}
	stringIndex.buildReverseIndex();
}

std::vector<std::pair<size_t, std::vector<size_t>>> ConsensusMaker::getHpcVariants(const size_t unitig, const size_t minCoverage)
{
	std::vector<std::pair<size_t, std::vector<size_t>>> result;
	for (size_t j = 0; j < simpleCounts[unitig].size(); j++)
	{
		auto found = find(unitig, j);
		size_t realI = std::get<0>(found);
		size_t realJ = std::get<1>(found);
		uint16_t compressed = compressedSequences[unitig].get(j);
		std::unordered_map<size_t, size_t> lengthCounts;
		if (simpleCounts[realI][realJ].second > 0)
		{
			lengthCounts[stringIndex.getString(compressed, simpleCounts[realI][realJ].first).size()] = simpleCounts[realI][realJ].second;
		}
		if (complexCounts.count(realI) == 1)
		{
			if (complexCounts.at(realI).count(realJ) == 1)
			{
				for (auto pair : complexCounts.at(realI).at(realJ))
				{
					uint32_t index = pair.first;
					uint32_t count = pair.second;
					lengthCounts[stringIndex.getString(compressed, index).size()] += count;
				}
			}
		}
		std::vector<size_t> hpcVariants;
		for (auto pair : lengthCounts)
		{
			if (pair.second >= minCoverage) hpcVariants.push_back(pair.first);
		}
		std::sort(hpcVariants.begin(), hpcVariants.end());
		if (hpcVariants.size() >= 2)
		{
			result.emplace_back(j, std::move(hpcVariants));
		}
	}
	return result;
}

uint16_t ConsensusMaker::getCompressed(const size_t unitig, const size_t offset) const
{
	return compressedSequences[unitig].get(offset);
}
