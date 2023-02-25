#include <unordered_map>
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
	totalLength = 0;
	for (size_t i = 0; i < unitigLengths.size(); i++)
	{
		totalLength += unitigLengths[i];
		assert(unitigLengths[i] >= 1);
		compressedSequences[i].resize(unitigLengths[i]);
		simpleCounts[i].resize(unitigLengths[i]);
	}
	stringIndex.init(maxCode());
}

void ConsensusMaker::addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap)
{
	size_t fromUnitig = std::numeric_limits<size_t>::max();
	size_t fromStartIndex = std::numeric_limits<size_t>::max();
	bool fromFw = true;
	size_t toUnitig = std::numeric_limits<size_t>::max();
	size_t toStartIndex = std::numeric_limits<size_t>::max();
	bool toFw = true;
	size_t matchLength = 0;
	size_t currentInsertUnitig = std::numeric_limits<size_t>::max();
	size_t currentInsertIndex = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < overlap; i++)
	{
		size_t fromIndex = unitigLength(from.first) - overlap + i;
		if (!from.second) fromIndex = overlap - i - 1;
		size_t toIndex = i;
		if (!to.second) toIndex = unitigLength(to.first) - 1 - i;
		auto fromTuple = find(from.first, fromIndex);
		auto toTuple = find(to.first, toIndex);
		assert(getParent(std::get<0>(fromTuple), std::get<1>(fromTuple)) == std::make_tuple(std::get<0>(fromTuple), std::get<1>(fromTuple), true));
		assert(getParent(std::get<0>(toTuple), std::get<1>(toTuple)) == std::make_tuple(std::get<0>(toTuple), std::get<1>(toTuple), true));
		if (!from.second) std::get<2>(fromTuple) = !std::get<2>(fromTuple);
		if (!to.second) std::get<2>(toTuple) = !std::get<2>(toTuple);
		size_t newFromUnitig = std::get<0>(fromTuple);
		size_t newFromPos = std::get<1>(fromTuple);
		bool newFromFw = std::get<2>(fromTuple);
		size_t newToUnitig = std::get<0>(toTuple);
		size_t newToPos = std::get<1>(toTuple);
		bool newToFw = std::get<2>(toTuple);
		if (fromUnitig == newFromUnitig && toUnitig == newToUnitig && (newFromFw == fromFw) && (newToFw == toFw))
		{
			if (newFromPos == fromStartIndex + (fromFw ? 1 : -1) * matchLength && newToPos == toStartIndex + (toFw ? 1 : -1) * matchLength)
			{
				if (newFromUnitig != newToUnitig || newFromPos != newToPos)
				{
					assert(currentInsertUnitig < parent.size());
					assert(currentInsertIndex < parent[currentInsertUnitig].size());
					std::get<1>(parent[currentInsertUnitig][currentInsertIndex]) += 1;
					if (!fromFw) std::get<0>(parent[currentInsertUnitig][currentInsertIndex]) -= 1;
					if (!toFw) std::get<3>(parent[currentInsertUnitig][currentInsertIndex]) -= 1;
					matchLength += 1;
					continue;
				}
			}
		}
		fromUnitig = newFromUnitig;
		fromStartIndex = newFromPos;
		fromFw = newFromFw;
		toUnitig = newToUnitig;
		toStartIndex = newToPos;
		toFw = newToFw;
		matchLength = 1;
		if (fromUnitig == toUnitig && fromStartIndex == toStartIndex)
		{
			matchLength = 0;
			fromUnitig = std::numeric_limits<size_t>::max();
			toUnitig = std::numeric_limits<size_t>::max();
			currentInsertUnitig = std::numeric_limits<size_t>::max();
		}
		if (matchLength >= 1)
		{
			size_t fromRealStart = fromStartIndex;
			if (!fromFw) fromRealStart -= matchLength - 1;
			size_t toRealStart = toStartIndex;
			if (!toFw) toRealStart -= matchLength - 1;
			assert(fromRealStart < unitigLength(fromUnitig));
			assert(toRealStart < unitigLength(toUnitig));
			assert(fromRealStart + matchLength <= unitigLength(fromUnitig));
			assert(toRealStart + matchLength <= unitigLength(toUnitig));
			parent[fromUnitig].emplace_back(fromRealStart, matchLength, toUnitig, toRealStart, fromFw == toFw);
			std::sort(parent[fromUnitig].begin(), parent[fromUnitig].end(), [](auto left, auto right) { return std::get<0>(left) < std::get<0>(right); });
			currentInsertUnitig = fromUnitig;
			currentInsertIndex = std::numeric_limits<size_t>::max();
			for (size_t j = 0; j < parent[fromUnitig].size(); j++)
			{
				if (std::get<0>(parent[fromUnitig][j]) == fromRealStart)
				{
					currentInsertIndex = j;
					break;
				}
			}
			assert(currentInsertIndex != std::numeric_limits<size_t>::max());
		}
	}
	for (size_t i = 0; i < overlap; i++)
	{
		size_t fromIndex = unitigLength(from.first) - overlap + i;
		if (!from.second) fromIndex = overlap - i - 1;
		size_t toIndex = i;
		if (!to.second) toIndex = unitigLength(to.first) - 1 - i;
		auto fromTuple = find(from.first, fromIndex);
		auto toTuple = find(to.first, toIndex);
		// if (!from.second) std::get<2>(fromTuple) = !std::get<2>(fromTuple);
		// if (!to.second) std::get<2>(toTuple) = !std::get<2>(toTuple);
		assert(std::get<0>(fromTuple) == std::get<0>(toTuple));
		assert(std::get<1>(fromTuple) == std::get<1>(toTuple));
		// assert((std::get<2>(fromTuple) == std::get<2>(toTuple)) == (from.second == to.second));
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
						if (count == maxCount)
						{
							if (std::get<2>(found))
							{
								if (stringIndex.getString(compressed, index) < stringIndex.getString(compressed, maxIndex)) continue;
							}
							else if (compressed == complement(compressed))
							{
								if (stringIndex.getString(compressed, index) < stringIndex.getString(compressed, maxIndex)) continue;
							}
							else
							{
								if (stringIndex.getString(complement(compressed), stringIndex.getReverseIndex(compressed, index)) < stringIndex.getString(complement(compressed), stringIndex.getReverseIndex(compressed, maxIndex))) continue;
							}
						}
						maxIndex = index;
						maxCount = count;
					}
				}
			}
			if (!std::get<2>(found))
			{
				if (compressed == complement(compressed))
				{
					maxIndex = stringIndex.getReverseIndex(compressed, maxIndex);
				}
				else
				{
					assert(maxIndex == stringIndex.getReverseIndex(compressed, maxIndex));
				}
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

std::tuple<size_t, size_t, bool> ConsensusMaker::find(size_t unitig, size_t index)
{
	auto result = getParent(unitig, index);
	size_t iterations = 0;
	while (true)
	{
		iterations += 1;
		auto nextResult = getParent(std::get<0>(result), std::get<1>(result));
		if (std::get<0>(nextResult) == std::get<0>(result) && std::get<1>(nextResult) == std::get<1>(result))
		{
			assert(std::get<2>(nextResult));
			break;
		}
		if (!std::get<2>(result)) std::get<2>(nextResult) = !std::get<2>(nextResult);
		result = nextResult;
		assert(iterations < totalLength);
	}
	return result;
}

std::tuple<size_t, size_t, bool> ConsensusMaker::getParent(size_t unitig, size_t index) const
{
	if (parent[unitig].size() == 0)
	{
		return std::make_tuple(unitig, index, true);
	}
	size_t startIndex = 0;
	size_t endIndex = parent[unitig].size();
	size_t mid = (endIndex+startIndex)/2;
	while (endIndex > startIndex+1)
	{
		size_t midStart = std::get<0>(parent[unitig][mid]);
		if (midStart > index)
		{
			endIndex = mid;
			mid = (endIndex+startIndex)/2;
			continue;
		}
		size_t midEnd = midStart + std::get<1>(parent[unitig][mid]);
		if (midEnd <= index)
		{
			startIndex = mid;
			mid = (endIndex+startIndex)/2;
			continue;
		}
		assert(midStart <= index && midEnd > index);
		break;
	}
	size_t midStart = std::get<0>(parent[unitig][mid]);
	size_t midEnd = midStart + std::get<1>(parent[unitig][mid]);
	if (midEnd <= index || midStart > index)
	{
		return std::make_tuple(unitig, index, true);
	}
	assert(mid == parent[unitig].size()-1 || std::get<0>(parent[unitig][mid+1]) >= midEnd);
	assert(mid == 0 || std::get<0>(parent[unitig][mid-1]) + std::get<1>(parent[unitig][mid-1]) <= midStart);
	size_t otherUnitig = std::get<2>(parent[unitig][mid]);
	size_t otherStart = std::get<3>(parent[unitig][mid]);
	bool fw = std::get<4>(parent[unitig][mid]);
	assert(otherUnitig != unitig || otherStart != midStart || !fw);
	size_t otherOffset;
	if (fw)
	{
		otherOffset = otherStart + (index - midStart);
	}
	else
	{
		otherOffset = otherStart + std::get<1>(parent[unitig][mid]) - 1 - (index - midStart);
	}
	return std::make_tuple(otherUnitig, otherOffset, fw);
}

void ConsensusMaker::findParentLinks()
{
	for (size_t i = 0; i < parent.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>> newParent;
		size_t unitig = std::numeric_limits<size_t>::max();
		size_t startPos = 0;
		size_t matchLength = 0;
		bool fw = true;
		if (longestLeftOverlap[i] + longestRightOverlap[i] < unitigLength(i))
		{
			for (size_t j = 0; j < longestLeftOverlap[i]; j++)
			{
				auto found = find(i, j);
				if (std::get<0>(found) == unitig && std::get<1>(found) == startPos + (fw ? 1 : -1) * matchLength && std::get<2>(found) == fw)
				{
					matchLength += 1;
					continue;
				}
				if (matchLength > 0)
				{
					if (unitig != i || startPos+matchLength != j)
					{
						size_t realStart = startPos;
						if (!fw) realStart -= matchLength - 1;
						newParent.emplace_back(j-matchLength, matchLength, unitig, realStart, fw);
					}
				}
				unitig = std::get<0>(found);
				startPos = std::get<1>(found);
				fw = std::get<2>(found);
				matchLength = 1;
			}
			if (matchLength > 0)
			{
				if (unitig != i || startPos+matchLength != longestLeftOverlap[i])
				{
					size_t realStart = startPos;
					if (!fw) realStart -= matchLength - 1;
					newParent.emplace_back(longestLeftOverlap[i]-matchLength, matchLength, unitig, realStart, fw);
				}
			}
			size_t off = unitigLength(i) - longestRightOverlap[i];
			unitig = std::numeric_limits<size_t>::max();
			startPos = 0;
			fw = true;
			matchLength = 0;
			for (size_t j = 0; j < longestRightOverlap[i]; j++)
			{
				auto found = find(i, off+j);
				if (std::get<0>(found) == unitig && std::get<1>(found) == startPos + (fw ? 1 : -1) * matchLength && std::get<2>(found) == fw)
				{
					matchLength += 1;
					continue;
				}
				if (matchLength > 0)
				{
					if (unitig != i || startPos+matchLength != off+j)
					{
						size_t realStart = startPos;
						if (!fw) realStart -= matchLength - 1;
						newParent.emplace_back(off+j-matchLength, matchLength, unitig, realStart, fw);
					}
				}
				unitig = std::get<0>(found);
				startPos = std::get<1>(found);
				fw = std::get<2>(found);
				matchLength = 1;
			}
			if (matchLength > 0)
			{
				if (unitig != i || startPos+matchLength != off+longestRightOverlap[i])
				{
					size_t realStart = startPos;
					if (!fw) realStart -= matchLength - 1;
					newParent.emplace_back(off+longestRightOverlap[i]-matchLength, matchLength, unitig, realStart, fw);
				}
			}
		}
		else
		{
			for (size_t j = 0; j < unitigLength(i); j++)
			{
				auto found = find(i, j);
				if (std::get<0>(found) == unitig && std::get<1>(found) == startPos + (fw ? 1 : -1) * matchLength && std::get<2>(found) == fw)
				{
					matchLength += 1;
					continue;
				}
				if (matchLength > 0)
				{
					if (unitig != i || startPos+matchLength != j)
					{
						size_t realStart = startPos;
						if (!fw) realStart -= matchLength - 1;
						newParent.emplace_back(j-matchLength, matchLength, unitig, realStart, fw);
					}
				}
				unitig = std::get<0>(found);
				startPos = std::get<1>(found);
				fw = std::get<2>(found);
				matchLength = 1;
			}
			if (matchLength > 0)
			{
				if (unitig != i || startPos+matchLength != unitigLength(i))
				{
					size_t realStart = startPos;
					if (!fw) realStart -= matchLength - 1;
					newParent.emplace_back(unitigLength(i)-matchLength, matchLength, unitig, realStart, fw);
				}
			}
		}
		std::swap(parent[i], newParent);
	}
}

void ConsensusMaker::prepareHpcVariants(const std::vector<bool>& checkUnitig)
{
	for (auto pair : needsComplementVerification)
	{
		auto found = find(pair.first, pair.second);
		if (!checkUnitig[std::get<0>(found)]) continue;
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
		std::vector<std::pair<size_t, size_t>> lengthCountsVec { lengthCounts.begin(), lengthCounts.end() };
		std::sort(lengthCountsVec.begin(), lengthCountsVec.end(), [](const std::pair<size_t, size_t>& left, const std::pair<size_t, size_t>& right) { return left.first < right.first; });
		bool currentHasHighCount = false;
		bool valid = true;
		std::vector<size_t> hpcVariants;
		size_t motifLength = codeMotifLength(compressed);
		for (size_t i = 0; i < lengthCountsVec.size(); i++)
		{
			if (i > 0 && lengthCountsVec[i].first > lengthCountsVec[i-1].first + motifLength)
			{
				if (!currentHasHighCount)
				{
					valid = false;
					break;
				}
				hpcVariants.push_back(lengthCountsVec[i-1].first);
				currentHasHighCount = false;
			}
			if (lengthCountsVec[i].second >= minCoverage) currentHasHighCount = true;
		}
		if (!currentHasHighCount) valid = false;
		if (!valid) continue;
		hpcVariants.push_back(lengthCountsVec.back().first);
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
