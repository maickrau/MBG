#include <algorithm>
#include <thread>
#include <chrono>
#include "HashList.h"

AdjacentMinimizerBucket::AdjacentMinimizerBucket() :
	data(),
	lastHash(0)
{
}

TwobitView AdjacentMinimizerBucket::getView(size_t coord1, size_t coord2, size_t size) const
{
	return TwobitView { data[coord1], coord2, coord2+size };
}

std::pair<size_t, size_t> AdjacentMinimizerBucket::addString(std::string_view str, HashType currentHash, HashType previousHash, size_t overlap)
{
	if (data.size() == 0 || lastHash == 0 || previousHash == 0 || previousHash != lastHash)
	{
		data.emplace_back(str);
		lastHash = currentHash;
		return std::make_pair(data.size()-1, 0);
	}
	assert(overlap < str.size());
	data.back().insert(str.begin() + overlap, str.end());
	lastHash = currentHash;
	assert(data.back().size() >= str.size());
	return std::make_pair(data.size()-1, data.back().size() - str.size());
}

AdjacentMinimizerBucket AdjacentMinimizerBucket::getReverseComplementStorage() const
{
	AdjacentMinimizerBucket result;
	result.data.resize(data.size());
	for (size_t i = 0; i < data.size(); i++)
	{
		result.data[i] = revCompRLE(data[i]);
	}
	return result;
}

std::pair<size_t, size_t> AdjacentMinimizerBucket::getRevCompLocation(size_t coord1, size_t coord2, size_t size) const
{
	assert(data[coord1].size() >= coord2 + size);
	return std::make_pair(coord1, data[coord1].size() - size - coord2);
}

AdjacentLengthBucket::AdjacentLengthBucket() :
	sums(),
	counts(),
	lastHash(0)
{
}

std::vector<uint16_t> AdjacentLengthBucket::getData(size_t coord1, size_t coord2, size_t size) const
{
	std::vector<uint16_t> result;
	result.resize(size, 0);
	assert(sums.size() == counts.size());
	assert(sums[coord1].size() == counts[coord1].size());
	for (size_t i = 0; i < size; i++)
	{
		result[i] = int(double(sums[coord1][coord2+i]) / double(counts[coord1][coord2+i]) + .5);
	}
	return result;
}

void AdjacentLengthBucket::setData(size_t coord1, size_t coord2, size_t value)
{
	sums[coord1][coord2] = value;
	counts[coord1][coord2] = 1;
}

std::pair<size_t, size_t> AdjacentLengthBucket::addData(const std::vector<uint16_t>& lens, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap)
{
	assert(end > start);
	assert(end <= lens.size());
	if (sums.size() == 0 || lastHash == 0 || previousHash == 0 || previousHash != lastHash)
	{
		sums.emplace_back(lens.begin() + start, lens.begin() + end);
		counts.emplace_back();
		counts.back().resize(end - start, 1);
		lastHash = currentHash;
		assert(sums.back().size() == counts.back().size());
		return std::make_pair(sums.size()-1, 0);
	}
	assert(overlap < lens.size());
	assert(end > start + overlap);
	sums.back().insert(sums.back().end(), lens.begin() + start + overlap, lens.begin() + end);
	counts.back().resize(counts.back().size() + end - (start + overlap), 1);
	assert(counts.back().size() == sums.back().size());
	lastHash = currentHash;
	assert(sums.back().size() >= end - start);
	return std::make_pair(sums.size()-1, sums.back().size() - (end - start));
}

void AdjacentLengthBucket::addCounts(const std::vector<uint16_t>& lens, bool fw, size_t start, size_t end, size_t coord1, size_t coord2)
{
	assert(coord1 < sums.size());
	assert(sums.size() == counts.size());
	assert(coord2 < sums[coord1].size());
	assert(sums[coord1].size() == counts[coord1].size());
	assert(coord2 + end - start <= sums[coord1].size());
	if (fw)
	{
		for (size_t i = 0; i < end - start; i++)
		{
			if (lens[start+i] < std::numeric_limits<uint16_t>::max() - sums[coord1][coord2+i] && counts[coord1][coord2 + i] < std::numeric_limits<uint8_t>::max() - 1)
			{
				sums[coord1][coord2 + i] += lens[start + i];
				counts[coord1][coord2 + i] += 1;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < end - start; i++)
		{
			if (lens[end - 1 - i] < std::numeric_limits<uint16_t>::max() - sums[coord1][coord2+i] && counts[coord1][coord2 + i] < std::numeric_limits<uint8_t>::max() - 1)
			{
				sums[coord1][coord2 + i] += lens[end - 1 - i];
				counts[coord1][coord2 + i] += 1;
			}
		}
	}
}

size_t AdjacentLengthBucket::size() const
{
	size_t total = 0;
	for (size_t i = 0; i < sums.size(); i++)
	{
		total += sums[i].size();
	}
	return total;
}

AdjacentMinimizerList::AdjacentMinimizerList(size_t numBuckets) :
buckets(numBuckets),
bucketMutexes(numBuckets)
{
}

TwobitView AdjacentMinimizerList::getView(size_t bucket, std::pair<size_t, size_t> coords, size_t size) const
{
	return buckets[bucket].getView(coords.first, coords.second, size);
}

size_t AdjacentMinimizerList::hashToBucket(HashType hash) const
{
	return hash % buckets.size();
}

std::pair<size_t, std::pair<size_t, size_t>> AdjacentMinimizerList::addString(std::string_view str, HashType currentHash, HashType previousHash, size_t overlap, uint64_t bucketHash)
{
	size_t bucket = hashToBucket(bucketHash);
	std::lock_guard<std::mutex> lock { bucketMutexes[bucket] };
	return std::make_pair(bucket, buckets[bucket].addString(str, currentHash, previousHash, overlap));
}

AdjacentMinimizerList AdjacentMinimizerList::getReverseComplementStorage() const
{
	AdjacentMinimizerList result { buckets.size() };
	for (size_t i = 0; i < buckets.size(); i++)
	{
		result.buckets[i] = buckets[i].getReverseComplementStorage();
	}
	return result;
}

std::pair<size_t, std::pair<size_t, size_t>> AdjacentMinimizerList::getRevCompLocation(size_t bucket, std::pair<size_t, size_t> coords, size_t size) const
{
	return std::make_pair(bucket, buckets[bucket].getRevCompLocation(coords.first, coords.second, size));
}

AdjacentLengthList::AdjacentLengthList(size_t numBuckets) :
buckets(numBuckets),
bucketMutexes(numBuckets)
{
}

std::vector<uint16_t> AdjacentLengthList::getData(size_t bucket, std::pair<size_t, size_t> coords, size_t size) const
{
	return buckets[bucket].getData(coords.first, coords.second, size);
}

void AdjacentLengthList::setData(size_t bucket, std::pair<size_t, size_t> coords, size_t value)
{
	buckets[bucket].setData(coords.first, coords.second, value);
}

std::pair<size_t, std::pair<size_t, size_t>> AdjacentLengthList::addData(const std::vector<uint16_t>& lens, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap, uint64_t bucketHash)
{
	size_t bucket = hashToBucket(bucketHash);
	std::lock_guard<std::mutex> lock { bucketMutexes[bucket] };
	return std::make_pair(bucket, buckets[bucket].addData(lens, start, end, currentHash, previousHash, overlap));
}

void AdjacentLengthList::addCounts(const std::vector<uint16_t>& lens, bool fw, size_t start, size_t end, size_t bucket, std::pair<size_t, size_t> coords)
{
	std::lock_guard<std::mutex> lock { bucketMutexes[bucket] };
	buckets[bucket].addCounts(lens, fw, start, end, coords.first, coords.second);
}

size_t AdjacentLengthList::size() const
{
	size_t sum = 0;
	for (size_t i = 0; i < buckets.size(); i++)
	{
		sum += buckets[i].size();
	}
	return sum;
}

size_t AdjacentLengthList::hashToBucket(HashType hash) const
{
	return hash % buckets.size();
}

HashList::HashList(size_t kmerSize, bool collapseRunLengths, size_t numBuckets) :
	hashCharacterLengths(numBuckets),
	hashSequences(numBuckets),
	hashSequencesRevComp(numBuckets),
	kmerSize(kmerSize),
	collapseRunLengths(collapseRunLengths)
{}

size_t HashList::numSequenceOverlaps() const
{
	size_t total = 0;
	for (size_t i = 0; i < sequenceOverlap.size(); i++)
	{
		total += sequenceOverlap[std::make_pair(i, true)].size();
		total += sequenceOverlap[std::make_pair(i, false)].size();
	}
	return total;
}

size_t HashList::getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	assert(edgeCoverage.at(from).count(to) == 1);
	return edgeCoverage.at(from).at(to);
}

size_t HashList::getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return sequenceOverlap.at(from).at(to);
}

void HashList::addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap)
{
	std::tie(from, to) = canon(from, to);
	std::lock_guard<std::mutex> lock { indexMutex };
	if (sequenceOverlap[from].count(to) == 1) return;
	sequenceOverlap[from][to] = overlap;
}

size_t HashList::size() const
{
	return hashSeqPtr.size();
}

std::vector<uint16_t> HashList::getHashCharacterLength(size_t index) const
{
	if (collapseRunLengths)
	{
		std::vector<uint16_t> result;
		result.resize(kmerSize, 1);
		return result;
	}
	return hashCharacterLengths.getData(hashCharacterLengthPtr[index].first, hashCharacterLengthPtr[index].second, kmerSize);
}

void HashList::addHashCharacterLength(const std::vector<uint16_t>& data, bool fw, size_t start, size_t end, size_t node, std::pair<size_t, std::pair<size_t, size_t>> position)
{
	if (collapseRunLengths) return;
	hashCharacterLengths.addCounts(data, fw, start, end, position.first, position.second);
}

TwobitView HashList::getHashSequenceRLE(size_t index) const
{
	return hashSequences.getView(hashSeqPtr[index].first, hashSeqPtr[index].second, kmerSize);
}

TwobitView HashList::getRevCompHashSequenceRLE(size_t index) const
{
	auto pos = hashSequences.getRevCompLocation(hashSeqPtr[index].first, hashSeqPtr[index].second, kmerSize);
	return hashSequencesRevComp.getView(pos.first, pos.second, kmerSize);
}

void HashList::buildReverseCompHashSequences()
{
	hashSequencesRevComp = hashSequences.getReverseComplementStorage();
}

std::pair<size_t, bool> HashList::getNodeOrNull(std::string_view sequence) const
{
	HashType fwHash = hash(sequence);
	auto found = hashToNode.find(fwHash);
	if (found == hashToNode.end())
	{
		return std::pair<size_t, bool> { std::numeric_limits<size_t>::max(), true };
	}
	return found->second;
}

void HashList::addEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	std::lock_guard<std::mutex> lock { indexMutex };
	edgeCoverage[from][to] += 1;
}

std::pair<std::pair<size_t, bool>, HashType> HashList::addNode(std::string_view sequence, std::string_view reverse, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap, uint64_t bucketHash)
{
	HashType fwHash = hash(sequence);
	start_without_sleep:
	if (false)
	{
		start_with_sleep: // ehhh just a bit gnarly
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
	std::pair<size_t, bool> index { std::numeric_limits<size_t>::max(), false };
	std::pair<size_t, std::pair<size_t, size_t>> position;
	size_t fwNode;
	{
		std::lock_guard<std::mutex> lock { indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			index = found->second;
			if (!collapseRunLengths)
			{
				position = hashCharacterLengthPtr[index.first];
				if (position.first == std::numeric_limits<size_t>::max())
				{
					goto start_with_sleep;
				}
			}
			coverage[found->second.first] += 1;
		}
	}
	if (index.first != std::numeric_limits<size_t>::max())
	{
		if (!collapseRunLengths) addHashCharacterLength(sequenceCharacterLength, index.second, seqCharLenStart, seqCharLenEnd, index.first, position);
		return std::make_pair(index, fwHash);
	}
	{
		std::lock_guard<std::mutex> lock { indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end()) goto start_without_sleep;
		assert(found == hashToNode.end());
		HashType bwHash = hash(reverse);
		found = hashToNode.find(bwHash);
		assert(hashToNode.find(bwHash) == hashToNode.end());
		fwNode = size();
		hashToNode[fwHash] = std::make_pair(fwNode, true);
		hashToNode[bwHash] = std::make_pair(fwNode, false);
		assert(coverage.size() == fwNode);
		assert(edgeCoverage.size() == fwNode);
		assert(sequenceOverlap.size() == fwNode);
		assert(hashSeqPtr.size() == fwNode);
		assert(collapseRunLengths || hashCharacterLengthPtr.size() == fwNode);
		coverage.emplace_back(1);
		edgeCoverage.emplace_back();
		sequenceOverlap.emplace_back();
		hashSeqPtr.emplace_back();
		if (!collapseRunLengths)
		{
			hashCharacterLengthPtr.emplace_back();
			hashCharacterLengthPtr.back().first = std::numeric_limits<size_t>::max();
		}
	}
	std::pair<size_t, std::pair<size_t, size_t>> seqResult;
	std::pair<size_t, std::pair<size_t, size_t>> lengthResult;
	seqResult = hashSequences.addString(sequence, fwHash, previousHash, overlap, bucketHash);
	if (!collapseRunLengths)
	{
		lengthResult = hashCharacterLengths.addData(sequenceCharacterLength, seqCharLenStart, seqCharLenEnd, fwHash, previousHash, overlap, bucketHash);
	}
	{
		std::lock_guard<std::mutex> lock { indexMutex };
		hashSeqPtr[fwNode] = seqResult;
		if (!collapseRunLengths) hashCharacterLengthPtr[fwNode] = lengthResult;
	}
	return std::make_pair(std::make_pair(fwNode, true), fwHash);
}

size_t HashList::getRunLength(size_t index, size_t offset) const
{
	auto coords = hashCharacterLengthPtr[index].second;
	coords.second += offset;
	return hashCharacterLengths.getData(hashCharacterLengthPtr[index].first, coords, 1)[0];
}

void HashList::setRunLength(size_t index, size_t offset, size_t length)
{
	auto coords = hashCharacterLengthPtr[index].second;
	coords.second += offset;
	hashCharacterLengths.setData(hashCharacterLengthPtr[index].first, coords, length);
}