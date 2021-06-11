#include <algorithm>
#include <thread>
#include <chrono>
#include "HashList.h"

HashList::HashList(size_t kmerSize) :
	kmerSize(kmerSize)
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
	return coverage.size();
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
	{
		std::lock_guard<std::mutex> lock { indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			coverage[found->second.first] += 1;
			return std::make_pair(found->second, fwHash);
		}
		assert(found == hashToNode.end());
	}
	HashType bwHash = hash(reverse);
	{
		std::lock_guard<std::mutex> lock { indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			coverage[found->second.first] += 1;
			return std::make_pair(found->second, fwHash);
		}
		assert(found == hashToNode.end());
		found = hashToNode.find(bwHash);
		assert(hashToNode.find(bwHash) == hashToNode.end());
		size_t fwNode = size();
		hashToNode[fwHash] = std::make_pair(fwNode, true);
		hashToNode[bwHash] = std::make_pair(fwNode, false);
		assert(coverage.size() == fwNode);
		assert(edgeCoverage.size() == fwNode);
		assert(sequenceOverlap.size() == fwNode);
		coverage.emplace_back(1);
		edgeCoverage.emplace_back();
		sequenceOverlap.emplace_back();
		return std::make_pair(std::make_pair(fwNode, true), fwHash);
	}
}
