#include <algorithm>
#include <thread>
#include <chrono>
#include "HashList.h"

HashList::HashList(size_t kmerSize) :
	kmerSize(kmerSize)
{
	indexMutex = std::make_shared<std::mutex>();
}

size_t HashList::numSequenceOverlaps() const
{
	size_t total = 0;
	for (size_t i = 0; i < sequenceOverlap.size(); i++)
	{
		total += sequenceOverlap.getValues(std::make_pair(i, true)).size();
		total += sequenceOverlap.getValues(std::make_pair(i, false)).size();
	}
	return total;
}

size_t HashList::getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return edgeCoverage.get(from, to);
}

std::vector<std::pair<std::pair<size_t, bool>, size_t>> HashList::getEdgeCoverages(std::pair<size_t, bool> from) const
{
	return edgeCoverage.getValues(from);
}

void HashList::setEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t coverage)
{
	std::tie(from, to) = canon(from, to);
	edgeCoverage.set(from, to, coverage);
}

size_t HashList::getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return sequenceOverlap.get(from, to);
}

bool HashList::hasSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return sequenceOverlap.hasValue(from, to);
}

std::vector<std::pair<std::pair<size_t, bool>, size_t>> HashList::getSequenceOverlaps(std::pair<size_t, bool> from) const
{
	return sequenceOverlap.getValues(from);
}

void HashList::addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap)
{
	std::tie(from, to) = canon(from, to);
	std::lock_guard<std::mutex> lock { *indexMutex };
	sequenceOverlap.set(from, to, overlap);
}

size_t HashList::size() const
{
	return coverage.size();
}

std::pair<size_t, bool> HashList::getNodeOrNull(VectorView<CharType> sequence) const
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
	std::lock_guard<std::mutex> lock { *indexMutex };
	if (!edgeCoverage.hasValue(from, to))
	{
		edgeCoverage.set(from, to, 1);
		return;
	}
	edgeCoverage.set(from, to, edgeCoverage.get(from, to) + 1);
}

std::pair<std::pair<size_t, bool>, HashType> HashList::addNode(VectorView<CharType> sequence, VectorView<CharType> reverse, HashType previousHash, size_t overlap, uint64_t bucketHash)
{
	HashType fwHash = hash(sequence);
	{
		std::lock_guard<std::mutex> lock { *indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			coverage.set(found->second.first, coverage.get(found->second.first)+1);
			return std::make_pair(found->second, fwHash);
		}
		assert(found == hashToNode.end());
	}
	HashType bwHash = hash(reverse);
	{
		std::lock_guard<std::mutex> lock { *indexMutex };
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			coverage.set(found->second.first, coverage.get(found->second.first)+1);
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

void HashList::resize(size_t size)
{
	coverage.resize(size, 0);
	sequenceOverlap.resize(size);
	edgeCoverage.resize(size);
}

void HashList::filter(const RankBitvector& kept)
{
	assert(kept.size() == size());
	if (kept.getRank(kept.size()-1) == size()-1) return;
	{
		LittleBigVector<uint8_t, size_t> newCoverage;
		newCoverage.resize(kept.size());
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			newCoverage.set(kept.getRank(i), coverage.get(i));
		}
		std::swap(coverage, newCoverage);
	}
	{
		phmap::flat_hash_map<HashType, std::pair<size_t, bool>> newHashToNode;
		for (auto pair : hashToNode)
		{
			if (!kept.get(pair.second.first)) continue;
			newHashToNode[pair.first] = std::make_pair(kept.getRank(pair.second.first), pair.second.second);
		}
		std::swap(hashToNode, newHashToNode);
	}
	{
		MostlySparse2DHashmap<uint8_t, size_t> newEdgeCoverage;
		newEdgeCoverage.resize(kept.size());
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			for (auto key : edgeCoverage.getValues(std::make_pair(i, true)))
			{
				if (!kept.get(key.first.first)) continue;
				newEdgeCoverage.set(std::make_pair(kept.getRank(i), true), std::make_pair(kept.getRank(key.first.first), key.first.second), key.second);
			}
			for (auto key : edgeCoverage.getValues(std::make_pair(i, false)))
			{
				if (!kept.get(key.first.first)) continue;
				newEdgeCoverage.set(std::make_pair(kept.getRank(i), false), std::make_pair(kept.getRank(key.first.first), key.first.second), key.second);
			}
		}
		std::swap(edgeCoverage, newEdgeCoverage);
	}
	{
		MostlySparse2DHashmap<uint16_t, size_t> newSequenceOverlap;
		newSequenceOverlap.resize(kept.size());
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			for (auto key : sequenceOverlap.getValues(std::make_pair(i, true)))
			{
				if (!kept.get(key.first.first)) continue;
				newSequenceOverlap.set(std::make_pair(kept.getRank(i), true), std::make_pair(kept.getRank(key.first.first), key.first.second), key.second);
			}
			for (auto key : sequenceOverlap.getValues(std::make_pair(i, false)))
			{
				if (!kept.get(key.first.first)) continue;
				newSequenceOverlap.set(std::make_pair(kept.getRank(i), false), std::make_pair(kept.getRank(key.first.first), key.first.second), key.second);
			}
		}
		std::swap(sequenceOverlap, newSequenceOverlap);
	}
}
