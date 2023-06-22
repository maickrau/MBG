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
	return getNodeOrNull(fwHash);
}

std::pair<size_t, bool> HashList::getNodeOrNull(HashType fwHash) const
{
	HashType bwHash = (fwHash << 64) + (fwHash >> 64);
	HashType canonHash = std::min(fwHash, bwHash);
	if (fwHash == bwHash)
	{
		throw PalindromicKmer {};
	}
	assert(fwHash != bwHash);
	bool fw = fwHash < bwHash;
	auto found = hashToNode.find(canonHash);
	if (found != hashToNode.end())
	{
		return std::make_pair(found->second, fw);
	}
	return std::pair<size_t, bool> { std::numeric_limits<size_t>::max(), true };
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

void HashList::addForcedEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	std::lock_guard<std::mutex> lock { *indexMutex };
	forcedEdges.addEdge(from, to);
	forcedEdges.addEdge(reverse(to), reverse(from));
}

std::pair<std::pair<size_t, bool>, HashType> HashList::addNode(VectorView<CharType> sequence, VectorView<CharType> reverse)
{
	HashType fwHash = hash(sequence, reverse);
	return std::make_pair(addNode(fwHash), fwHash);
}

std::pair<size_t, bool> HashList::addNode(HashType fwHash)
{
	HashType bwHash = (fwHash << 64) + (fwHash >> 64);
	// this is a true assertion but commented out just for performance
	// assert(bwHash == hash(reverse));
	HashType canonHash = std::min(fwHash, bwHash);
	assert(fwHash != bwHash);
	bool fw = fwHash < bwHash;
	{
		std::lock_guard<std::mutex> lock { *indexMutex };
		auto found = hashToNode.find(canonHash);
		if (found != hashToNode.end())
		{
			coverage.set(found->second, coverage.get(found->second)+1);
			auto node = std::make_pair(found->second, fw);
			return node;
		}
		assert(found == hashToNode.end());
		size_t fwNode = size();
		hashToNode[canonHash] = fwNode;
		assert(coverage.size() == fwNode);
		assert(edgeCoverage.size() == fwNode);
		assert(sequenceOverlap.size() == fwNode);
		coverage.emplace_back(1);
		edgeCoverage.emplace_back();
		sequenceOverlap.emplace_back();
		forcedEdges.emplace_back();
		forced.emplace_back(false);
		return std::make_pair(fwNode, fw);
	}
}

std::pair<size_t, bool> HashList::addNodeForced(HashType fwHash)
{
	auto result = addNode(fwHash);
	forced[result.first] = true;
	return result;
}

void HashList::resize(size_t size)
{
	coverage.resize(size, 0);
	sequenceOverlap.resize(size);
	edgeCoverage.resize(size);
	forcedEdges.resize(size);
	forced.resize(size, false);
}

void HashList::filter(const RankBitvector& kept)
{
	if (kept.size() == 0) return;
	assert(kept.size() == size());
	size_t newSize = kept.getRank(kept.size()-1) + (kept.get(kept.size()-1) ? 1 : 0);
	if (newSize == size()) return;
	{
		LittleBigVector<uint8_t, size_t> newCoverage;
		newCoverage.resize(newSize);
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			newCoverage.set(kept.getRank(i), coverage.get(i));
		}
		std::swap(coverage, newCoverage);
	}
	{
		SparseEdgeContainer newForcedEdges;
		newForcedEdges.resize(newSize);
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			for (auto edge : forcedEdges.getEdges(std::make_pair(i, true)))
			{
				assert(kept.get(edge.first));
				newForcedEdges.addEdge(std::make_pair(kept.getRank(i), true), std::make_pair(kept.getRank(edge.first), edge.second));
			}
			for (auto edge : forcedEdges.getEdges(std::make_pair(i, false)))
			{
				assert(kept.get(edge.first));
				newForcedEdges.addEdge(std::make_pair(kept.getRank(i), false), std::make_pair(kept.getRank(edge.first), edge.second));
			}
		}
		std::swap(forcedEdges, newForcedEdges);
	}
	{
		std::vector<bool> newForced;
		newForced.resize(newSize, false);
		for (size_t i = 0; i < kept.size(); i++)
		{
			if (!kept.get(i)) continue;
			if (!forced[i]) continue;
			newForced[kept.getRank(i)] = true;
		}
		std::swap(forced, newForced);
	}
	{
		phmap::flat_hash_map<HashType, size_t> newHashToNode;
		for (auto pair : hashToNode)
		{
			if (!kept.get(pair.second)) continue;
			newHashToNode[pair.first] = kept.getRank(pair.second);
		}
		std::swap(hashToNode, newHashToNode);
	}
	{
		MostlySparse2DHashmap<uint8_t, size_t> newEdgeCoverage;
		newEdgeCoverage.resize(newSize);
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
		newSequenceOverlap.resize(newSize);
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

std::pair<size_t, bool> HashList::getHashNode(HashType fwHash) const
{
	HashType bwHash = (fwHash << 64) + (fwHash >> 64);
	HashType canonHash = std::min(fwHash, bwHash);
	assert(fwHash != bwHash);
	bool fw = fwHash < bwHash;
	auto node = hashToNode.at(canonHash);
	auto result = std::make_pair(node, fw);
	return result;
}

std::vector<size_t> HashList::sortByHash()
{
	std::vector<HashType> hashes;
	hashes.reserve(hashToNode.size());
	for (auto pair : hashToNode)
	{
		hashes.emplace_back(pair.first);
	}
	std::sort(hashes.begin(), hashes.end());
	std::vector<size_t> mapping;
	mapping.resize(hashes.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < hashes.size(); i++)
	{
		assert(mapping[hashToNode.at(hashes[i])] == std::numeric_limits<size_t>::max());
		mapping[hashToNode.at(hashes[i])] = i;
	}
	for (size_t i = 0; i < hashes.size(); i++)
	{
		assert(mapping[i] != std::numeric_limits<size_t>::max());
		assert(mapping[i] < hashes.size());
	}
	{
		phmap::flat_hash_map<HashType, size_t> newHashToNode;
		for (size_t i = 0; i < hashes.size(); i++)
		{
			newHashToNode[hashes[i]] = i;
		}
		std::swap(hashToNode, newHashToNode);
	}
	{
		LittleBigVector<uint8_t, size_t> newCoverage;
		newCoverage.resize(coverage.size());
		for (size_t i = 0; i < coverage.size(); i++)
		{
			newCoverage.set(mapping[i], coverage.get(i));
		}
		std::swap(coverage, newCoverage);
	}
	{
		MostlySparse2DHashmap<uint8_t, size_t> newEdgeCoverage;
		newEdgeCoverage.resize(edgeCoverage.size());
		for (size_t i = 0; i < edgeCoverage.size(); i++)
		{
			auto values = edgeCoverage.getValues(std::make_pair(i, true));
			for (auto value : values)
			{
				std::pair<size_t, bool> from = std::make_pair(i, true);
				std::pair<size_t, bool> to = value.first;
				from.first = mapping[from.first];
				to.first = mapping[to.first];
				size_t coverage = value.second;
				auto key = canon(from, to);
				newEdgeCoverage.set(key.first, key.second, coverage);
			}
			values = edgeCoverage.getValues(std::make_pair(i, false));
			for (auto value : values)
			{
				std::pair<size_t, bool> from = std::make_pair(i, false);
				std::pair<size_t, bool> to = value.first;
				from.first = mapping[from.first];
				to.first = mapping[to.first];
				size_t coverage = value.second;
				auto key = canon(from, to);
				newEdgeCoverage.set(key.first, key.second, coverage);
			}
		}
		std::swap(edgeCoverage, newEdgeCoverage);
	}
	{
		MostlySparse2DHashmap<uint16_t, size_t> newSequenceOverlap;
		newSequenceOverlap.resize(sequenceOverlap.size());
		for (size_t i = 0; i < sequenceOverlap.size(); i++)
		{
			auto values = sequenceOverlap.getValues(std::make_pair(i, true));
			for (auto value : values)
			{
				std::pair<size_t, bool> from = std::make_pair(i, true);
				std::pair<size_t, bool> to = value.first;
				from.first = mapping[from.first];
				to.first = mapping[to.first];
				size_t overlap = value.second;
				auto key = canon(from, to);
				newSequenceOverlap.set(key.first, key.second, overlap);
			}
			values = sequenceOverlap.getValues(std::make_pair(i, false));
			for (auto value : values)
			{
				std::pair<size_t, bool> from = std::make_pair(i, false);
				std::pair<size_t, bool> to = value.first;
				from.first = mapping[from.first];
				to.first = mapping[to.first];
				size_t overlap = value.second;
				auto key = canon(from, to);
				newSequenceOverlap.set(key.first, key.second, overlap);
			}
		}
		std::swap(sequenceOverlap, newSequenceOverlap);
	}
	return mapping;
}

void HashList::clear()
{
	decltype(hashToNode) tmp;
	decltype(edgeCoverage) tmp2;
	decltype(sequenceOverlap) tmp3;
	decltype(coverage) tmp4;
	decltype(forced) tmp5;
	decltype(forcedEdges) tmp6;
	std::swap(hashToNode, tmp);
	std::swap(edgeCoverage, tmp2);
	std::swap(sequenceOverlap, tmp3);
	std::swap(coverage, tmp4);
	std::swap(forced, tmp5);
	std::swap(forcedEdges, tmp6);
}
