#include <unordered_map>
#include <unordered_set>
#include "BluntGraph.h"
#include "VectorWithDirection.h"
#include "HPCConsensus.h"
#include "ErrorMaskHelper.h"
#include "UnitigHelper.h"

template <typename F>
void iterateEdges(const HashList& hashlist, const UnitigGraph& unitigs, const size_t kmerSize, F callback)
{
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : unitigs.edges[fw])
		{
			size_t overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, fw, to);
			size_t coverage = unitigs.edgeCoverage(fw, to);
			callback(fw, to, overlap, coverage);
		}
		for (auto to : unitigs.edges[bw])
		{
			size_t overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, bw, to);
			size_t coverage = unitigs.edgeCoverage(bw, to);
			callback(bw, to, overlap, coverage);
		}
	}
}

BluntGraph::BluntGraph(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const size_t kmerSize)
{
	VectorWithDirection<size_t> maxOverlap = getMaxOverlaps(hashlist, unitigs, kmerSize);
	initializeNodes(unitigSequences, stringIndex, maxOverlap, unitigs);
	initializeEdgesAndEdgenodes(hashlist, unitigs, maxOverlap, unitigSequences, stringIndex, kmerSize);
	assert(nodes.size() == nodeAvgCoverage.size());
}

void BluntGraph::initializeEdgesAndEdgenodes(const HashList& hashlist, const UnitigGraph& unitigs, const VectorWithDirection<size_t>& maxOverlap, const std::vector<CompressedSequenceType>& unbluntSequences, const StringIndex& stringIndex, const size_t kmerSize)
{
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> processedCanons;
	processedCanons.resize(unitigs.unitigs.size());
	iterateEdges(hashlist, unitigs, kmerSize, [this, &maxOverlap, &unbluntSequences, &stringIndex, &processedCanons](const std::pair<size_t, bool> from, const std::pair<size_t, bool> to, const size_t overlap, const size_t coverage) {
		auto canonDir = canon(from, to);
		if (processedCanons[canonDir.first].count(canonDir.second) == 1) return;
		processedCanons[canonDir.first].insert(canonDir.second);
		size_t fromClipped = (maxOverlap[from] + 1) / 2;
		size_t toClipped = (maxOverlap[reverse(to)] + 1) / 2;
		if (fromClipped + toClipped == overlap)
		{
			edges.emplace_back(from.first, from.second, to.first, to.second, coverage);
			return;
		}
		assert(fromClipped + toClipped > overlap);
		size_t missingSeqSize = fromClipped + toClipped - overlap;
		CompressedSequenceType missingSeq = unbluntSequences[from.first];
		if (!from.second)
		{
			missingSeq = missingSeq.revComp(stringIndex);
		}
		size_t startIndex = missingSeq.compressedSize() - fromClipped;
		assert(startIndex < missingSeq.compressedSize());
		assert(startIndex + missingSeqSize <= missingSeq.compressedSize());
		missingSeq = missingSeq.substr(startIndex, missingSeqSize);
		assert(missingSeq.compressedSize() > 0);
		assert(missingSeq.compressedSize() == missingSeqSize);
		size_t newNodeIndex = nodes.size();
		assert(nodeAvgCoverage.size() == newNodeIndex);
		nodes.push_back(missingSeq.getExpandedSequence(stringIndex));
		nodeAvgCoverage.push_back(coverage);
		edges.emplace_back(from.first, from.second, newNodeIndex, true, coverage);
		edges.emplace_back(newNodeIndex, true, to.first, to.second, coverage);
	});
}

void BluntGraph::initializeNodes(const std::vector<CompressedSequenceType>& unbluntSequences, const StringIndex& stringIndex, const VectorWithDirection<size_t>& maxOverlap, const UnitigGraph& unitigs)
{
	nodes.resize(unitigs.unitigs.size());
	nodeAvgCoverage.resize(unitigs.unitigs.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t leftClip = (maxOverlap[std::make_pair(i, false)] + 1) / 2;
		size_t rightClip = (maxOverlap[std::make_pair(i, true)] + 1) / 2;
		assert(leftClip + rightClip < unbluntSequences[i].compressedSize());
		CompressedSequenceType sequence = unbluntSequences[i].substr(leftClip, unbluntSequences[i].compressedSize() - leftClip - rightClip);
		assert(sequence.compressedSize() >= 1);
		nodes[i] = sequence.getExpandedSequence(stringIndex);
		nodeAvgCoverage[i] = unitigs.averageCoverage(i);
	}
}

VectorWithDirection<size_t> BluntGraph::getMaxOverlaps(const HashList& hashlist, const UnitigGraph& unitigs, const size_t kmerSize) const
{
	VectorWithDirection<size_t> maxOverlap;
	maxOverlap.resize(unitigs.unitigs.size(), 0);
	iterateEdges(hashlist, unitigs, kmerSize, [&maxOverlap](const std::pair<size_t, bool> from, const std::pair<size_t, bool> to, const size_t overlap, const size_t coverage) {
		maxOverlap[from] = std::max(maxOverlap[from], overlap);
		maxOverlap[reverse(to)] = std::max(maxOverlap[reverse(to)], overlap);
	});
	return maxOverlap;
}