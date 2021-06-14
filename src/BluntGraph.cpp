#include <unordered_map>
#include <unordered_set>
#include "BluntGraph.h"
#include "VectorWithDirection.h"
#include "HPCConsensus.h"

template <typename F>
void iterateEdges(const HashList& hashlist, const UnitigGraph& unitigs, F callback)
{
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : unitigs.edges[fw])
		{
			std::pair<size_t, bool> last = unitigs.unitigs[i].back();
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t overlap = hashlist.getOverlap(last, first);
			size_t coverage = unitigs.edgeCoverage(fw, to);
			callback(fw, to, overlap, coverage);
		}
		for (auto to : unitigs.edges[bw])
		{
			std::pair<size_t, bool> last = reverse(unitigs.unitigs[i][0]);
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t overlap = hashlist.getOverlap(last, first);
			size_t coverage = unitigs.edgeCoverage(bw, to);
			callback(bw, to, overlap, coverage);
		}
	}
}

BluntGraph::BluntGraph(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::pair<std::string, std::vector<uint8_t>>>& unitigSequences)
{
	VectorWithDirection<size_t> maxOverlap = getMaxOverlaps(hashlist, unitigs);
	initializeNodes(unitigSequences, maxOverlap, unitigs);
	initializeEdgesAndEdgenodes(hashlist, unitigs, maxOverlap, unitigSequences);
	assert(nodes.size() == nodeAvgCoverage.size());
}

void BluntGraph::initializeEdgesAndEdgenodes(const HashList& hashlist, const UnitigGraph& unitigs, const VectorWithDirection<size_t>& maxOverlap, const std::vector<std::pair<std::string, std::vector<uint8_t>>>& unbluntSequences)
{
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> processedCanons;
	processedCanons.resize(unitigs.unitigs.size());
	iterateEdges(hashlist, unitigs, [this, &maxOverlap, &unbluntSequences, &processedCanons](const std::pair<size_t, bool> from, const std::pair<size_t, bool> to, const size_t overlap, const size_t coverage) {
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
		std::string missingSeq = unbluntSequences[from.first].first;
		std::vector<uint8_t> missingLength = unbluntSequences[from.first].second;
		if (!from.second)
		{
			missingSeq = revCompRLE(missingSeq);
			std::reverse(missingLength.begin(), missingLength.end());
		}
		size_t startIndex = missingSeq.size() - fromClipped;
		assert(startIndex < missingSeq.size());
		assert(startIndex + missingSeqSize <= missingSeq.size());
		assert(missingSeq.size() == missingLength.size());
		missingSeq = missingSeq.substr(startIndex, missingSeqSize);
		missingLength = std::vector<uint8_t> { missingLength.begin() + startIndex, missingLength.begin() + startIndex + missingSeqSize };
		assert(missingSeq.size() == missingLength.size());
		assert(missingSeq.size() > 0);
		size_t newNodeIndex = nodes.size();
		assert(nodeAvgCoverage.size() == newNodeIndex);
		nodes.push_back(getExpandedSequence(missingSeq, missingLength));
		nodeAvgCoverage.push_back(coverage);
		edges.emplace_back(from.first, from.second, newNodeIndex, true, coverage);
		edges.emplace_back(newNodeIndex, true, to.first, to.second, coverage);
	});
}

void BluntGraph::initializeNodes(const std::vector<std::pair<std::string, std::vector<uint8_t>>>& unbluntSequences, const VectorWithDirection<size_t>& maxOverlap, const UnitigGraph& unitigs)
{
	nodes.resize(unitigs.unitigs.size());
	nodeAvgCoverage.resize(unitigs.unitigs.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t leftClip = (maxOverlap[std::make_pair(i, false)] + 1) / 2;
		size_t rightClip = (maxOverlap[std::make_pair(i, true)] + 1) / 2;
		std::string sequence = unbluntSequences[i].first;
		std::vector<uint8_t> lengths = unbluntSequences[i].second;
		sequence = sequence.substr(leftClip);
		sequence = sequence.substr(0, sequence.size() - rightClip);
		lengths = std::vector<uint8_t> { lengths.begin() + leftClip, lengths.end() - rightClip };
		nodes[i] = getExpandedSequence(sequence, lengths);
		nodeAvgCoverage[i] = unitigs.averageCoverage(i);
	}
}

VectorWithDirection<size_t> BluntGraph::getMaxOverlaps(const HashList& hashlist, const UnitigGraph& unitigs) const
{
	VectorWithDirection<size_t> maxOverlap;
	maxOverlap.resize(unitigs.unitigs.size(), 0);
	iterateEdges(hashlist, unitigs, [&maxOverlap](const std::pair<size_t, bool> from, const std::pair<size_t, bool> to, const size_t overlap, const size_t coverage) {
		maxOverlap[from] = std::max(maxOverlap[from], overlap);
		maxOverlap[reverse(to)] = std::max(maxOverlap[reverse(to)], overlap);
	});
	return maxOverlap;
}