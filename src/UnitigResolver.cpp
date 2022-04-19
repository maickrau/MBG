#include <iostream>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <phmap.h>
#include "MBGCommon.h"
#include "UnitigResolver.h"
#include "RankBitvector.h"
#include "ReadHelper.h"
#include "BigVectorSet.h"

#define assertPrintReads(expression, graph, paths, node) {if (!(expression)) {printReads(graph, paths, node);} assert(expression);}

class PathGroup
{
public:
	class Read
	{
	public:
		size_t readNameIndex;
		std::vector<size_t> readPoses;
		size_t expandedReadPosStart;
		size_t expandedReadPosEnd;
		size_t leftClip;
		size_t rightClip;
		size_t readLength;
		size_t readLengthHPC;
	};
	std::vector<std::pair<size_t, bool>> path;
	std::vector<Read> reads;
};

class ResolvableUnitigGraph
{
public:
	ResolvableUnitigGraph(const HashList& hashlist, const size_t kmerSize) :
	lastTippableChecked(0),
	hashlist(hashlist),
	kmerSize(kmerSize)
	{
	}
	std::vector<std::vector<std::pair<size_t, bool>>> unitigs;
	std::vector<size_t> unitigLeftClipBp;
	std::vector<size_t> unitigRightClipBp;
	VectorWithDirection<phmap::flat_hash_set<std::pair<size_t, bool>>> edges;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> overlaps;
	std::vector<bool> unitigRemoved;
	std::vector<phmap::flat_hash_set<uint32_t>> readsCrossingNode;
	std::vector<size_t> everTippable;
	size_t lastTippableChecked;
	mutable std::vector<size_t> precalcedUnitigLengths;
	size_t getBpOverlap(const std::pair<size_t, bool> from, const std::pair<size_t, bool> to) const
	{
		size_t kmerOverlap = overlaps.at(canon(from, to));
		size_t result;
		if (kmerOverlap == 0)
		{
			std::pair<size_t, bool> fromKmer = from.second ? unitigs[from.first].back() : reverse(unitigs[from.first][0]);
			std::pair<size_t, bool> toKmer = to.second ? unitigs[to.first][0] : reverse(unitigs[to.first].back());
			result = hashlist.getOverlap(fromKmer, toKmer);
		}
		else
		{
			result = kmerOverlap * kmerSize;
			for (size_t i = 1; i < kmerOverlap; i++)
			{
				std::pair<size_t, bool> oldKmer = unitigs[to.first][i-1];
				std::pair<size_t, bool> newKmer = unitigs[to.first][i];
				if (!to.second)
				{
					oldKmer = reverse(unitigs[to.first][unitigs[to.first].size()-i]);
					newKmer = reverse(unitigs[to.first][unitigs[to.first].size()-i-1]);
				}
				result -= hashlist.getOverlap(oldKmer, newKmer);
			}
		}
		size_t fromClip = from.second ? unitigRightClipBp[from.first] : unitigLeftClipBp[from.first];
		size_t toClip = to.second ? unitigLeftClipBp[to.first] : unitigRightClipBp[to.first];
		assert(result >= fromClip + toClip);
		result -= fromClip + toClip;
		return result;
	}
	size_t unitigLength(size_t i) const
	{
		if (i >= precalcedUnitigLengths.size())
		{
			assert(precalcedUnitigLengths.size() < unitigs.size());
			precalcedUnitigLengths.resize(unitigs.size(), 0);
		}
		assert(i < precalcedUnitigLengths.size());
		if (precalcedUnitigLengths[i] == 0)
		{
			size_t result = unitigs[i].size() * kmerSize;
			for (size_t j = 1; j < unitigs[i].size(); j++)
			{
				assert(result > hashlist.getOverlap(unitigs[i][j-1], unitigs[i][j]));
				result -= hashlist.getOverlap(unitigs[i][j-1], unitigs[i][j]);
			}
			assert(result >= kmerSize);
			assert(result > unitigLeftClipBp[i] + unitigRightClipBp[i]);
			result -= unitigLeftClipBp[i] + unitigRightClipBp[i];
			assert(result >= kmerSize);
			precalcedUnitigLengths[i] = result;
		}
		assert(precalcedUnitigLengths[i] >= kmerSize);
		return precalcedUnitigLengths[i];
	}
	const HashList& hashlist;
	std::vector<std::string> readNames;
	double averageCoverage;
private:
	const size_t kmerSize;
};

void printReads(const ResolvableUnitigGraph& graph, const std::vector<PathGroup>& paths, size_t node)
{
	std::unordered_set<size_t> closeReads;
	std::unordered_set<size_t> semiCloseReads;
	std::cerr << "around center node " << node << std::endl;
	for (auto readi : graph.readsCrossingNode[node])
	{
		for (auto read : paths[readi].reads)
		{
			if (closeReads.count(read.readNameIndex) == 1) continue;
			std::cerr << "Read close to assertion: " << graph.readNames[read.readNameIndex] << std::endl;
			closeReads.insert(read.readNameIndex);
		}
	}
	for (auto edge : graph.edges[std::make_pair(node, true)])
	{
		std::cerr << "around edge node " << edge.first << std::endl;
		for (auto readi : graph.readsCrossingNode[edge.first])
		{
			for (auto read : paths[readi].reads)
			{
				if (closeReads.count(read.readNameIndex) == 1 || semiCloseReads.count(read.readNameIndex) == 1) continue;
				std::cerr << "Read semi-close to assertion: " << graph.readNames[read.readNameIndex] << std::endl;
				semiCloseReads.insert(read.readNameIndex);
			}
		}
	}
	for (auto edge : graph.edges[std::make_pair(node, false)])
	{
		std::cerr << "around edge node " << edge.first << std::endl;
		for (auto readi : graph.readsCrossingNode[edge.first])
		{
			for (auto read : paths[readi].reads)
			{
				if (closeReads.count(read.readNameIndex) == 1 || semiCloseReads.count(read.readNameIndex) == 1) continue;
				std::cerr << "Read semi-close to assertion: " << graph.readNames[read.readNameIndex] << std::endl;
				semiCloseReads.insert(read.readNameIndex);
			}
		}
	}
}

void compact(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& paths, std::vector<size_t>& queueNodes)
{
	{
		std::vector<phmap::flat_hash_set<uint32_t>> tmp;
		std::swap(resolvableGraph.readsCrossingNode, tmp);
	}
	RankBitvector kept { resolvableGraph.unitigs.size() };
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		kept.set(i, !resolvableGraph.unitigRemoved[i]);
	}
	kept.buildRanks();
	size_t newSize = kept.getRank(kept.size()-1) + (kept.get(kept.size()-1) ? 1 : 0);
	for (size_t i = resolvableGraph.everTippable.size()-1; i < resolvableGraph.everTippable.size(); i--)
	{
		if (!kept.get(resolvableGraph.everTippable[i]))
		{
			std::swap(resolvableGraph.everTippable[i], resolvableGraph.everTippable.back());
			resolvableGraph.everTippable.pop_back();
			continue;
		}
		resolvableGraph.everTippable[i] = kept.getRank(resolvableGraph.everTippable[i]);
	}
	if (resolvableGraph.lastTippableChecked < resolvableGraph.unitigs.size())
	{
		resolvableGraph.lastTippableChecked = kept.getRank(resolvableGraph.lastTippableChecked);
	}
	else
	{
		resolvableGraph.lastTippableChecked = newSize;
	}
	{
		std::vector<size_t> newPrecalcedUnitigLengths;
		newPrecalcedUnitigLengths.resize(newSize, 0);
		for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			newPrecalcedUnitigLengths[newIndex] = resolvableGraph.precalcedUnitigLengths[i];
		}
		std::swap(newPrecalcedUnitigLengths, resolvableGraph.precalcedUnitigLengths);
	}
	resolvableGraph.unitigRemoved.resize(newSize);
	for (size_t i = 0; i < resolvableGraph.unitigRemoved.size(); i++)
	{
		resolvableGraph.unitigRemoved[i] = false;
	}
	{
		std::vector<std::vector<std::pair<size_t, bool>>> newUnitigs;
		std::vector<size_t> newUnitigLeftClipBp;
		std::vector<size_t> newUnitigRightClipBp;
		newUnitigs.resize(newSize);
		newUnitigRightClipBp.resize(newSize);
		newUnitigLeftClipBp.resize(newSize);
		for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			assert(newIndex < newUnitigs.size());
			assert(newIndex < newUnitigLeftClipBp.size());
			assert(newIndex < newUnitigRightClipBp.size());
			std::swap(newUnitigs[newIndex], resolvableGraph.unitigs[i]);
			newUnitigLeftClipBp[newIndex] = resolvableGraph.unitigLeftClipBp[i];
			newUnitigRightClipBp[newIndex] = resolvableGraph.unitigRightClipBp[i];
		}
		std::swap(resolvableGraph.unitigs, newUnitigs);
		std::swap(resolvableGraph.unitigRightClipBp, newUnitigRightClipBp);
		std::swap(resolvableGraph.unitigLeftClipBp, newUnitigLeftClipBp);
	}
	{
		VectorWithDirection<phmap::flat_hash_set<std::pair<size_t, bool>>> newEdges;
		newEdges.resize(newSize);
		for (size_t i = 0; i < resolvableGraph.edges.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			assert(newIndex < newEdges.size());
			for (auto target : resolvableGraph.edges[std::make_pair(i, true)])
			{
				assert(kept.get(target.first));
				assert(resolvableGraph.edges[reverse(target)].count(reverse(std::make_pair(i, true))) == 1);
				newEdges[std::make_pair(newIndex, true)].emplace(kept.getRank(target.first), target.second);
			}
			for (auto target : resolvableGraph.edges[std::make_pair(i, false)])
			{
				assert(kept.get(target.first));
				assert(resolvableGraph.edges[reverse(target)].count(reverse(std::make_pair(i, false))) == 1);
				newEdges[std::make_pair(newIndex, false)].emplace(kept.getRank(target.first), target.second);
			}
		}
		std::swap(resolvableGraph.edges, newEdges);
	}
	{
		phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newOverlaps;
		for (auto pair : resolvableGraph.overlaps)
		{
			auto from = pair.first.first;
			auto to = pair.first.second;
			if (!kept.get(from.first) || !kept.get(to.first)) continue;
			auto overlap = pair.second;
			from.first = kept.getRank(from.first);
			to.first = kept.getRank(to.first);
			assert(from.first < resolvableGraph.unitigs.size());
			assert(to.first < resolvableGraph.unitigs.size());
			newOverlaps[canon(from, to)] = overlap;
		}
		std::swap(resolvableGraph.overlaps, newOverlaps);
	}
	for (size_t i = paths.size()-1; i < paths.size(); i--)
	{
		if (paths[i].path.size() == 0 || paths[i].reads.size() == 0)
		{
			std::swap(paths[i], paths.back());
			paths.pop_back();
			continue;
		}
	}
	std::sort(paths.begin(), paths.end(), [](const PathGroup& left, const PathGroup& right) { return left.path < right.path; });
	for (size_t i = paths.size()-1; i > 0; i--)
	{
		if (paths[i].path != paths[i-1].path) continue;
		paths[i-1].reads.insert(paths[i-1].reads.end(), paths[i].reads.begin(), paths[i].reads.end());
		std::swap(paths[i], paths.back());
		paths.pop_back();
	}
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].path.size(); j++)
		{
			assert(kept.get(paths[i].path[j].first));
			paths[i].path[j].first = kept.getRank(paths[i].path[j].first);
			assert(paths[i].path[j].first < resolvableGraph.unitigs.size());
		}
	}
	for (size_t i = queueNodes.size()-1; i < queueNodes.size(); i--)
	{
		if (!kept.get(queueNodes[i]))
		{
			std::swap(queueNodes[i], queueNodes.back());
			queueNodes.pop_back();
			continue;
		}
		queueNodes[i] = kept.getRank(queueNodes[i]);
	}
	resolvableGraph.readsCrossingNode.resize(newSize);
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (auto pair : paths[i].path)
		{
			resolvableGraph.readsCrossingNode[pair.first].insert(i);
		}
	}
	for (size_t i = 0; i < resolvableGraph.unitigRemoved.size(); i++)
	{
		assert(!resolvableGraph.unitigRemoved[i]);
	}
}

class UnitigLengthComparer
{
public:
	UnitigLengthComparer(const ResolvableUnitigGraph& resolvableGraph) :
	resolvableGraph(resolvableGraph)
	{
	}
	bool operator()(const size_t left, const size_t right) const
	{
		return resolvableGraph.unitigLength(left) > resolvableGraph.unitigLength(right);
	}
private:
	const ResolvableUnitigGraph& resolvableGraph;
};

ResolvableUnitigGraph getUnitigs(const UnitigGraph& initial, size_t minCoverage, const HashList& hashlist, const size_t kmerSize, const bool keepGaps)
{
	ResolvableUnitigGraph result { hashlist, kmerSize };
	result.unitigs.resize(initial.unitigs.size());
	result.unitigLeftClipBp.resize(initial.unitigs.size(), 0);
	result.unitigRightClipBp.resize(initial.unitigs.size(), 0);
	result.unitigRemoved.resize(initial.unitigs.size(), false);
	result.edges.resize(initial.unitigs.size());
	result.readsCrossingNode.resize(initial.unitigs.size());
	std::vector<std::pair<size_t, bool>> checkKeepTips;
	for (size_t i = 0; i < initial.unitigs.size(); i++)
	{
		result.unitigs[i].insert(result.unitigs[i].end(), initial.unitigs[i].begin(), initial.unitigs[i].end());
		std::pair<size_t, bool> fw { i, true };
		bool keepCheck = false;
		for (auto pair : initial.edgeCov[fw])
		{
			if (pair.second < minCoverage)
			{
				keepCheck = true;
				continue;
			}
			auto canonpair = canon(fw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[fw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(fw));
		}
		if (keepCheck) checkKeepTips.push_back(fw);
		std::pair<size_t, bool> bw { i, false };
		keepCheck = false;
		for (auto pair : initial.edgeCov[bw])
		{
			if (pair.second < minCoverage)
			{
				keepCheck = true;
				continue;
			}
			auto canonpair = canon(bw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[bw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(bw));
		}
		if (keepCheck) checkKeepTips.push_back(bw);
	}
	if (keepGaps)
	{
		std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> addEdges;
		for (auto tip : checkKeepTips)
		{
			if (result.edges[tip].size() != 0) continue;
			for (auto pair : initial.edgeCov[tip])
			{
				if (result.edges[reverse(pair.first)].size() != 0) continue;
				auto canonpair = canon(tip, pair.first);
				addEdges.emplace_back(canonpair);
			}
		}
		for (auto pair : addEdges)
		{
			result.overlaps[pair] = 0;
			result.edges[pair.first].emplace(pair.second);
			result.edges[reverse(pair.second)].emplace(reverse(pair.first));
		}
	}
	double allkmerSum = 0;
	double allkmerDivisor = 0;
	double longkmerSum = 0;
	double longkmerDivisor = 0;
	for (size_t i = 0; i < result.unitigs.size(); i++)
	{
		double sumHere = 0;
		for (size_t j = 0; j < result.unitigs[i].size(); j++)
		{
			sumHere += hashlist.coverage.get(result.unitigs[i][j].first);
		}
		allkmerSum += sumHere;
		allkmerDivisor += result.unitigs[i].size();
		if (result.unitigLength(i) > 100000)
		{
			longkmerSum += sumHere;
			longkmerDivisor += result.unitigs[i].size();
		}
	}
	if (longkmerDivisor > 0)
	{
		result.averageCoverage = longkmerSum / longkmerDivisor;
	}
	else
	{
		assert(longkmerDivisor > 0 || longkmerSum == 0);
		result.averageCoverage = 0;
		if (longkmerDivisor > 0)
		{
			result.averageCoverage = longkmerSum / longkmerDivisor;
		}
	}
	return result;
}

size_t getNumberOfHashes(const ResolvableUnitigGraph& resolvableGraph, size_t leftClip, size_t rightClip, const std::vector<std::pair<size_t, bool>>& path)
{
	if (path.size() == 0) return 0;
	size_t result = 0;
	assert(path.size() < 2 || leftClip < resolvableGraph.unitigs[path[0].first].size());
	assert(path.size() < 2 || rightClip < resolvableGraph.unitigs[path.back().first].size());
	assert(path.size() != 1 || leftClip + rightClip < resolvableGraph.unitigs[path[0].first].size());
	for (size_t i = 0; i < path.size(); i++)
	{
		result += resolvableGraph.unitigs[path[i].first].size();
		if (i > 0)
		{
			size_t overlap = resolvableGraph.overlaps.at(canon(path[i-1], path[i]));
			assert(resolvableGraph.unitigs[path[i].first].size() >= overlap);
			result -= overlap;
		}
	}
	assert(result > leftClip + rightClip);
	result -= leftClip;
	result -= rightClip;
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPath>> resolvableToUnitigs(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths)
{
	{
		std::vector<std::vector<size_t>> coverages;
		coverages.resize(resolvableGraph.unitigs.size());
		for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
		{
			coverages[i].resize(resolvableGraph.unitigs[i].size(), 0);
		}
		for (const auto& path : readPaths)
		{
			if (path.path.size() == 0) continue;
			for (size_t i = 0; i < path.path.size(); i++)
			{
				assert(path.path[i].first < resolvableGraph.unitigs.size());
				assert(!resolvableGraph.unitigRemoved[path.path[i].first]);
			}
			for (size_t i = 1; i < path.path.size(); i++)
			{
				assert(resolvableGraph.edges[path.path[i-1]].count(path.path[i]) == 1);
				assert(resolvableGraph.edges[reverse(path.path[i])].count(reverse(path.path[i-1])) == 1);
			}
			size_t pathHashCount = getNumberOfHashes(resolvableGraph, 0, 0, path.path);
			if (path.path.size() == 1)
			{
				for (const auto& read : path.reads)
				{
					assert(read.leftClip + read.rightClip + read.readPoses.size() == pathHashCount);
					assert(read.leftClip + read.rightClip < resolvableGraph.unitigs[path.path[0].first].size());
					for (size_t i = read.leftClip; i < resolvableGraph.unitigs[path.path[0].first].size() - read.rightClip; i++)
					{
						size_t index = i;
						if (!path.path[0].second) index = resolvableGraph.unitigs[path.path[0].first].size() - i - 1;
						coverages[path.path[0].first][index] += 1;
					}
				}
				continue;
			}
			assert(path.path.size() >= 2);
			for (size_t i = 1; i < path.path.size()-1; i++)
			{
				for (size_t j = 0; j < resolvableGraph.unitigs[path.path[i].first].size(); j++)
				{
					coverages[path.path[i].first][j] += path.reads.size();
				}
			}
			for (const auto& read : path.reads)
			{
				assert(read.leftClip + read.rightClip + read.readPoses.size() == pathHashCount);
				assert(read.leftClip < resolvableGraph.unitigs[path.path[0].first].size());
				for (size_t i = read.leftClip; i < resolvableGraph.unitigs[path.path[0].first].size(); i++)
				{
					size_t index = i;
					if (!path.path[0].second) index = resolvableGraph.unitigs[path.path[0].first].size() - i - 1;
					coverages[path.path[0].first][index] += 1;
				}
				assert(read.rightClip < resolvableGraph.unitigs[path.path.back().first].size());
				for (size_t i = 0; i < resolvableGraph.unitigs[path.path.back().first].size() - read.rightClip; i++)
				{
					size_t index = i;
					if (!path.path.back().second) index = resolvableGraph.unitigs[path.path.back().first].size() - i - 1;
					coverages[path.path.back().first][index] += 1;
				}
			}
		}
		for (size_t i = 0; i < coverages.size(); i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			for (size_t j = 0; j < coverages[i].size(); j++)
			{
				assertPrintReads(coverages[i][j] > 0, resolvableGraph, readPaths, i);
			}
		}
	}
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		assert(path.reads.size() > 0);
		size_t pathHashCount = getNumberOfHashes(resolvableGraph, 0, 0, path.path);
		for (const auto& read : path.reads)
		{
			assert(read.leftClip + read.rightClip + read.readPoses.size() == pathHashCount);
		}
	}
	UnitigGraph result;
	RankBitvector newIndex { resolvableGraph.unitigs.size() };
	assert(resolvableGraph.unitigs.size() == resolvableGraph.unitigRemoved.size());
	assert(resolvableGraph.unitigs.size() == resolvableGraph.edges.size());
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		newIndex.set(i, !resolvableGraph.unitigRemoved[i]);
	}
	newIndex.buildRanks();
	const size_t newSize = newIndex.getRank(newIndex.size()-1) + (newIndex.get(newIndex.size()-1) ? 1 : 0);
	result.unitigs.resize(newSize);
	result.leftClip.resize(newSize);
	result.rightClip.resize(newSize);
	result.unitigCoverage.resize(newSize);
	result.edges.resize(newSize);
	result.edgeCov.resize(newSize);
	result.edgeOvlp.resize(newSize);
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (!newIndex.get(i)) continue;
		const size_t unitig = newIndex.getRank(i);
		assert(unitig < newSize);
		result.unitigs[unitig].insert(result.unitigs[unitig].end(), resolvableGraph.unitigs[i].begin(), resolvableGraph.unitigs[i].end());
		for (size_t j = 0; j < resolvableGraph.unitigs[i].size(); j++)
		{
			result.unitigCoverage[unitig].push_back(0);
		}
		result.leftClip[unitig] = resolvableGraph.unitigLeftClipBp[i];
		result.rightClip[unitig] = resolvableGraph.unitigRightClipBp[i];
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> newfw { newIndex.getRank(i), true };
		assert(newfw.first < newSize);
		for (auto edge : resolvableGraph.edges[fw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			assert(newEdge.first < newSize);
			result.edges[newfw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newfw));
			result.edgeCoverage(newfw, newEdge) = 0;
			result.edgeOverlap(newfw, newEdge) = resolvableGraph.overlaps.at(canon(fw, edge));
		}
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newbw { newIndex.getRank(i), false };
		assert(newbw.first < newSize);
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			assert(newEdge.first < newSize);
			result.edges[newbw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newbw));
			result.edgeCoverage(newbw, newEdge) = 0;
			result.edgeOverlap(newbw, newEdge) = resolvableGraph.overlaps.at(canon(bw, edge));
		}
	}
	std::vector<ReadPath> resultReads;
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		assert(path.reads.size() > 0);
		std::vector<std::pair<size_t, bool>> fixPath = path.path;
		for (size_t i = 0; i < fixPath.size(); i++)
		{
			assert(newIndex.get(fixPath[i].first));
			fixPath[i].first = newIndex.getRank(fixPath[i].first);
			assert(fixPath[i].first < newSize);
		}
		assert(fixPath.size() > 0);
		if (fixPath.size() == 1)
		{
			for (const auto& read : path.reads)
			{
				for (size_t i = read.leftClip; i < result.unitigCoverage[fixPath[0].first].size() - read.rightClip; i++)
				{
					size_t index = i;
					if (!fixPath[0].second) index = result.unitigCoverage[fixPath[0].first].size() - 1 - i;
					result.unitigCoverage[fixPath[0].first][index] += 1;
				}
			}
		}
		else
		{
			assert(fixPath.size() >= 2);
			for (size_t i = 1; i < fixPath.size()-1; i++)
			{
				for (size_t j = 0; j < result.unitigCoverage[fixPath[i].first].size(); j++)
				{
					result.unitigCoverage[fixPath[i].first][j] += path.reads.size();
				}
			}
			for (const auto& read : path.reads)
			{
				for (size_t i = read.leftClip; i < result.unitigCoverage[fixPath[0].first].size(); i++)
				{
					size_t index = i;
					if (!fixPath[0].second) index = result.unitigCoverage[fixPath[0].first].size()-1-i;
					result.unitigCoverage[fixPath[0].first][index] += 1;
				}
			}
			for (const auto& read : path.reads)
			{
				for (size_t i = 0; i < result.unitigCoverage[fixPath.back().first].size() - read.rightClip; i++)
				{
					size_t index = i;
					if (!fixPath.back().second) index = result.unitigCoverage[fixPath.back().first].size()-1-i;
					result.unitigCoverage[fixPath.back().first][index] += 1;
				}
			}
		}
		for (size_t j = 1; j < fixPath.size(); j++)
		{
			result.edgeCoverage(fixPath[j-1], fixPath[j]) += path.reads.size();
		}
		for (const auto& read : path.reads)
		{
			resultReads.emplace_back();
			resultReads.back().path = fixPath;
			resultReads.back().readName = resolvableGraph.readNames[read.readNameIndex];
			resultReads.back().readPoses = read.readPoses;
			resultReads.back().expandedReadPosStart = read.expandedReadPosStart;
			resultReads.back().expandedReadPosEnd = read.expandedReadPosEnd;
			resultReads.back().leftClip = read.leftClip;
			resultReads.back().rightClip = read.rightClip;
			resultReads.back().readLength = read.readLength;
			resultReads.back().readLengthHPC = read.readLengthHPC;
		}
	}
	for (size_t i = 0; i < result.unitigCoverage.size(); i++)
	{
		for (size_t j = 0; j < result.unitigCoverage[i].size(); j++)
		{
			assert(result.unitigCoverage[i][j] > 0);
		}
	}
	return std::make_pair(result, resultReads);
}

std::vector<std::pair<size_t, bool>> extend(const ResolvableUnitigGraph& resolvableGraph, const std::pair<size_t, bool> start)
{
	std::pair<size_t, bool> pos = start;
	std::vector<std::pair<size_t, bool>> result;
	result.emplace_back(pos);
	assert(!resolvableGraph.unitigRemoved[start.first]);
	while (true)
	{
		assert(!resolvableGraph.unitigRemoved[pos.first]);
		if (resolvableGraph.edges[pos].size() != 1) break;
		auto newPos = *resolvableGraph.edges[pos].begin();
		auto revNewPos = reverse(newPos);
		if (resolvableGraph.edges[revNewPos].size() != 1) break;
		if (newPos == start)
		{
			// include the start to note that it's circular
			result.emplace_back(newPos);
			break;
		}
		if (newPos.first == pos.first) break;
		pos = newPos;
		result.emplace_back(pos);
	}
	return result;
}

void erasePath(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t i)
{
	assert(readPaths[i].path.size() > 0);
	for (auto node : readPaths[i].path)
	{
		if (resolvableGraph.readsCrossingNode[node.first].count(i) == 1) resolvableGraph.readsCrossingNode[node.first].erase(i);
	}
	readPaths[i].path.clear();
	readPaths[i].reads.clear();
}

void addPath(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, PathGroup&& newPath)
{
	assert(newPath.path.size() > 0);
	assert(newPath.reads.size() > 0);
	for (const auto& read : newPath.reads)
	{
		// assert(read.leftClip + read.rightClip + read.readPoses.size() == getNumberOfHashes(resolvableGraph, 0, 0, newPath.path));
		if (newPath.path.size() == 1)
		{
			assert(resolvableGraph.unitigs[newPath.path[0].first].size() > read.leftClip + read.rightClip);
		}
		else
		{
			assert(resolvableGraph.unitigs[newPath.path[0].first].size() > read.leftClip);
			assert(resolvableGraph.unitigs[newPath.path.back().first].size() > read.rightClip);
		}
	}
	for (size_t i = 0; i < newPath.path.size(); i++)
	{
		assert(!resolvableGraph.unitigRemoved[newPath.path[i].first]);
		if (i > 0)
		{
			assert(resolvableGraph.edges[newPath.path[i-1]].count(newPath.path[i]) == 1);
			assert(resolvableGraph.edges[reverse(newPath.path[i])].count(reverse(newPath.path[i-1])) == 1);
		}
	}
	size_t newIndex = readPaths.size();
	for (auto pos : newPath.path)
	{
		resolvableGraph.readsCrossingNode[pos.first].insert(newIndex);
	}
	readPaths.emplace_back(std::move(newPath));
}

void replacePathNodes(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::vector<std::pair<std::vector<std::pair<size_t, bool>>, size_t>>& newUnitigs)
{
	if (newUnitigs.size() == 0) return;
	phmap::flat_hash_map<size_t, std::pair<size_t, size_t>> unitigIndex;
	std::vector<std::vector<size_t>> leftClip;
	std::vector<std::vector<size_t>> rightClip;
	leftClip.resize(newUnitigs.size());
	rightClip.resize(newUnitigs.size());
	for (size_t i = 0; i < newUnitigs.size(); i++)
	{
		if (newUnitigs[i].first.size() == 1) continue;
		assert(newUnitigs[i].first.size() >= 2);
		for (size_t j = 0; j < newUnitigs[i].first.size(); j++)
		{
			assert(unitigIndex.count(newUnitigs[i].first[j].first) == 0);
			unitigIndex[newUnitigs[i].first[j].first] = std::make_pair(i, j);
		}
		leftClip[i].resize(newUnitigs[i].first.size());
		rightClip[i].resize(newUnitigs[i].first.size());
		assert(resolvableGraph.unitigs[newUnitigs[i].first[0].first].size() <= resolvableGraph.unitigs[newUnitigs[i].second].size());
		assert(resolvableGraph.unitigLength(newUnitigs[i].first[0].first) < resolvableGraph.unitigLength(newUnitigs[i].second));
		size_t leftClipSum = 0;
		size_t rightClipSum = resolvableGraph.unitigs[newUnitigs[i].second].size() - resolvableGraph.unitigs[newUnitigs[i].first[0].first].size();
		leftClip[i][0] = leftClipSum;
		rightClip[i][0] = rightClipSum;
		assert(leftClipSum + rightClipSum + resolvableGraph.unitigs[newUnitigs[i].first[0].first].size() == resolvableGraph.unitigs[newUnitigs[i].second].size());
		for (size_t j = 1; j < newUnitigs[i].first.size(); j++)
		{
			size_t overlap = resolvableGraph.overlaps[canon(newUnitigs[i].first[j-1], newUnitigs[i].first[j])];
			assert(resolvableGraph.unitigs[newUnitigs[i].first[j-1].first].size() >= overlap);
			assert(resolvableGraph.unitigs[newUnitigs[i].first[j].first].size() >= overlap);
			assert(resolvableGraph.unitigs[newUnitigs[i].first[j].first].size() - overlap <= rightClipSum);
			leftClipSum = leftClipSum + resolvableGraph.unitigs[newUnitigs[i].first[j-1].first].size() - overlap;
			rightClipSum = rightClipSum - resolvableGraph.unitigs[newUnitigs[i].first[j].first].size() + overlap;
			leftClip[i][j] = leftClipSum;
			rightClip[i][j] = rightClipSum;
			assert(leftClipSum + rightClipSum + resolvableGraph.unitigs[newUnitigs[i].first[j].first].size() == resolvableGraph.unitigs[newUnitigs[i].second].size());
		}
		assert(rightClipSum == 0);
	}
	if (unitigIndex.size() == 0) return;
	phmap::flat_hash_set<size_t> relevantReads;
	for (size_t i = 0; i < newUnitigs.size(); i++)
	{
		if (newUnitigs[i].first.size() >= 2)
		{
			for (auto pos : newUnitigs[i].first)
			{
				relevantReads.insert(resolvableGraph.readsCrossingNode[pos.first].begin(), resolvableGraph.readsCrossingNode[pos.first].end());
			}
		}
	}
	for (const size_t i : relevantReads)
	{
		PathGroup newPath;
		newPath.path.reserve(readPaths[i].path.size());
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			auto inUnitig = unitigIndex.find(readPaths[i].path[j].first);
			if (inUnitig == unitigIndex.end())
			{
				newPath.path.emplace_back(readPaths[i].path[j]);
				continue;
			}
			size_t unitigindex = inUnitig->second.first;
			size_t posindex = inUnitig->second.second;
			if (j == 0)
			{
				bool fw = newUnitigs[unitigindex].first[posindex].second;
				if (!readPaths[i].path[j].second) fw = !fw;
				newPath.path.emplace_back(newUnitigs[unitigindex].second, fw);
			}
			else if (readPaths[i].path[j] == newUnitigs[unitigindex].first[0])
			{
				newPath.path.emplace_back(newUnitigs[unitigindex].second, true);
			}
			else if (readPaths[i].path[j] == reverse(newUnitigs[unitigindex].first.back()))
			{
				newPath.path.emplace_back(newUnitigs[unitigindex].second, false);
			}
		}
		size_t extraLeftClip = 0;
		size_t extraRightClip = 0;
		auto firstInUnitig = unitigIndex.find(readPaths[i].path[0].first);
		if (firstInUnitig != unitigIndex.end())
		{
			size_t unitigindex = firstInUnitig->second.first;
			size_t posindex = firstInUnitig->second.second;
			bool fw = newUnitigs[unitigindex].first[posindex].second;
			if (!readPaths[i].path[0].second) fw = !fw;
			if (fw)
			{
				extraLeftClip += leftClip[unitigindex][posindex];
			}
			else
			{
				extraLeftClip += rightClip[unitigindex][posindex];
			}
		}
		auto lastInUnitig = unitigIndex.find(readPaths[i].path.back().first);
		if (lastInUnitig != unitigIndex.end())
		{
			size_t unitigindex = lastInUnitig->second.first;
			size_t posindex = lastInUnitig->second.second;
			bool fw = newUnitigs[unitigindex].first[posindex].second;
			if (!readPaths[i].path.back().second) fw = !fw;
			if (fw)
			{
				extraRightClip += rightClip[unitigindex][posindex];
			}
			else
			{
				extraRightClip += leftClip[unitigindex][posindex];
			}
		}
		std::swap(newPath.reads, readPaths[i].reads);
		if (extraLeftClip > 0 || extraRightClip > 0)
		{
			for (auto& read : newPath.reads)
			{
				read.leftClip += extraLeftClip;
				read.rightClip += extraRightClip;
			}
		}
		addPath(resolvableGraph, readPaths, std::move(newPath));
		erasePath(resolvableGraph, readPaths, i);
	}
}

// todo maybe fix? or does it matter?
void cutRemovedEdgesFromPaths(ResolvableUnitigGraph& graph, std::vector<ReadPath>& readPaths)
{
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		if (readPaths[i].path.size() == 0) continue;
		bool remove = false;
		size_t lastStart = 0;
		for (size_t j = 1; j < readPaths[i].path.size(); j++)
		{
			if (graph.edges[readPaths[i].path[j-1]].count(readPaths[i].path[j]) == 1) continue;
			remove = true;
		}
		if (!remove) continue;
		std::vector<size_t> pathStartPoses;
		std::vector<size_t> pathEndPoses;
		size_t pos = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			pathStartPoses.push_back(pos);
			pos += graph.unitigs[readPaths[i].path[j].first].size();
			pathEndPoses.push_back(pos);
		}
		for (size_t j = 1; j < readPaths[i].path.size(); j++)
		{
			if (graph.edges[readPaths[i].path[j-1]].count(readPaths[i].path[j]) == 1) continue;
			ReadPath newPath;
			newPath.readName = readPaths[i].readName;
			newPath.readLength = readPaths[i].readLength;
			newPath.readLengthHPC = readPaths[i].readLengthHPC;
			newPath.path.insert(newPath.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.begin() + j);
			newPath.leftClip = 0;
			if (lastStart == 0) newPath.leftClip = readPaths[i].leftClip;
			newPath.rightClip = 0;
			size_t wantedStart = 0;
			if (pathStartPoses[lastStart] > readPaths[i].leftClip) wantedStart = pathStartPoses[lastStart] - readPaths[i].leftClip;
			assert(readPaths[i].leftClip < pathEndPoses[j-1]);
			size_t wantedEnd = pathEndPoses[j-1] - readPaths[i].leftClip;
			assert(wantedEnd < readPaths[i].readPoses.size());
			newPath.readPoses.insert(newPath.readPoses.end(), readPaths[i].readPoses.begin() + wantedStart, readPaths[i].readPoses.begin() + wantedEnd);
			readPaths.emplace_back(std::move(newPath));
			lastStart = j;
		}
		assert(lastStart != 0);
		ReadPath newPath;
		newPath.readName = readPaths[i].readName;
		newPath.readLength = readPaths[i].readLength;
		newPath.readLengthHPC = readPaths[i].readLengthHPC;
		newPath.path.insert(newPath.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.end());
		newPath.leftClip = 0;
		if (lastStart == 0) newPath.leftClip = readPaths[i].leftClip;
		newPath.rightClip = readPaths[i].rightClip;
		size_t wantedStart = 0;
		if (pathStartPoses[lastStart] > readPaths[i].leftClip) wantedStart = pathStartPoses[lastStart] - readPaths[i].leftClip;
		assert(readPaths[i].leftClip < pathEndPoses.back());
		assert(pathEndPoses.back() > readPaths[i].leftClip + readPaths[i].rightClip);
		size_t wantedEnd = pathEndPoses.back() - readPaths[i].leftClip - readPaths[i].rightClip;
		assert(wantedEnd == readPaths[i].readPoses.size());
		newPath.readPoses.insert(newPath.readPoses.end(), readPaths[i].readPoses.begin() + wantedStart, readPaths[i].readPoses.begin() + wantedEnd);
		readPaths.emplace_back(std::move(newPath));
		readPaths[i].path.clear();
	}
}

std::vector<std::pair<size_t, bool>> getUnitigPath(const ResolvableUnitigGraph& resolvableGraph, const size_t unitig)
{
	std::vector<std::pair<size_t, bool>> fwExtension = extend(resolvableGraph, std::make_pair(unitig, true));
	std::vector<std::pair<size_t, bool>> bwExtension = extend(resolvableGraph, std::make_pair(unitig, false));
	assert(bwExtension.size() >= 1);
	assert(fwExtension.size() >= 1);
	std::vector<std::pair<size_t, bool>> newUnitig = revCompPath(bwExtension);
	if (newUnitig.size() >= 2 && newUnitig[0] == newUnitig.back())
	{
		// circular
		assert(newUnitig.size() == fwExtension.size());
		for (size_t i = 0; i < fwExtension.size(); i++)
		{
			assert(newUnitig[i] == fwExtension[i]);
		}
		newUnitig.pop_back();
	}
	else
	{
		assert(newUnitig.back() == std::make_pair(unitig, true));
		newUnitig.insert(newUnitig.end(), fwExtension.begin()+1, fwExtension.end());
	}
	assert(newUnitig.size() >= 1);
	return newUnitig;
}

size_t replaceUnitigPath(ResolvableUnitigGraph& resolvableGraph, const std::vector<std::pair<size_t, bool>>& newUnitig)
{
	size_t newIndex = resolvableGraph.unitigs.size();
	size_t leftClip, rightClip;
	if (newUnitig[0].second)
	{
		leftClip = resolvableGraph.unitigLeftClipBp[newUnitig[0].first];
	}
	else
	{
		leftClip = resolvableGraph.unitigRightClipBp[newUnitig[0].first];
	}
	if (newUnitig.back().second)
	{
		rightClip = resolvableGraph.unitigRightClipBp[newUnitig.back().first];
	}
	else
	{
		rightClip = resolvableGraph.unitigLeftClipBp[newUnitig.back().first];
	}
	resolvableGraph.unitigLeftClipBp.push_back(leftClip);
	resolvableGraph.unitigRightClipBp.push_back(rightClip);
	resolvableGraph.unitigs.emplace_back();
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	for (size_t i = 0; i < newUnitig.size(); i++)
	{
		assert(i == 0 || resolvableGraph.edges[reverse(newUnitig[i])].size() == 1);
		assert(i == newUnitig.size()-1 || resolvableGraph.edges[newUnitig[i]].size() == 1);
		std::vector<std::pair<size_t, bool>> add = resolvableGraph.unitigs[newUnitig[i].first];
		if (!newUnitig[i].second) add = revCompPath(add);
		size_t overlap = 0;
		if (i > 0) overlap = resolvableGraph.overlaps.at(canon(newUnitig[i-1], newUnitig[i]));
		resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin()+overlap, add.end());
	}
	phmap::flat_hash_set<size_t> nodesInUnitig;
	for (auto pos : newUnitig)
	{
		assert(nodesInUnitig.count(pos.first) == 0);
		nodesInUnitig.insert(pos.first);
	}
	std::pair<size_t, bool> bw { newIndex, false };
	std::pair<size_t, bool> fw { newIndex, true };
	bool removeSelfWeirdoLoop = false;
	for (auto edge : resolvableGraph.edges[reverse(newUnitig[0])])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(bw, bw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			assert(edge.first != newUnitig[0].first);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
		else if (edge == newUnitig[0])
		{
			resolvableGraph.overlaps[canon(bw, fw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(fw);
			resolvableGraph.edges[reverse(fw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			removeSelfWeirdoLoop = true;
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(bw, edge)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			assert(edge.first != newUnitig[0].first);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
	}
	if (removeSelfWeirdoLoop)
	{
		assert(resolvableGraph.edges[reverse(newUnitig[0])].count(newUnitig[0]) == 1);
		resolvableGraph.edges[reverse(newUnitig[0])].erase(newUnitig[0]);
	}
	removeSelfWeirdoLoop = false;
	for (auto edge : resolvableGraph.edges[newUnitig.back()])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(fw, bw)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			assert(edge.first != newUnitig[0].first);
			removeSelfWeirdoLoop = true;
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(fw, edge)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			assert(edge.first != newUnitig[0].first);
			resolvableGraph.edges[reverse(edge)].erase(reverse(newUnitig.back()));
		}
	}
	if (removeSelfWeirdoLoop)
	{
		assert(resolvableGraph.edges[newUnitig.back()].count(reverse(newUnitig.back())) == 1);
		resolvableGraph.edges[newUnitig.back()].erase(reverse(newUnitig.back()));
	}
	for (size_t i = 0; i < newUnitig.size(); i++)
	{
		resolvableGraph.unitigRemoved[newUnitig[i].first] = true;
		resolvableGraph.edges[newUnitig[i]].clear();
		resolvableGraph.edges[reverse(newUnitig[i])].clear();
	}
	return newIndex;
}

std::vector<std::pair<std::vector<std::pair<size_t, bool>>, size_t>> unitigifySet(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::vector<size_t>& unitigifiable)
{
	phmap::flat_hash_set<size_t> inUnitig;
	std::vector<std::pair<std::vector<std::pair<size_t, bool>>, size_t>> result;
	for (auto node : unitigifiable)
	{
		if (resolvableGraph.unitigRemoved[node]) continue;
		if (inUnitig.count(node) == 1) continue;
		auto newUnitig = getUnitigPath(resolvableGraph, node);
		if (newUnitig.size() >= 2)
		{
			for (auto node : newUnitig)
			{
				assert(inUnitig.count(node.first) == 0);
				inUnitig.insert(node.first);
			}
		}
		result.emplace_back(std::move(newUnitig), 0);
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		if (result[i].first.size() == 1)
		{
			result[i].second = result[i].first[0].first;
			assert(result[i].first[0].second);
		}
		else
		{
			assert(result[i].first.size() >= 2);
			size_t newIndex = replaceUnitigPath(resolvableGraph, result[i].first);
			result[i].second = newIndex;
		}
	}
	replacePathNodes(resolvableGraph, readPaths, result);
	return result;
}

std::vector<std::pair<size_t, bool>> unitigifyOne(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t unitig)
{
	std::vector<std::pair<size_t, bool>> newUnitig = getUnitigPath(resolvableGraph, unitig);
	if (newUnitig.size() == 1) return newUnitig;
	size_t newIndex = replaceUnitigPath(resolvableGraph, newUnitig);
	std::vector<std::pair<std::vector<std::pair<size_t, bool>>, size_t>> unitigifyThis;
	unitigifyThis.emplace_back(std::move(newUnitig), newIndex);
	replacePathNodes(resolvableGraph, readPaths, unitigifyThis);
	return newUnitig;
}

void unitigifyAll(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths)
{
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		unitigifyOne(resolvableGraph, readPaths, i);
	}
}

void unresolveRecursively(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, phmap::flat_hash_set<size_t>& unresolvables, const size_t node)
{
	unresolvables.emplace(node);
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		if (resolvables.count(edge.first) == 0) continue;
		if (unresolvables.count(edge.first) == 1) continue;
		unresolveRecursively(resolvableGraph, resolvables, unresolvables, edge.first);
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		if (resolvables.count(edge.first) == 0) continue;
		if (unresolvables.count(edge.first) == 1) continue;
		unresolveRecursively(resolvableGraph, resolvables, unresolvables, edge.first);
	}
}

void createEdgeNode(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, const size_t kmerSize, phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes, const phmap::flat_hash_set<size_t>& resolvables, const phmap::flat_hash_set<size_t>& unresolvables, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	size_t newIndex = resolvableGraph.unitigs.size();
	newEdgeNodes[std::make_pair(from, to)] = newIndex;
	resolvableGraph.unitigs.emplace_back();
	std::vector<std::pair<size_t, bool>> add = resolvableGraph.unitigs[from.first];
	size_t leftClipBp = resolvableGraph.unitigLeftClipBp[from.first];
	size_t rightClipBp = resolvableGraph.unitigRightClipBp[from.first];
	if (!from.second)
	{
		add = revCompPath(add);
		std::swap(leftClipBp, rightClipBp);
	}
	resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin(), add.end());
	size_t overlapIncrement = 0;
	if (resolvables.count(to.first) == 1 && unresolvables.count(to.first) == 0)
	{
		size_t start = resolvableGraph.overlaps.at(canon(from, to));
		assert(start <= resolvableGraph.unitigs[to.first].size());
		add = resolvableGraph.unitigs[to.first];
		if (!to.second) add = revCompPath(add);
		resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin()+start, add.end());
		if (to.second)
		{
			rightClipBp = resolvableGraph.unitigRightClipBp[to.first];
		}
		else
		{
			rightClipBp = resolvableGraph.unitigLeftClipBp[to.first];
		}
	}
	else
	{
		size_t availableIncrease = resolvableGraph.unitigLength(to.first) - resolvableGraph.getBpOverlap(from, to) - 1;
		assert(availableIncrease >= 1);
		if (rightClipBp == 0)
		{
			size_t start = resolvableGraph.overlaps.at(canon(from, to));
			assert(start < resolvableGraph.unitigs[to.first].size());
			add = resolvableGraph.unitigs[to.first];
			if (!to.second) add = revCompPath(add);
			overlapIncrement = 1;
			for (size_t i = 0; i < start; i++)
			{
				assert(resolvableGraph.unitigs.back()[resolvableGraph.unitigs.back().size() - start+i] == add[i]);
			}
			assert(add.size() > start);
			rightClipBp = kmerSize - hashlist.getOverlap(resolvableGraph.unitigs.back().back(), add[start]) - 1;
			assert(rightClipBp < kmerSize);
			if (rightClipBp > availableIncrease - 1)
			{
				rightClipBp -= availableIncrease - 1;
			}
			else if (rightClipBp > 1)
			{
				rightClipBp = 1;
			}
			resolvableGraph.unitigs.back().emplace_back(add[start]);
		}
		else
		{
			assert(rightClipBp < kmerSize);
			if (rightClipBp > availableIncrease)
			{
				rightClipBp -= availableIncrease;
			}
			else
			{
				rightClipBp = 0;
			}
		}
	}
	resolvableGraph.unitigRightClipBp.push_back(rightClipBp);
	resolvableGraph.unitigLeftClipBp.push_back(leftClipBp);
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	if (resolvables.count(to.first) == 0 || unresolvables.count(to.first) == 1)
	{
		resolvableGraph.edges[std::make_pair(newIndex, true)].emplace(to);
		resolvableGraph.edges[reverse(to)].emplace(std::make_pair(newIndex, false));
		size_t overlap = resolvableGraph.overlaps.at(canon(from, to));
		assert(resolvableGraph.unitigs[newIndex].size() >= overlap + overlapIncrement);
		assert(resolvableGraph.unitigs[to.first].size() >= overlap + overlapIncrement);
		resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = overlap + overlapIncrement;
		// resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = resolvableGraph.unitigs[from.first].size();
		assert(resolvableGraph.getBpOverlap(std::make_pair(newIndex, true), to) < resolvableGraph.unitigLength(to.first));
	}
	assert(resolvableGraph.unitigLength(newIndex) > resolvableGraph.unitigLength(from.first));
	if (resolvableGraph.unitigs[newIndex].size() >= 2)
	{
		assert(resolvableGraph.unitigLeftClipBp[newIndex] < kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][0], resolvableGraph.unitigs[newIndex][1]));
		assert(resolvableGraph.unitigRightClipBp[newIndex] < kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-2], resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-1]));
	}
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getRawTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage)
{
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> tripletCoverage;
	for (const size_t i : resolvableGraph.readsCrossingNode[node])
	{
		if (readPaths[i].path.size() == 0) continue;
		if (readPaths[i].path.size() == 1) continue;
		if (readPaths[i].reads.size() == 0) continue;
		if (readPaths[i].path[0].first == node && resolvableGraph.edges[reverse(readPaths[i].path[0])].size() == 0)
		{
			size_t coverageHere = readPaths[i].reads.size();
			if (coverageHere > 0)
			{
				std::pair<size_t, bool> left { std::numeric_limits<size_t>::max(), true };
				std::pair<size_t, bool> right { std::numeric_limits<size_t>::max(), true };
				if (readPaths[i].path[0].second)
				{
					right = readPaths[i].path[1];
				}
				else
				{
					left = reverse(readPaths[i].path[1]);
				}
				tripletCoverage[std::make_pair(left, right)] += coverageHere;
			}
		}
		if (readPaths[i].path.back().first == node && resolvableGraph.edges[reverse(readPaths[i].path.back())].size() == 0)
		{
			size_t coverageHere = readPaths[i].reads.size();
			if (coverageHere > 0)
			{
				std::pair<size_t, bool> left { std::numeric_limits<size_t>::max(), true };
				std::pair<size_t, bool> right { std::numeric_limits<size_t>::max(), true };
				if (readPaths[i].path[0].second)
				{
					left = readPaths[i].path[readPaths[i].path.size()-2];
				}
				else
				{
					right = reverse(readPaths[i].path[readPaths[i].path.size()-2]);
				}
				tripletCoverage[std::make_pair(left, right)] += coverageHere;
			}
		}
		for (size_t j = 1; j < readPaths[i].path.size()-1; j++)
		{
			if (readPaths[i].path[j].first != node) continue;
			std::pair<size_t, bool> left;
			std::pair<size_t, bool> right;
			if (readPaths[i].path[j].second)
			{
				left = readPaths[i].path[j-1];
				right = readPaths[i].path[j+1];
			}
			else
			{
				left = reverse(readPaths[i].path[j+1]);
				right = reverse(readPaths[i].path[j-1]);
			}
			tripletCoverage[std::make_pair(left, right)] += readPaths[i].reads.size();
			assert(!resolvableGraph.unitigRemoved[node]);
			assert(!resolvableGraph.unitigRemoved[left.first]);
			assert(!resolvableGraph.unitigRemoved[right.first]);
			assert(resolvableGraph.edges[left].count(std::make_pair(node, true)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, false)].count(reverse(left)) == 1);
			assert(resolvableGraph.edges[reverse(right)].count(std::make_pair(node, false)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, true)].count(right) == 1);
		}
	}
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> coveredTriplets;
	for (auto pair : tripletCoverage)
	{
		if (pair.second < minCoverage) continue;
		coveredTriplets.push_back(pair.first);
	}
	return coveredTriplets;
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getReadSupportedTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional)
{
	auto coveredTriplets = getRawTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
	if (unconditional && coveredTriplets.size() >= 2) return coveredTriplets;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredInNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		if (pair.first.first != std::numeric_limits<size_t>::max()) coveredInNeighbors.emplace(reverse(pair.first));
		if (pair.second.first != std::numeric_limits<size_t>::max()) coveredOutNeighbors.emplace(pair.second);
	}
	assert(coveredInNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, false)].size());
	assert(coveredOutNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, true)].size());
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> empty;
	if (coveredInNeighbors.size() < resolvableGraph.edges[std::make_pair(node, false)].size()) return empty;
	if (coveredOutNeighbors.size() < resolvableGraph.edges[std::make_pair(node, true)].size()) return empty;
	return coveredTriplets;
}

double getCoverage(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const size_t unitig)
{
	double result = 0;
	for (auto i : resolvableGraph.readsCrossingNode[unitig])
	{
		assert(readPaths[i].path.size() > 0);
		assert(readPaths[i].reads.size() > 0);
		if (readPaths[i].path.size() == 1)
		{
			assert(readPaths[i].path[0].first == unitig);
			for (const auto& read : readPaths[i].reads)
			{
				assert(read.leftClip + read.rightClip < resolvableGraph.unitigs[readPaths[i].path[0].first].size());
				result += (double)(resolvableGraph.unitigs[readPaths[i].path[0].first].size() - read.leftClip - read.rightClip) / (double)(resolvableGraph.unitigs[readPaths[i].path[0].first].size());
			}
		}
		else
		{
			for (size_t j = 0; j < readPaths[i].path.size(); j++)
			{
				if (readPaths[i].path[j].first != unitig) continue;
				if (j == 0)
				{
					for (const auto& read : readPaths[i].reads)
					{
						assert(read.leftClip < resolvableGraph.unitigs[unitig].size());
						result += (double)(resolvableGraph.unitigs[unitig].size() - read.leftClip) / (double)(resolvableGraph.unitigs[unitig].size());
					}
				}
				if (j > 0 && j < readPaths[i].path.size()-1) result += readPaths[i].reads.size();
				if (j == readPaths[i].path.size()-1)
				{
					for (const auto& read : readPaths[i].reads)
					{
						assert(read.rightClip < resolvableGraph.unitigs[unitig].size());
						result += (double)(resolvableGraph.unitigs[unitig].size() - read.rightClip) / (double)(resolvableGraph.unitigs[unitig].size());
					}
				}
			}
		}
	}
	return result;
}

size_t getAnchorSize(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	size_t result = 0;
	for (auto readi : resolvableGraph.readsCrossingNode[from.first])
	{
		if (resolvableGraph.readsCrossingNode[to.first].count(readi) == 0) continue;
		for (size_t i = 1; i < readPaths[readi].path.size(); i++)
		{
			if (readPaths[readi].path[i-1] == from && readPaths[readi].path[i] == to)
			{
				if (i < readPaths[readi].path.size()-1)
				{
					return resolvableGraph.unitigs[to.first].size();
				}
				for (const auto& path : readPaths[readi].reads)
				{
					assert(path.rightClip < resolvableGraph.unitigs[to.first].size());
					result = std::max(result, resolvableGraph.unitigs[to.first].size() - path.rightClip);
				}
			}
			if (readPaths[readi].path[i-1] == reverse(to) && readPaths[readi].path[i] == reverse(from))
			{
				if (i > 0)
				{
					return resolvableGraph.unitigs[to.first].size();
				}
				for (const auto& path : readPaths[readi].reads)
				{
					assert(path.leftClip < resolvableGraph.unitigs[to.first].size());
					result = std::max(result, resolvableGraph.unitigs[to.first].size() - path.leftClip);
				}
			}
		}
	}
	return result;
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getGuessworkTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional)
{
	size_t fwDegree = resolvableGraph.edges[std::make_pair(node, true)].size();
	size_t bwDegree = resolvableGraph.edges[std::make_pair(node, false)].size();
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> empty;
	if (fwDegree == 1 || bwDegree == 1) return empty;
	double nodeCoverage = getCoverage(resolvableGraph, readPaths, node);
	size_t nodeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
	if (nodeCopyCount < 2) return empty;
	size_t outneighborCopyCountSum = 0;
	size_t inneighborCopyCountSum = 0;
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> outneighborCopyCounts;
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> inneighborCopyCounts;
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		double nodeCoverage = getCoverage(resolvableGraph, readPaths, edge.first);
		size_t edgeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
		outneighborCopyCountSum += edgeCopyCount;
		outneighborCopyCounts.emplace_back(edge, edgeCopyCount);
		if (resolvables.count(edge.first) == 1) return empty;
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		double nodeCoverage = getCoverage(resolvableGraph, readPaths, edge.first);
		size_t edgeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
		inneighborCopyCountSum += edgeCopyCount;
		inneighborCopyCounts.emplace_back(edge, edgeCopyCount);
		if (resolvables.count(edge.first) == 1) return empty;
	}
	if (outneighborCopyCountSum != nodeCopyCount) return empty;
	if (inneighborCopyCountSum != nodeCopyCount) return empty;
	auto coveredTriplets = getRawTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
	if (coveredTriplets.size() == 0) return empty;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredInNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		if (pair.first.first != std::numeric_limits<size_t>::max()) coveredInNeighbors.emplace(reverse(pair.first));
		if (pair.second.first != std::numeric_limits<size_t>::max()) coveredOutNeighbors.emplace(pair.second);
	}
	assert(coveredInNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, false)].size());
	assert(coveredOutNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, true)].size());
	// todo: this could be handled but it's hard to implement
	if (resolvableGraph.edges[std::make_pair(node, false)].size() >= coveredInNeighbors.size() + 2 || resolvableGraph.edges[std::make_pair(node, true)].size() >= coveredOutNeighbors.size() + 2) return empty;
	size_t uncoveredOutCopycounts = 0;
	size_t uncoveredInCopycounts = 0;
	for (auto pair : outneighborCopyCounts)
	{
		if (coveredOutNeighbors.count(pair.first) == 0) uncoveredOutCopycounts += pair.second;
	}
	for (auto pair : inneighborCopyCounts)
	{
		if (coveredInNeighbors.count(pair.first) == 0) uncoveredInCopycounts += pair.second;
	}
	if (uncoveredInCopycounts != uncoveredOutCopycounts) return empty;
	if (uncoveredOutCopycounts == 0) return coveredTriplets;
	if (uncoveredOutCopycounts >= 4) return empty;
	size_t shortestOutAnchor = resolvableGraph.unitigs[node].size();
	size_t shortestInAnchor = resolvableGraph.unitigs[node].size();
	for (auto pair : outneighborCopyCounts)
	{
		if (coveredOutNeighbors.count(pair.first) == 1) continue;
		if (pair.second == 0) continue;
		size_t anchorSize = getAnchorSize(resolvableGraph, readPaths, reverse(pair.first), std::make_pair(node, false));
		shortestOutAnchor = std::min(shortestOutAnchor, anchorSize);
	}
	for (auto pair : inneighborCopyCounts)
	{
		if (coveredOutNeighbors.count(pair.first) == 1) continue;
		if (pair.second == 0) continue;
		size_t anchorSize = getAnchorSize(resolvableGraph, readPaths, reverse(pair.first), std::make_pair(node, true));
		shortestInAnchor = std::min(shortestInAnchor, anchorSize);
	}
	if (shortestInAnchor + shortestOutAnchor <= resolvableGraph.unitigs[node].size()) return empty;
	size_t addedGuesses = 0;
	for (auto inpair : inneighborCopyCounts)
	{
		if (inpair.second == 0) continue;
		if (coveredInNeighbors.count(inpair.first) == 1) continue;
		for (auto outpair : outneighborCopyCounts)
		{
			if (outpair.second == 0) continue;
			if (coveredOutNeighbors.count(outpair.first) == 1) continue;
			coveredTriplets.emplace_back(reverse(inpair.first), outpair.first);
			addedGuesses += 1;
		}
	}
	assert(addedGuesses > 0);
	return coveredTriplets;
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getValidTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional, bool guesswork)
{
	auto triplets = getReadSupportedTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional);
	if (triplets.size() == 0 && guesswork)
	{
		triplets = getGuessworkTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional);;
	}
	std::sort(triplets.begin(), triplets.end(), [](std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> left, std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> right) {
		auto leftComp = canon(left.first, left.second);
		auto rightComp = canon(right.first, right.second);
		assert(leftComp != rightComp);
		if (leftComp.first.first < rightComp.first.first) return true;
		if (leftComp.first.first > rightComp.first.first) return false;
		if (!leftComp.first.second && rightComp.first.second) return true;
		if (leftComp.first.second && !rightComp.first.second) return false;
		if (leftComp.second.first < rightComp.second.first) return true;
		if (leftComp.second.first > rightComp.second.first) return false;
		if (!leftComp.second.second && rightComp.second.second) return true;
		if (leftComp.second.second && !rightComp.second.second) return false;
		assert(false);
		return false;
	});
	return triplets;
}

void replacePaths(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const BigVectorSet& actuallyResolvables, const phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes)
{
	phmap::flat_hash_set<size_t> relevantReads;
	for (const auto node : actuallyResolvables)
	{
		relevantReads.insert(resolvableGraph.readsCrossingNode[node].begin(), resolvableGraph.readsCrossingNode[node].end());
	}
	for (const size_t i : relevantReads)
	{
		PathGroup newPath;
		std::vector<size_t> nodePosStarts;
		std::vector<size_t> nodePosEnds;
		size_t runningKmerStartPos = 0;
		size_t runningKmerEndPos = 0;
		std::vector<size_t> breakFromInvalidEdge;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			runningKmerStartPos = runningKmerEndPos;
			runningKmerEndPos += resolvableGraph.unitigs[readPaths[i].path[j].first].size();
			size_t overlap = 0;
			if (j > 0)
			{
				overlap = resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
				runningKmerEndPos -= overlap;
			}
			if (!actuallyResolvables.get(readPaths[i].path[j].first))
			{
				newPath.path.push_back(readPaths[i].path[j]);
				size_t start = 0;
				if (j > 0)
				{
					start = runningKmerStartPos;
					// assert(start == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j }));
					assert(start >= overlap);
					start -= overlap;
				}
				nodePosStarts.push_back(start);
				size_t end = runningKmerEndPos;
				// assert(end == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 }));
				nodePosEnds.push_back(end);
				continue;
			}
			if (j > 0 && j < readPaths[i].path.size()-1 && !actuallyResolvables.get(readPaths[i].path[j-1].first) && !actuallyResolvables.get(readPaths[i].path[j+1].first))
			{
				if (newEdgeNodes.count(std::make_pair(reverse(readPaths[i].path[j]), reverse(readPaths[i].path[j-1]))) == 0 && newEdgeNodes.count(std::make_pair(readPaths[i].path[j-1], readPaths[i].path[j])) == 0)
				{
					if (newEdgeNodes.count(std::make_pair(readPaths[i].path[j], readPaths[i].path[j+1])) == 0 && newEdgeNodes.count(std::make_pair(reverse(readPaths[i].path[j+1]), reverse(readPaths[i].path[j]))) == 0)
					{
						breakFromInvalidEdge.push_back(newPath.path.size());
						continue;
					}
				}
			}
			if (j > 0 && newEdgeNodes.count(std::make_pair(reverse(readPaths[i].path[j]), reverse(readPaths[i].path[j-1]))) == 1)
			{
				size_t start = runningKmerStartPos;
				// assert(start == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j }));
				if (!actuallyResolvables.get(readPaths[i].path[j-1].first))
				{
					assert(start >= overlap);
					start -= overlap;
				}
				else
				{
					assert(start >= resolvableGraph.unitigs[readPaths[i].path[j-1].first].size());
					start -= resolvableGraph.unitigs[readPaths[i].path[j-1].first].size();
				}
				nodePosStarts.push_back(start);
				size_t end = runningKmerEndPos;
				// assert(end == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 }));
				nodePosEnds.push_back(end);
				newPath.path.emplace_back(newEdgeNodes.at(std::make_pair(reverse(readPaths[i].path[j]), reverse(readPaths[i].path[j-1]))), false);
			}
			if (j < readPaths[i].path.size()-1 && newEdgeNodes.count(std::make_pair(readPaths[i].path[j], readPaths[i].path[j+1])) == 1)
			{
				size_t start = 0;
				if (j > 0)
				{
					start = runningKmerStartPos;
					// assert(start == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j }));
					assert(start >= overlap);
					start -= overlap;
				}
				nodePosStarts.push_back(start);
				size_t end = runningKmerEndPos;
				// assert(end == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 }));
				if (actuallyResolvables.get(readPaths[i].path[j+1].first))
				{
					end += resolvableGraph.unitigs[readPaths[i].path[j+1].first].size();
					end -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j], readPaths[i].path[j+1]));
					// assert(end == getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 2 }));
				}
				nodePosEnds.push_back(end);
				newPath.path.emplace_back(newEdgeNodes.at(std::make_pair(readPaths[i].path[j], readPaths[i].path[j+1])), true);
			}
		}
		std::reverse(breakFromInvalidEdge.begin(), breakFromInvalidEdge.end());
		size_t kmerPathLength = runningKmerEndPos;
		assert(kmerPathLength == getNumberOfHashes(resolvableGraph, 0, 0, readPaths[i].path));
		if (newPath.path.size() == 0)
		{
			erasePath(resolvableGraph, readPaths, i);
			continue;
		}
		std::swap(readPaths[i].reads, newPath.reads);
		assert(nodePosStarts.size() == nodePosEnds.size());
		assert(nodePosEnds.size() == newPath.path.size());
		assert(nodePosStarts[0] >= 0);
		assert(nodePosEnds.back() <= kmerPathLength);
		if (nodePosStarts[0] != 0 || nodePosEnds.back() != kmerPathLength)
		{
			size_t startRemove = nodePosStarts[0];
			for (auto& read : newPath.reads)
			{
				if (nodePosStarts[0] > read.leftClip)
				{
					read.readPoses.erase(read.readPoses.begin(), read.readPoses.begin() + nodePosStarts[0] - read.leftClip);
					read.leftClip = 0;
				}
				else
				{
					assert(read.leftClip >= nodePosStarts[0]);
					read.leftClip -= nodePosStarts[0];
				}
				if (nodePosEnds.back() < kmerPathLength - read.rightClip)
				{
					size_t extraClip = (kmerPathLength - read.rightClip) - nodePosEnds.back();
					read.readPoses.erase(read.readPoses.begin() + read.readPoses.size() - extraClip, read.readPoses.begin() + read.readPoses.size());
					read.rightClip = 0;
				}
				else
				{
					assert(read.rightClip >= (kmerPathLength - nodePosEnds.back()));
					read.rightClip -= (kmerPathLength - nodePosEnds.back());
				}
				assert(nodePosEnds.back() - nodePosStarts[0] == read.readPoses.size() + read.leftClip + read.rightClip);
			}
			for (size_t j = 0; j < nodePosStarts.size(); j++)
			{
				assert(nodePosStarts[j] >= startRemove);
				nodePosStarts[j] -= startRemove;
				assert(nodePosEnds[j] >= startRemove);
				nodePosEnds[j] -= startRemove;
			}
		}
		assert(nodePosStarts[0] == 0);
		size_t lastStart = 0;
		for (size_t j = 1; j < newPath.path.size(); j++)
		{
			assert(resolvableGraph.edges[newPath.path[j-1]].count(newPath.path[j]) == resolvableGraph.edges[reverse(newPath.path[j])].count(reverse(newPath.path[j-1])));
			assert(breakFromInvalidEdge.size() == 0 || breakFromInvalidEdge.back() >= j);
			if (resolvableGraph.edges[newPath.path[j-1]].count(newPath.path[j]) == 0 || (breakFromInvalidEdge.size() > 0 && breakFromInvalidEdge.back() == j))
			{
				if (breakFromInvalidEdge.size() > 0 && breakFromInvalidEdge.back() == j)
				{
					breakFromInvalidEdge.pop_back();
				}
				PathGroup path;
				path.path.insert(path.path.end(), newPath.path.begin() + lastStart, newPath.path.begin() + j);
				path.reads.reserve(newPath.reads.size());
				for (const auto& read : newPath.reads)
				{
					size_t posesStart = nodePosStarts[lastStart];
					size_t posesEnd = nodePosEnds[j-1];
					size_t leftClipRemove = 0;
					size_t rightClipRemove = 0;
					if (posesEnd < read.leftClip + read.readPoses.size())
					{
						rightClipRemove = read.rightClip;
						posesEnd = read.readPoses.size() - (read.leftClip + read.readPoses.size() - posesEnd);
					}
					else if (posesEnd < read.leftClip + read.rightClip + read.readPoses.size())
					{
						rightClipRemove = (read.leftClip + read.rightClip + read.readPoses.size()) - posesEnd;
						posesEnd = read.readPoses.size();
					}
					if (posesStart > read.leftClip)
					{
						leftClipRemove = read.leftClip;
						posesStart = posesStart - read.leftClip;
					}
					else if (posesStart > 0)
					{
						leftClipRemove = posesStart;
						posesStart = 0;
					}
					assert(posesStart < posesEnd);
					assert(posesEnd <= read.readPoses.size());
					path.reads.emplace_back();
					path.reads.back().readPoses.insert(path.reads.back().readPoses.end(), read.readPoses.begin() + posesStart, read.readPoses.begin() + posesEnd);
					path.reads.back().leftClip = read.leftClip;
					path.reads.back().rightClip = read.rightClip;
					assert(path.reads.back().leftClip >= leftClipRemove);
					path.reads.back().leftClip -= leftClipRemove;
					assert(path.reads.back().rightClip >= rightClipRemove);
					path.reads.back().rightClip -= rightClipRemove;
					path.reads.back().readNameIndex = read.readNameIndex;
					path.reads.back().readLength = read.readLength;
					path.reads.back().readLengthHPC = read.readLengthHPC;
				}
				addPath(resolvableGraph, readPaths, std::move(path));
				lastStart = j;
			}
		}
		PathGroup path;
		path.path.insert(path.path.end(), newPath.path.begin() + lastStart, newPath.path.end());
		path.reads.reserve(newPath.reads.size());
		for (const auto& read : newPath.reads)
		{
			size_t posesStart = nodePosStarts[lastStart];
			size_t posesEnd = nodePosEnds.back();
			size_t leftClipRemove = 0;
			size_t rightClipRemove = 0;
			assert(posesEnd <= read.leftClip + read.rightClip + read.readPoses.size());
			if (posesEnd < read.leftClip + read.readPoses.size())
			{
				rightClipRemove = read.rightClip;
				posesEnd = read.readPoses.size() - (read.leftClip + read.readPoses.size() - posesEnd);
			}
			else
			{
				rightClipRemove = (read.leftClip + read.rightClip + read.readPoses.size()) - posesEnd;
				posesEnd = read.readPoses.size();
			}
			if (posesStart > read.leftClip)
			{
				leftClipRemove = read.leftClip;
				posesStart = posesStart - read.leftClip;
			}
			else if (posesStart > 0)
			{
				leftClipRemove = posesStart;
				posesStart = 0;
			}
			assert(posesStart < posesEnd);
			assert(posesEnd <= read.readPoses.size());
			path.reads.emplace_back();
			path.reads.back().readPoses.insert(path.reads.back().readPoses.end(), read.readPoses.begin() + posesStart, read.readPoses.begin() + posesEnd);
			path.reads.back().leftClip = read.leftClip;
			path.reads.back().rightClip = read.rightClip;
			assert(path.reads.back().leftClip >= leftClipRemove);
			path.reads.back().leftClip -= leftClipRemove;
			assert(path.reads.back().rightClip >= rightClipRemove);
			path.reads.back().rightClip -= rightClipRemove;
			path.reads.back().readNameIndex = read.readNameIndex;
			path.reads.back().readLength = read.readLength;
			path.reads.back().readLengthHPC = read.readLengthHPC;
		}
		addPath(resolvableGraph, readPaths, std::move(path));
		erasePath(resolvableGraph, readPaths, i);
	}
}

struct ResolutionResult
{
public:
	size_t nodesResolved;
	size_t nodesAdded;
	phmap::flat_hash_set<size_t> maybeUnitigifiable;
	phmap::flat_hash_map<std::pair<size_t, bool>, size_t> maybeTrimmable;
};

ResolutionResult resolve(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, const size_t kmerSize, std::vector<PathGroup>& readPaths, const phmap::flat_hash_set<size_t>& resolvables, const size_t minCoverage, const bool unconditional, const bool guesswork)
{
	static BigVectorSet actuallyResolvables;
	ResolutionResult result;
	result.nodesResolved = 0;
	result.nodesAdded = 0;
	phmap::flat_hash_set<size_t> unresolvables;
	for (auto node : resolvables)
	{
		for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
		{
			if (edge.first == node && !edge.second)
			{
				std::cout << "unresolvable even length exact palindrome frozen, id: " << node << std::endl;
				unresolvables.insert(node);
			}
		}
		for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
		{
			if (edge.first == node && edge.second)
			{
				std::cout << "unresolvable even length exact palindrome frozen, id: " << node << std::endl;
				unresolvables.insert(node);
			}
		}
	}
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork);
		if (triplets.size() == 0)
		{
			unresolvables.insert(node);
			// unresolveRecursively(resolvableGraph, resolvables, unresolvables, node);
		}
	}
	std::vector<size_t> check;
	check.insert(check.end(), resolvables.begin(), resolvables.end());
	while (true)
	{
		std::vector<size_t> removeThese;
		phmap::flat_hash_set<size_t> newCheck;
		for (auto node : check)
		{
			if (unresolvables.count(node) == 1) continue;
			auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork);
			for (auto triplet : triplets)
			{
				if (triplet.first.first == std::numeric_limits<size_t>::max()) continue;
				if (triplet.second.first == std::numeric_limits<size_t>::max()) continue;
				if (resolvables.count(triplet.first.first) == 0 || unresolvables.count(triplet.first.first) == 1)
				{
					if (resolvables.count(triplet.second.first) == 0 || unresolvables.count(triplet.second.first) == 1)
					{
						if (resolvableGraph.unitigLength(triplet.first.first) == resolvableGraph.getBpOverlap(std::make_pair(node, false), reverse(triplet.first))+1)
						{
							if (resolvableGraph.unitigLength(triplet.second.first) == resolvableGraph.getBpOverlap(std::make_pair(node, true), triplet.second)+1)
							{
								assert(unresolvables.count(node) == 0);
								unresolvables.insert(node);
								removeThese.emplace_back(node);
								for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
								{
									if (resolvables.count(edge.first) == 1 && unresolvables.count(edge.first) == 0) newCheck.insert(edge.first);
								}
								for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
								{
									if (resolvables.count(edge.first) == 1 && unresolvables.count(edge.first) == 0) newCheck.insert(edge.first);
								}
								break;
							}
						}
					}
				}
			}
		}
		check.clear();
		check.insert(check.end(), newCheck.begin(), newCheck.end());
		if (newCheck.size() == 0) break;
	}
	if (unresolvables.size() == resolvables.size()) return result;
	actuallyResolvables.resize(resolvableGraph.unitigs.size());
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		actuallyResolvables.set(node);
	}
	assert(actuallyResolvables.activeSize() == resolvables.size() - unresolvables.size());
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newEdgeNodes;
	for (auto node : actuallyResolvables)
	{
		assert(getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork).size() > 0);
		std::unordered_set<std::pair<size_t, bool>> fwCovered;
		std::unordered_set<std::pair<size_t, bool>> bwCovered;
		for (auto pair : getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork))
		{
			fwCovered.insert(pair.second);
			bwCovered.insert(reverse(pair.first));
		}
		std::pair<size_t, bool> pos { node, true };
		std::vector<std::pair<size_t, bool>> fwEdges { resolvableGraph.edges[pos].begin(), resolvableGraph.edges[pos].end() };
		std::vector<std::pair<size_t, bool>> bwEdges { resolvableGraph.edges[reverse(pos)].begin(), resolvableGraph.edges[reverse(pos)].end() };
		std::sort(fwEdges.begin(), fwEdges.end());
		std::sort(bwEdges.begin(), bwEdges.end());
		for (auto edge : fwEdges)
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if (fwCovered.count(edge) == 0)
			{
				result.maybeUnitigifiable.insert(edge.first);
				result.maybeTrimmable[reverse(edge)] = std::max(resolvableGraph.overlaps.at(canon(pos, edge)), result.maybeTrimmable[reverse(edge)]);
				continue;
			}
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1) continue;
			assert(actuallyResolvables.get(edge.first) || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, hashlist, kmerSize, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
			if (bwEdges.size() == 0)
			{
				result.maybeUnitigifiable.insert(newEdgeNodes.at(std::make_pair(pos, edge)));
				result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(pos, edge)), false)] = resolvableGraph.unitigs[pos.first].size();
			}
		}
		pos = std::make_pair(node, false);
		for (auto edge : bwEdges)
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if (bwCovered.count(edge) == 0)
			{
				result.maybeUnitigifiable.insert(edge.first);
				result.maybeTrimmable[reverse(edge)] = std::max(resolvableGraph.overlaps.at(canon(pos, edge)), result.maybeTrimmable[reverse(edge)]);
				continue;
			}
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1) continue;
			assert(actuallyResolvables.get(edge.first) || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, hashlist, kmerSize, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
			if (fwEdges.size() == 0)
			{
				result.maybeUnitigifiable.insert(newEdgeNodes.at(std::make_pair(pos, edge)));
				result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(pos, edge)), false)] = resolvableGraph.unitigs[pos.first].size();
			}
		}
	}
	if (unconditional)
	{
		for (auto pair : newEdgeNodes)
		{
			auto fromnode = pair.first.first;
			auto tonode = pair.first.second;
			if (!actuallyResolvables.get(fromnode.first) || !actuallyResolvables.get(tonode.first)) continue;
			bool fromHasTo = false;
			bool toHasFrom = false;
			for (auto triplet : getValidTriplets(resolvableGraph, resolvables, readPaths, fromnode.first, minCoverage, unconditional, guesswork))
			{
				if (fromnode.second && triplet.second == tonode)
				{
					fromHasTo = true;
				}
				if (!fromnode.second && reverse(triplet.first) == tonode)
				{
					fromHasTo = true;
				}
			}
			for (auto triplet : getValidTriplets(resolvableGraph, resolvables, readPaths, tonode.first, minCoverage, unconditional, guesswork))
			{
				if (tonode.second && triplet.second == reverse(fromnode))
				{
					toHasFrom = true;
				}
				if (!tonode.second && reverse(triplet.first) == reverse(fromnode))
				{
					toHasFrom = true;
				}
			}
			assert(fromHasTo || toHasFrom);
			if (fromHasTo && toHasFrom) continue;
			if (fromHasTo)
			{
				assert(!toHasFrom);
				result.maybeTrimmable[std::make_pair(pair.second, true)] = resolvableGraph.unitigs[tonode.first].size();
			}
			if (toHasFrom)
			{
				assert(!fromHasTo);
				result.maybeTrimmable[std::make_pair(pair.second, false)] = resolvableGraph.unitigs[fromnode.first].size();
			}
		}
	}
	for (auto node : actuallyResolvables)
	{
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork);
		assert(triplets.size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto triplet : triplets)
		{
			const std::pair<size_t, bool> before = triplet.first;
			const std::pair<size_t, bool> after = triplet.second;
			assert(before.first != std::numeric_limits<size_t>::max() || after.first != std::numeric_limits<size_t>::max());
			if (before.first == std::numeric_limits<size_t>::max())
			{
				assert(after.first != std::numeric_limits<size_t>::max());
				assert(newEdgeNodes.count(std::make_pair(pos, after)) == 1 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1 || !actuallyResolvables.get(after.first));
				if (newEdgeNodes.count(std::make_pair(pos, after)) == 1)
				{
					result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(pos, after)), false)] = resolvableGraph.unitigs[node].size();
				}
				else if (newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1)
				{
					result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(reverse(after), reverse(pos))), true)] = resolvableGraph.unitigs[node].size();
				}
				else
				{
					assert(!actuallyResolvables.get(after.first));
					assert(!resolvableGraph.unitigRemoved[after.first]);
					result.maybeTrimmable[reverse(after)] = resolvableGraph.unitigs[after.first].size();
				}
				continue;
			}
			if (after.first == std::numeric_limits<size_t>::max())
			{
				assert(before.first != std::numeric_limits<size_t>::max());
				assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1 || newEdgeNodes.count(std::make_pair(before, pos)) == 1 || !actuallyResolvables.get(before.first));
				if (newEdgeNodes.count(std::make_pair(before, pos)) == 1)
				{
					result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(before, pos)), true)] = resolvableGraph.unitigs[node].size();
				}
				else if (newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1)
				{
					result.maybeTrimmable[std::make_pair(newEdgeNodes.at(std::make_pair(reverse(pos), reverse(before))), false)] = resolvableGraph.unitigs[node].size();
				}
				else
				{
					assert(!actuallyResolvables.get(before.first));
					assert(!resolvableGraph.unitigRemoved[before.first]);
					result.maybeTrimmable[before] = resolvableGraph.unitigs[before.first].size();
				}
				continue;
			}
			std::pair<size_t, bool> leftNode = before;
			std::pair<size_t, bool> rightNode = after;
			size_t overlap = resolvableGraph.unitigs[node].size();
			bool foundLeft = false;
			bool foundRight = false;
			assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 0 || newEdgeNodes.count(std::make_pair(before, pos)) == 0);
			assert(newEdgeNodes.count(std::make_pair(pos, after)) == 0 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(before, pos)) == 0);
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(pos), reverse(before))), false);
				foundLeft = true;
			}
			else if (newEdgeNodes.count(std::make_pair(before, pos)) == 1)
			{
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(before, pos)), true);
				foundLeft = true;
			}
			else
			{
				assert(resolvables.count(before.first) == 0 || unresolvables.count(before.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(pos, after)) == 1 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1);
				assert(resolvableGraph.unitigLength(before.first) == resolvableGraph.getBpOverlap(before, pos)+1);
				overlap = resolvableGraph.overlaps.at(canon(before, pos));
			}
			if (newEdgeNodes.count(std::make_pair(pos, after)) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(pos, after)), true);
				foundRight = true;
			}
			else if (newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1)
			{
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(after), reverse(pos))), false);
				foundRight = true;
			}
			else
			{
				assert(resolvables.count(after.first) == 0 || unresolvables.count(after.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1 || newEdgeNodes.count(std::make_pair(before, pos)) == 1);
				assert(resolvableGraph.unitigLength(after.first) == resolvableGraph.getBpOverlap(pos, after)+1);
				overlap = resolvableGraph.overlaps.at(canon(pos, after));
			}
			result.maybeUnitigifiable.insert(leftNode.first);
			result.maybeUnitigifiable.insert(rightNode.first);
			assert(foundLeft || foundRight);
			assert(!resolvableGraph.unitigRemoved[leftNode.first]);
			assert(!resolvableGraph.unitigRemoved[rightNode.first]);
			assert(resolvables.count(leftNode.first) == 0 || unresolvables.count(leftNode.first) == 1);
			assert(resolvables.count(rightNode.first) == 0 || unresolvables.count(rightNode.first) == 1);
			assert(resolvableGraph.edges[leftNode].count(rightNode) == 0);
			assert(resolvableGraph.edges[reverse(rightNode)].count(reverse(leftNode)) == 0);
			resolvableGraph.edges[leftNode].emplace(rightNode);
			resolvableGraph.edges[reverse(rightNode)].emplace(reverse(leftNode));
			resolvableGraph.overlaps[canon(leftNode, rightNode)] = overlap;
			auto leftSeq = resolvableGraph.unitigs[leftNode.first];
			if (!leftNode.second) leftSeq = revCompPath(leftSeq);
			auto rightSeq = resolvableGraph.unitigs[rightNode.first];
			if (!rightNode.second) rightSeq = revCompPath(rightSeq);
			assert(leftSeq.size() >= overlap);
			assert(rightSeq.size() >= overlap);
			for (size_t i = 0; i < overlap; i++)
			{
				assert(leftSeq[leftSeq.size()-overlap+i] == rightSeq[i]);
			}
			assert(resolvableGraph.edges[leftNode].count(rightNode) == 1);
			assert(resolvableGraph.edges[reverse(rightNode)].count(reverse(leftNode)) == 1);
			assert(resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(leftNode.first));
			assert(resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(rightNode.first));
		}
	}
	for (auto node : actuallyResolvables)
	{
		resolvableGraph.unitigRemoved[node] = true;
		std::pair<size_t, bool> fw { node, true };
		std::pair<size_t, bool> bw { node, false };
		for (auto edge : resolvableGraph.edges[fw])
		{
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(fw)) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(fw));
		}
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(bw)) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(bw));
		}
		resolvableGraph.edges[fw].clear();
		resolvableGraph.edges[bw].clear();
	}
	replacePaths(resolvableGraph, readPaths, actuallyResolvables, newEdgeNodes);
	actuallyResolvables.clear();
	assert(resolvables.size() > unresolvables.size());
	result.nodesResolved = resolvables.size() - unresolvables.size();
	result.nodesAdded = newEdgeNodes.size();
	return result;
}

void checkValidity(const ResolvableUnitigGraph& graph, const std::vector<PathGroup>& readPaths, const size_t kmerSize)
{
	return;
	assert(graph.unitigs.size() == graph.edges.size());
	assert(graph.unitigs.size() == graph.unitigRightClipBp.size());
	assert(graph.unitigs.size() == graph.unitigLeftClipBp.size());
	assert(graph.unitigs.size() == graph.unitigRemoved.size());
	assert(graph.unitigs.size() == graph.readsCrossingNode.size());
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		if (graph.unitigRemoved[i])
		{
			assert(graph.edges[std::make_pair(i, true)].size() == 0);
			assert(graph.edges[std::make_pair(i, false)].size() == 0);
			continue;
		}
		if (graph.unitigs[i].size() >= 2)
		{
			assert(graph.unitigLeftClipBp[i] < kmerSize - graph.hashlist.getOverlap(graph.unitigs[i][0], graph.unitigs[i][1]));
			assert(graph.unitigRightClipBp[i] < kmerSize - graph.hashlist.getOverlap(graph.unitigs[i][graph.unitigs[i].size()-2], graph.unitigs[i][graph.unitigs[i].size()-1]));
		}
		assert(graph.unitigLeftClipBp[i] < kmerSize);
		assert(graph.unitigRightClipBp[i] < kmerSize);
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges[fw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(fw)) == 1);
			assert(graph.edges[fw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2 || edge.first == i);
			assert(graph.getBpOverlap(fw, edge) < graph.unitigLength(i));
			assert(graph.getBpOverlap(fw, edge) < graph.unitigLength(edge.first));
		}
		for (auto edge : graph.edges[bw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(bw)) == 1);
			assert(graph.edges[bw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2 || edge.first == i);
			assert(graph.getBpOverlap(bw, edge) < graph.unitigLength(i));
			assert(graph.getBpOverlap(bw, edge) < graph.unitigLength(edge.first));
		}
	}
	for (const auto& path : readPaths)
	{
		for (size_t i = 0; i < path.path.size(); i++)
		{
			assert(!graph.unitigRemoved[path.path[i].first]);
			if (i > 0) assert(graph.edges[path.path[i-1]].count(path.path[i]) == 1);
		}
		for (const auto& read : path.reads)
		{
			if (read.readPoses.size() > 0)
			{
				for (size_t i = 1; i < read.readPoses.size(); i++)
				{
					assert(read.readPoses[i-1] < read.readPoses[i]);
				}
				assert(read.readPoses.back() + kmerSize <= read.readLengthHPC);
			}
		}
	}
	std::vector<std::vector<size_t>> coverages;
	coverages.resize(graph.unitigs.size());
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		coverages[i].resize(graph.unitigs[i].size(), 0);
	}
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			assert(path.path[i].first < graph.unitigs.size());
			assert(!graph.unitigRemoved[path.path[i].first]);
		}
		for (size_t i = 1; i < path.path.size(); i++)
		{
			assert(graph.edges[path.path[i-1]].count(path.path[i]) == 1);
			assert(graph.edges[reverse(path.path[i])].count(reverse(path.path[i-1])) == 1);
		}
		size_t pathHashCount = getNumberOfHashes(graph, 0, 0, path.path);
		if (path.path.size() == 1)
		{
			for (const auto& read : path.reads)
			{
				assert(read.leftClip + read.rightClip + read.readPoses.size() == pathHashCount);
				assert(read.leftClip + read.rightClip < graph.unitigs[path.path[0].first].size());
				for (size_t i = read.leftClip; i < graph.unitigs[path.path[0].first].size() - read.rightClip; i++)
				{
					size_t index = i;
					if (!path.path[0].second) index = graph.unitigs[path.path[0].first].size() - i - 1;
					coverages[path.path[0].first][index] += 1;
				}
			}
			continue;
		}
		assert(path.path.size() >= 2);
		for (size_t i = 1; i < path.path.size()-1; i++)
		{
			for (size_t j = 0; j < graph.unitigs[path.path[i].first].size(); j++)
			{
				coverages[path.path[i].first][j] += path.reads.size();
			}
		}
		for (const auto& read : path.reads)
		{
			assert(read.leftClip + read.rightClip + read.readPoses.size() == pathHashCount);
			assert(read.leftClip < graph.unitigs[path.path[0].first].size());
			for (size_t i = read.leftClip; i < graph.unitigs[path.path[0].first].size(); i++)
			{
				size_t index = i;
				if (!path.path[0].second) index = graph.unitigs[path.path[0].first].size() - i - 1;
				coverages[path.path[0].first][index] += 1;
			}
			assert(read.rightClip < graph.unitigs[path.path.back().first].size());
			for (size_t i = 0; i < graph.unitigs[path.path.back().first].size() - read.rightClip; i++)
			{
				size_t index = i;
				if (!path.path.back().second) index = graph.unitigs[path.path.back().first].size() - i - 1;
				coverages[path.path.back().first][index] += 1;
			}
		}
	}
	for (size_t i = 0; i < coverages.size(); i++)
	{
		if (graph.unitigRemoved[i]) continue;
		for (size_t j = 0; j < coverages[i].size(); j++)
		{
			assertPrintReads(coverages[i][j] > 0, graph, readPaths, i);
		}
	}
}

size_t getEdgeCoverage(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	std::vector<size_t> relevantReads;
	for (auto read : resolvableGraph.readsCrossingNode[from.first])
	{
		if (resolvableGraph.readsCrossingNode[to.first].count(read) == 1) relevantReads.emplace_back(read);
	}
	size_t result = 0;
	for (size_t i : relevantReads)
	{
		assert(readPaths[i].path.size() > 0);
		assert(readPaths[i].reads.size() > 0);
		for (size_t j = 1; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j-1] == from && readPaths[i].path[j] == to)
			{
				result += readPaths[i].reads.size();
			}
			else if (readPaths[i].path[j-1] == reverse(to) && readPaths[i].path[j] == reverse(from))
			{
				result += readPaths[i].reads.size();
			}
		}
	}
	return result;
}

void removeNode(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, size_t node)
{
	assert(!resolvableGraph.unitigRemoved[node]);
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		if (edge.first == node) continue;
		assert(resolvableGraph.edges[reverse(edge)].count(std::make_pair(node, false)) == 1);
		resolvableGraph.edges[reverse(edge)].erase(std::make_pair(node, false));
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		if (edge.first == node) continue;
		assert(resolvableGraph.edges[reverse(edge)].count(std::make_pair(node, true)) == 1);
		resolvableGraph.edges[reverse(edge)].erase(std::make_pair(node, true));
	}
	resolvableGraph.edges[std::make_pair(node, true)].clear();
	resolvableGraph.edges[std::make_pair(node, false)].clear();
	for (const size_t i : resolvableGraph.readsCrossingNode[node])
	{
		assert(readPaths[i].path.size() != 0);
		if (readPaths[i].reads.size() == 0) continue;
		std::vector<size_t> nodePosStarts;
		std::vector<size_t> nodePosEnds;
		size_t pathKmerLength = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (j > 0) pathKmerLength -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
			nodePosStarts.push_back(pathKmerLength);
			pathKmerLength += resolvableGraph.unitigs[readPaths[i].path[j].first].size();
			nodePosEnds.push_back(pathKmerLength);
		}
		assert(nodePosStarts[0] == 0);
		assert(nodePosEnds.back() == pathKmerLength);
		for (const auto& read : readPaths[i].reads)
		{
			assert(pathKmerLength == read.readPoses.size() + read.leftClip + read.rightClip);
		}
		size_t lastStart = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j].first == node)
			{
				if (j == lastStart)
				{
					lastStart = j+1;
					continue;
				}
				PathGroup path;
				path.path.insert(path.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.begin() + j);
				path.reads.reserve(readPaths[i].reads.size());
				for (const auto& read : readPaths[i].reads)
				{
					size_t posesStart = nodePosStarts[lastStart];
					size_t posesEnd = nodePosEnds[j-1];
					size_t leftClipRemove = 0;
					size_t rightClipRemove = 0;
					if (posesEnd < read.leftClip + read.readPoses.size())
					{
						rightClipRemove = read.rightClip;
						posesEnd = read.readPoses.size() - (read.leftClip + read.readPoses.size() - posesEnd);
					}
					else
					{
						assert(posesEnd <= read.leftClip + read.rightClip + read.readPoses.size());
						rightClipRemove = (read.leftClip + read.rightClip + read.readPoses.size()) - posesEnd;
						posesEnd = read.readPoses.size();
					}
					if (posesStart > read.leftClip)
					{
						leftClipRemove = read.leftClip;
						posesStart = posesStart - read.leftClip;
					}
					else
					{
						leftClipRemove = posesStart;
						posesStart = 0;
					}
					assert(posesStart < posesEnd);
					assert(posesEnd <= read.readPoses.size());
					path.reads.emplace_back();
					path.reads.back().readPoses.insert(path.reads.back().readPoses.end(), read.readPoses.begin() + posesStart, read.readPoses.begin() + posesEnd);
					path.reads.back().leftClip = read.leftClip;
					path.reads.back().rightClip = read.rightClip;
					assert(path.reads.back().leftClip >= leftClipRemove);
					path.reads.back().leftClip -= leftClipRemove;
					assert(path.reads.back().rightClip >= rightClipRemove);
					path.reads.back().rightClip -= rightClipRemove;
					path.reads.back().readNameIndex = read.readNameIndex;
					path.reads.back().readLength = read.readLength;
					path.reads.back().readLengthHPC = read.readLengthHPC;
				}
				addPath(resolvableGraph, readPaths, std::move(path));
				lastStart = j+1;
			}
		}
		assert(lastStart != 0);
		if (lastStart < readPaths[i].path.size())
		{
			PathGroup path;
			path.path.insert(path.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.end());
			path.reads.reserve(readPaths[i].reads.size());
			for (const auto& read : readPaths[i].reads)
			{
				size_t posesStart = nodePosStarts[lastStart];
				size_t posesEnd = nodePosEnds.back();
				size_t leftClipRemove = 0;
				size_t rightClipRemove = 0;
				assert(posesEnd <= read.leftClip + read.rightClip + read.readPoses.size());
				if (posesEnd < read.leftClip + read.readPoses.size())
				{
					rightClipRemove = read.rightClip;
					posesEnd = read.readPoses.size() - (read.leftClip + read.readPoses.size() - posesEnd);
				}
				else
				{
					assert(posesEnd <= read.leftClip + read.rightClip + read.readPoses.size());
					rightClipRemove = (read.leftClip + read.rightClip + read.readPoses.size()) - posesEnd;
					posesEnd = read.readPoses.size();
				}
				if (posesStart > read.leftClip)
				{
					leftClipRemove = read.leftClip;
					posesStart = posesStart - read.leftClip;
				}
				else
				{
					leftClipRemove = posesStart;
					posesStart = 0;
				}
				assert(posesStart < posesEnd);
				assert(posesEnd <= read.readPoses.size());
				path.reads.emplace_back();
				path.reads.back().readPoses.insert(path.reads.back().readPoses.end(), read.readPoses.begin() + posesStart, read.readPoses.begin() + posesEnd);
				path.reads.back().leftClip = read.leftClip;
				path.reads.back().rightClip = read.rightClip;
				assert(path.reads.back().leftClip >= leftClipRemove);
				path.reads.back().leftClip -= leftClipRemove;
				assert(path.reads.back().rightClip >= rightClipRemove);
				path.reads.back().rightClip -= rightClipRemove;
				path.reads.back().readNameIndex = read.readNameIndex;
				path.reads.back().readLength = read.readLength;
				path.reads.back().readLengthHPC = read.readLengthHPC;
			}
			addPath(resolvableGraph, readPaths, std::move(path));
		}
		erasePath(resolvableGraph, readPaths, i);
	}
	assert(resolvableGraph.readsCrossingNode[node].size() == 0);
	resolvableGraph.unitigRemoved[node] = true;
}

struct UntippingResult
{
	size_t nodesRemoved;
	phmap::flat_hash_set<size_t> maybeUnitigifiable;
};

void tryRemoveTip(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const size_t kmerSize, const size_t i, UntippingResult& result)
{
	assert(!resolvableGraph.unitigRemoved[i]);
	size_t maxEdgeCoverage = 0;
	size_t fwHasSafeEdge = false;
	size_t fwHasSafeNode = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)])
	{
		maxEdgeCoverage = std::max(maxEdgeCoverage, getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, true), edge));
		if (getCoverage(resolvableGraph, readPaths, edge.first) >= minSafeCoverage) fwHasSafeNode = true;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) fwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, true)].size() > 0 && (!fwHasSafeNode || !fwHasSafeEdge)) return;
	bool bwHasSafeNode = false;
	bool bwHasSafeEdge = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)])
	{
		maxEdgeCoverage = std::max(maxEdgeCoverage, getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, false), edge));
		if (getCoverage(resolvableGraph, readPaths, edge.first) >= minSafeCoverage) bwHasSafeNode = true;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) bwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, false)].size() > 0 && (!bwHasSafeNode || !bwHasSafeEdge)) return;
	if (maxEdgeCoverage > maxRemovableCoverage) return;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)]) result.maybeUnitigifiable.insert(edge.first);
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)]) result.maybeUnitigifiable.insert(edge.first);
	removeNode(resolvableGraph, readPaths, i);
	result.nodesRemoved += 1;
}

UntippingResult removeLowCoverageTips(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const size_t kmerSize, const phmap::flat_hash_set<size_t>& maybeUntippable)
{
	for (size_t i = resolvableGraph.lastTippableChecked; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = getCoverage(resolvableGraph, readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	for (auto i : maybeUntippable)
	{
		if (i >= resolvableGraph.lastTippableChecked) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = getCoverage(resolvableGraph, readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	resolvableGraph.lastTippableChecked = resolvableGraph.unitigs.size();
	UntippingResult result;
	result.nodesRemoved = 0;
	for (size_t index = resolvableGraph.everTippable.size()-1; index < resolvableGraph.everTippable.size(); index--)
	{
		size_t i = resolvableGraph.everTippable[index];
		if (resolvableGraph.unitigRemoved[i])
		{
			std::swap(resolvableGraph.everTippable[index], resolvableGraph.everTippable.back());
			resolvableGraph.everTippable.pop_back();
			continue;
		}
		tryRemoveTip(resolvableGraph, readPaths, hashlist, maxRemovableCoverage, minSafeCoverage, maxRemovableLength, kmerSize, i, result);
	}
	return result;
}

bool canTrimRecursive(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> pos, const size_t trimAmount)
{
	assert(trimAmount > 0);
	assert(resolvableGraph.unitigs[pos.first].size() > trimAmount);
	for (auto edge : resolvableGraph.edges[pos])
	{
		if (resolvableGraph.overlaps.at(canon(pos, edge)) < trimAmount) return false;
	}
	for (auto edge : resolvableGraph.edges[reverse(pos)])
	{
		size_t overlap = resolvableGraph.overlaps.at(canon(reverse(pos), edge));
		if (overlap >= resolvableGraph.unitigs[pos.first].size() - trimAmount)
		{
			// todo: this might actually be trimmable, but it's hard so don't try for now
			if (overlap == resolvableGraph.unitigs[pos.first].size()) return false;
			assertPrintReads(overlap < resolvableGraph.unitigs[pos.first].size(), resolvableGraph, readPaths, pos.first);
			size_t trimThere = overlap - (resolvableGraph.unitigs[pos.first].size() - trimAmount) + 1;
			assert(trimThere <= trimAmount);
			assert(trimThere < resolvableGraph.unitigs[edge.first].size());
			assert(trimThere > 0);
			if (!canTrimRecursive(resolvableGraph, readPaths, reverse(edge), trimThere)) return false;
		}
	}
	size_t maxReadTrim = resolvableGraph.unitigs[pos.first].size();
	for (const size_t i : resolvableGraph.readsCrossingNode[pos.first])
	{
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j].first != pos.first) continue;
			if (j == 0 && readPaths[i].path[j] == reverse(pos))
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					maxReadTrim = std::min(maxReadTrim, readPaths[i].reads[k].leftClip);
				}
			}
			else if (j == readPaths[i].path.size()-1 && readPaths[i].path[j] == pos)
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					maxReadTrim = std::min(maxReadTrim, readPaths[i].reads[k].rightClip);
				}
			}
			else
			{
				maxReadTrim = 0;
			}
		}
		if (maxReadTrim == 0) break;
	}
	if (maxReadTrim < trimAmount) return false;
	return true;
}

void trimEndRecursive(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> pos, const size_t trimAmount, const size_t kmerSize)
{
	if (resolvableGraph.precalcedUnitigLengths.size() > pos.first) resolvableGraph.precalcedUnitigLengths[pos.first] = 0;
	assert(trimAmount > 0);
	assert(resolvableGraph.unitigs[pos.first].size() > trimAmount);
	for (auto edge : resolvableGraph.edges[pos])
	{
		assert(resolvableGraph.overlaps.at(canon(pos, edge)) >= trimAmount);
		resolvableGraph.overlaps[canon(pos, edge)] -= trimAmount;
	}
	for (auto edge : resolvableGraph.edges[reverse(pos)])
	{
		size_t overlap = resolvableGraph.overlaps.at(canon(reverse(pos), edge));
		if (overlap >= resolvableGraph.unitigs[pos.first].size() - trimAmount)
		{
			assert(overlap < resolvableGraph.unitigs[pos.first].size());
			size_t trimThere = overlap - (resolvableGraph.unitigs[pos.first].size() - trimAmount) + 1;
			assert(trimThere <= trimAmount);
			assert(trimThere < resolvableGraph.unitigs[edge.first].size());
			assert(trimThere > 0);
			trimEndRecursive(resolvableGraph, readPaths, reverse(edge), trimThere, kmerSize);
			assert(resolvableGraph.overlaps.at(canon(reverse(pos), edge)) == resolvableGraph.unitigs[pos.first].size() - trimAmount - 1);
		}
		assert(resolvableGraph.overlaps.at(canon(reverse(pos), edge)) < resolvableGraph.unitigs[pos.first].size() - trimAmount);
	}
	if (pos.second)
	{
		assert(resolvableGraph.unitigRightClipBp[pos.first] < kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[pos.first][resolvableGraph.unitigs[pos.first].size()-2], resolvableGraph.unitigs[pos.first].back()));
		resolvableGraph.unitigs[pos.first].erase(resolvableGraph.unitigs[pos.first].end() - trimAmount, resolvableGraph.unitigs[pos.first].end());
		resolvableGraph.unitigRightClipBp[pos.first] = 0;
	}
	else
	{
		assert(resolvableGraph.unitigLeftClipBp[pos.first] < kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[pos.first][0], resolvableGraph.unitigs[pos.first][1]));
		resolvableGraph.unitigs[pos.first].erase(resolvableGraph.unitigs[pos.first].begin(), resolvableGraph.unitigs[pos.first].begin() + trimAmount);
		resolvableGraph.unitigLeftClipBp[pos.first] = 0;
	}
	for (const size_t i : resolvableGraph.readsCrossingNode[pos.first])
	{
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j].first != pos.first) continue;
			if (j == 0 && readPaths[i].path[j] == reverse(pos))
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					assert(readPaths[i].reads[k].leftClip >= trimAmount);
					readPaths[i].reads[k].leftClip -= trimAmount;
				}
			}
			else if (j == readPaths[i].path.size()-1 && readPaths[i].path[j] == pos)
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					assert(readPaths[i].reads[k].rightClip >= trimAmount);
					readPaths[i].reads[k].rightClip -= trimAmount;
				}
			}
			else
			{
				assert(false);
			}
		}
	}
}

void maybeTrim(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t kmerSize, std::pair<size_t, bool> pos, size_t maxTrim)
{
	if (maxTrim == 0) return;
	if (resolvableGraph.edges[pos].size() > 0) return;
	assertPrintReads(!resolvableGraph.unitigRemoved[pos.first], resolvableGraph, readPaths, pos.first);
	size_t maxReadTrim = resolvableGraph.unitigs[pos.first].size();
	assertPrintReads(resolvableGraph.readsCrossingNode[pos.first].size() > 0, resolvableGraph, readPaths, pos.first);
	for (const size_t i : resolvableGraph.readsCrossingNode[pos.first])
	{
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j].first != pos.first) continue;
			if (j == 0 && readPaths[i].path[j] == reverse(pos))
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					maxReadTrim = std::min(maxReadTrim, readPaths[i].reads[k].leftClip);
				}
			}
			else if (j == readPaths[i].path.size()-1 && readPaths[i].path[j] == pos)
			{
				for (size_t k = 0; k < readPaths[i].reads.size(); k++)
				{
					maxReadTrim = std::min(maxReadTrim, readPaths[i].reads[k].rightClip);
				}
			}
			else
			{
				maxReadTrim = 0;
			}
		}
		if (maxReadTrim == 0) break;
	}
	if (maxReadTrim == 0) return;
	assertPrintReads(maxReadTrim <= maxTrim, resolvableGraph, readPaths, pos.first);
	if (!canTrimRecursive(resolvableGraph, readPaths, pos, maxReadTrim))
	{
		std::cout << "removed trim problem node " << pos.first << std::endl;
		std::vector<std::pair<std::pair<size_t, bool>, size_t>> alsoTrim;
		for (auto edge : resolvableGraph.edges[reverse(pos)])
		{
			assertPrintReads(resolvableGraph.edges[reverse(edge)].count(pos) == 1, resolvableGraph, readPaths, pos.first);
			assertPrintReads(resolvableGraph.readsCrossingNode[edge.first].size() > 0, resolvableGraph, readPaths, pos.first);
			size_t overlap = resolvableGraph.overlaps.at(canon(reverse(pos), edge));
			if (overlap == 0) continue;
			alsoTrim.emplace_back(reverse(edge), overlap);
		}
		removeNode(resolvableGraph, readPaths, pos.first);
		for (auto pair : alsoTrim)
		{
			if (resolvableGraph.unitigRemoved[pair.first.first]) continue;
			maybeTrim(resolvableGraph, readPaths, kmerSize, pair.first, pair.second);
		}
		return;
	}
	trimEndRecursive(resolvableGraph, readPaths, pos, maxReadTrim, kmerSize);
}

void resolveRound(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool guesswork)
{
	checkValidity(resolvableGraph, readPaths, kmerSize);
	std::priority_queue<size_t, std::vector<size_t>, UnitigLengthComparer> queue { UnitigLengthComparer { resolvableGraph } };
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		queue.emplace(i);
	}
	size_t lastTopSize = 0;
	size_t nodesRemoved = 0;
	while (queue.size() > 0)
	{
		size_t topSize = resolvableGraph.unitigLength(queue.top());
		if (topSize >= maxResolveLength) break;
		// assert(topSize >= lastTopSize);
		lastTopSize = topSize;
		phmap::flat_hash_set<size_t> resolvables;
		while (queue.size() > 0 && resolvableGraph.unitigLength(queue.top()) == topSize)
		{
			if (!resolvableGraph.unitigRemoved[queue.top()])
			{
				if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 2 || resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 2)
				{
					resolvables.emplace(queue.top());
				}
				else if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 1 && resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 1)
				{
					assert(resolvableGraph.edges[std::make_pair(queue.top(), true)].size() == 1);
					assert(resolvableGraph.edges[std::make_pair(queue.top(), false)].size() == 1);
					if (resolvableGraph.unitigLength(resolvableGraph.edges[std::make_pair(queue.top(), true)].begin()->first) == topSize || resolvableGraph.unitigLength(resolvableGraph.edges[std::make_pair(queue.top(), false)].begin()->first) == topSize)
					{
						if (resolvableGraph.edges[std::make_pair(queue.top(), true)].begin()->first != queue.top())
						{
							resolvables.emplace(queue.top());
						}
					}
				}
			}
			queue.pop();
		}
		if (resolvables.size() == 0) continue;
		assert(resolvables.size() > 0);
		size_t oldSize = resolvableGraph.unitigs.size();
		checkValidity(resolvableGraph, readPaths, kmerSize);
		std::cerr << "try resolve k=" << topSize;
		auto resolutionResult = resolve(resolvableGraph, hashlist, kmerSize, readPaths, resolvables, minCoverage, topSize < maxUnconditionalResolveLength, guesswork);
		size_t newSize = resolvableGraph.unitigs.size();
		std::cerr << ", replaced " << resolutionResult.nodesResolved << " nodes with " << resolutionResult.nodesAdded << " nodes";
		nodesRemoved += resolutionResult.nodesResolved;
		if (resolutionResult.nodesResolved == 0)
		{
			assert(resolutionResult.maybeTrimmable.size() == 0);
			std::cerr << std::endl;
			continue;
		}
		size_t unitigified = 0;
		size_t unitigifiedTo = 0;
		std::vector<size_t> maybeUnitigifiable { resolutionResult.maybeUnitigifiable.begin(), resolutionResult.maybeUnitigifiable.end() };
		std::sort(maybeUnitigifiable.begin(), maybeUnitigifiable.end());
		phmap::flat_hash_map<size_t, std::pair<size_t, bool>> belongsToUnitig;
		{
			auto unitigifiedHere = unitigifySet(resolvableGraph, readPaths, maybeUnitigifiable);
			for (size_t i = 0; i < unitigifiedHere.size(); i++)
			{
				assert(unitigifiedHere[i].first.size() > 0);
				if (unitigifiedHere[i].first.size() > 1)
				{
					for (auto unitig : unitigifiedHere[i].first)
					{
						belongsToUnitig[unitig.first] = std::make_pair(unitigifiedHere[i].second, unitig.second);
					}
					nodesRemoved += unitigifiedHere[i].first.size();
					unitigified += unitigifiedHere[i].first.size();
					unitigifiedTo += 1;
				}
				else
				{
					assert(unitigifiedHere[i].first[0] == std::make_pair(unitigifiedHere[i].second, true));
					assert(!resolvableGraph.unitigRemoved[unitigifiedHere[i].second]);
					queue.emplace(unitigifiedHere[i].second);
				}
			}
		}
		if (unitigified > 0)
		{
			std::cerr << ", unitigified " << unitigified << " nodes to " << unitigifiedTo << " nodes";
		}
		std::cerr << std::endl;
		if (resolutionResult.maybeTrimmable.size() > 0)
		{
			std::vector<std::pair<std::pair<size_t, bool>, size_t>> maybeTrimmable;
			for (auto pair : resolutionResult.maybeTrimmable)
			{
				std::pair<size_t, bool> trimIndex = pair.first;
				if (belongsToUnitig.count(trimIndex.first) == 1)
				{
					bool fw = trimIndex.second;
					trimIndex = belongsToUnitig.at(trimIndex.first);
					if (!fw) trimIndex.second = !trimIndex.second;
				}
				if (resolvableGraph.unitigRemoved[trimIndex.first]) continue;
				maybeTrimmable.emplace_back(trimIndex, pair.second);
			}
			std::sort(maybeTrimmable.begin(), maybeTrimmable.end(), [](const std::pair<std::pair<size_t, bool>, size_t>& left, const std::pair<std::pair<size_t, bool>, size_t>& right)
			{
				if (left.first < right.first) return true;
				if (left.first > right.first) return false;
				assert(left.first == right.first);
				if (left.second > right.second) return true;
				if (left.second < right.second) return false;
				assert(left.second == right.second);
				return false;
			});
			for (auto pair : maybeTrimmable)
			{
				assertPrintReads(!resolvableGraph.unitigRemoved[pair.first.first], resolvableGraph, readPaths, pair.first.first);
			}
			for (auto pair : maybeTrimmable)
			{
				if (resolvableGraph.unitigRemoved[pair.first.first]) continue;
				maybeTrim(resolvableGraph, readPaths, kmerSize, pair.first, pair.second);
			}
		}
		auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, kmerSize, resolutionResult.maybeUnitigifiable);
		resolutionResult.maybeUnitigifiable.insert(removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end());
		auto removed2 = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 2, 5, 10000, kmerSize, resolutionResult.maybeUnitigifiable);
		nodesRemoved += removed.nodesRemoved + removed2.nodesRemoved;
		if (removed.nodesRemoved + removed2.nodesRemoved > 0)
		{
			std::cerr << "removed " << removed.nodesRemoved + removed2.nodesRemoved << " tips" << std::endl;
			removed.maybeUnitigifiable.insert(removed2.maybeUnitigifiable.begin(), removed2.maybeUnitigifiable.end());
			std::vector<size_t> maybeUnitigifiable2 { removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end() };
			auto unitigifiedHere = unitigifySet(resolvableGraph, readPaths, maybeUnitigifiable2);
			for (size_t i = 0; i < unitigifiedHere.size(); i++)
			{
				assert(unitigifiedHere[i].first.size() > 0);
				if (unitigifiedHere[i].first.size() >= 2)
				{
					nodesRemoved += unitigifiedHere[i].first.size();
				}
				else
				{
					assert(unitigifiedHere[i].first[0] == std::make_pair(unitigifiedHere[i].second, true));
					queue.emplace(unitigifiedHere[i].second);
				}
			}
		}
		checkValidity(resolvableGraph, readPaths, kmerSize);
		for (size_t i = oldSize; i < resolvableGraph.unitigs.size(); i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			queue.emplace(i);
		}
		if (nodesRemoved > resolvableGraph.unitigs.size() / 2 + 50)
		{
			std::vector<size_t> queueNodes;
			while (queue.size() > 0)
			{
				auto top = queue.top();
				queue.pop();
				if (resolvableGraph.unitigRemoved[top]) continue;
				queueNodes.emplace_back(top);
			}
			size_t oldSize = resolvableGraph.unitigs.size();
			compact(resolvableGraph, readPaths, queueNodes);
			nodesRemoved = 0;
			size_t newSize = resolvableGraph.unitigs.size();
			std::cerr << "compacted from " << oldSize << " to " << newSize << std::endl;
			checkValidity(resolvableGraph, readPaths, kmerSize);
			for (auto node : queueNodes)
			{
				assert(node < resolvableGraph.unitigs.size());
				assert(!resolvableGraph.unitigRemoved[node]);
				queue.emplace(node);
			}
		}
	}
}

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, std::vector<ReadPath>& rawReadPaths, const ReadpartIterator& partIterator, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool keepGaps, const bool guesswork)
{
	auto resolvableGraph = getUnitigs(initial, minCoverage, hashlist, kmerSize, keepGaps);
	std::cerr << rawReadPaths.size() << " raw read paths" << std::endl;
	// todo maybe fix? or does it matter?
	cutRemovedEdgesFromPaths(resolvableGraph, rawReadPaths);
	std::sort(rawReadPaths.begin(), rawReadPaths.end(), [](const ReadPath& left, const ReadPath& right) { return left.path < right.path; });
	std::vector<PathGroup> readPaths;
	for (size_t i = 0; i < rawReadPaths.size(); i++)
	{
		if (rawReadPaths[i].path.size() == 0) continue;
		if (readPaths.size() == 0 || readPaths.back().path != rawReadPaths[i].path)
		{
			readPaths.emplace_back();
			readPaths.back().path = rawReadPaths[i].path;
			for (auto pos : rawReadPaths[i].path)
			{
				resolvableGraph.readsCrossingNode[pos.first].insert(readPaths.size()-1);
			}
		}
		assert(readPaths.size() > 0);
		assert(readPaths.back().path.size() > 0);
		assert(readPaths.back().path == rawReadPaths[i].path);
		resolvableGraph.readNames.emplace_back(rawReadPaths[i].readName);
		readPaths.back().reads.emplace_back();
		readPaths.back().reads.back().readNameIndex = resolvableGraph.readNames.size()-1;
		readPaths.back().reads.back().readPoses = rawReadPaths[i].readPoses;
		readPaths.back().reads.back().expandedReadPosStart = rawReadPaths[i].expandedReadPosStart;
		readPaths.back().reads.back().expandedReadPosEnd = rawReadPaths[i].expandedReadPosEnd;
		readPaths.back().reads.back().leftClip = rawReadPaths[i].leftClip;
		readPaths.back().reads.back().rightClip = rawReadPaths[i].rightClip;
		readPaths.back().reads.back().readLength = rawReadPaths[i].readLength;
		readPaths.back().reads.back().readLengthHPC = rawReadPaths[i].readLengthHPC;
	}
	std::cerr << rawReadPaths.size() << " raw read paths" << std::endl;
	std::cerr << readPaths.size() << " read paths" << std::endl;
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		assert(resolvableGraph.readsCrossingNode[i].size() >= 1);
	}
	unitigifyAll(resolvableGraph, readPaths);
	auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, kmerSize, phmap::flat_hash_set<size_t> {});
	if (removed.nodesRemoved > 0)
	{
		std::cerr << "removed " << removed.nodesRemoved << " tips" << std::endl;
		unitigifyAll(resolvableGraph, readPaths);
	}
	resolveRound(resolvableGraph, readPaths, hashlist, minCoverage, kmerSize, maxResolveLength, maxUnconditionalResolveLength, guesswork);
	resolveRound(resolvableGraph, readPaths, hashlist, 1, kmerSize, maxResolveLength, maxUnconditionalResolveLength, guesswork);
	checkValidity(resolvableGraph, readPaths, kmerSize);
	return resolvableToUnitigs(resolvableGraph, readPaths);
}
