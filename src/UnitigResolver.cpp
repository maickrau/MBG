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
#include "Node.h"

#define assertPrintReads(expression, graph, paths, node) {if (!(expression)) {printReads(graph, paths, node);} assert(expression);}

struct ReadPathInfo
{
	std::vector<uint32_t> readPoses;
	size_t expandedReadPosStart;
	size_t expandedReadPosEnd;
	size_t readLength;
	size_t readLengthHPC;
};

struct ResolveTriplet
{
	ResolveTriplet() = default;
	ResolveTriplet(std::pair<size_t, bool> left, std::pair<size_t, bool> right, size_t coverage) :
		left(left),
		right(right),
		coverage(coverage)
	{
	}
	std::pair<size_t, bool> left;
	std::pair<size_t, bool> right;
	size_t coverage;
};

class PathGroup
{
public:
	class Read
	{
	public:
		size_t readNameIndex;
		size_t readInfoIndex;
		size_t readPosZeroOffset;
		size_t readPosStartIndex;
		size_t readPosEndIndex;
		size_t leftClip;
		size_t rightClip;
	};
	std::vector<std::pair<size_t, bool>> path;
	std::vector<Read> reads;
};

class ReadCrosserIterator
{
public:
	ReadCrosserIterator(size_t index, std::vector<std::pair<uint32_t, uint32_t>>& vec, const std::vector<PathGroup>& paths) :
		index(index),
		vec(vec),
		paths(paths)
	{
	}
	bool operator==(const ReadCrosserIterator& other) const
	{
		return index == other.index;
	}
	bool operator!=(const ReadCrosserIterator& other) const
	{
		return !(*this == other);
	}
	std::pair<uint32_t, uint32_t> operator*() const
	{
		return vec[index];
	}
	ReadCrosserIterator& operator++()
	{
		while (true)
		{
			if (index == 0)
			{
				index = std::numeric_limits<size_t>::max();
				break;
			}
			index -= 1;
			if (paths[vec[index].first].path.size() == 0 || paths[vec[index].first].reads.size() == 0)
			{
				std::swap(vec[index], vec.back());
				vec.pop_back();
				continue;
			}
			else
			{
				break;
			}
		}
		return *this;
	}
private:
	size_t index;
	std::vector<std::pair<uint32_t, uint32_t>>& vec;
	const std::vector<PathGroup>& paths;
};

class ReadCrosserIteratorHelper
{
public:
	ReadCrosserIteratorHelper(std::vector<std::pair<uint32_t, uint32_t>>& vec, const std::vector<PathGroup>& paths) :
		vec(vec),
		paths(paths)
	{
	}
	ReadCrosserIterator begin()
	{
		return ReadCrosserIterator { vec.size()-1, vec, paths };
	}
	ReadCrosserIterator end()
	{
		return ReadCrosserIterator { std::numeric_limits<size_t>::max(), vec, paths };
	}
private:
	std::vector<std::pair<uint32_t, uint32_t>>& vec;
	const std::vector<PathGroup>& paths;
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
	mutable std::vector<std::vector<std::pair<uint32_t, uint32_t>>> readsCrossingNode;
	std::vector<size_t> everTippable;
	size_t lastTippableChecked;
	mutable std::vector<size_t> precalcedUnitigLengths;
	mutable std::vector<double> precalcedUnitigCoverages;
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
	double getCoverage(const std::vector<PathGroup>& readPaths, size_t unitig) const
	{
		if (unitig >= precalcedUnitigCoverages.size())
		{
			assert(precalcedUnitigCoverages.size() < unitigs.size());
			precalcedUnitigCoverages.resize(unitigs.size(), 0);
		}
		assert(unitig < precalcedUnitigCoverages.size());
		if (precalcedUnitigCoverages[unitig] == 0)
		{
			precalcedUnitigCoverages[unitig] = calculateCoverage(readPaths, unitig);
		}
		assert(precalcedUnitigCoverages[unitig] != 0);
		return precalcedUnitigCoverages[unitig];
	}
	const HashList& hashlist;
	std::vector<std::string> readNames;
	double averageCoverage;
	double calculateCoverage(const std::vector<PathGroup>& readPaths, size_t unitig) const
	{
		double result = 0;
		for (const std::pair<uint32_t, uint32_t> pospair : iterateCrossingReads(unitig, readPaths))
		{
			const size_t i = pospair.first;
			const size_t j = pospair.second;
			assert(readPaths[i].path.size() > 0);
			assert(readPaths[i].reads.size() > 0);
			if (readPaths[i].path.size() == 1)
			{
				assert(readPaths[i].path[0].first == unitig);
				for (const auto& read : readPaths[i].reads)
				{
					assert(read.leftClip + read.rightClip < unitigs[readPaths[i].path[0].first].size());
					result += (double)(unitigs[unitig].size() - read.leftClip - read.rightClip) / (double)(unitigs[unitig].size());
				}
			}
			else
			{
				if (j == 0)
				{
					for (const auto& read : readPaths[i].reads)
					{
						assert(read.leftClip < unitigs[unitig].size());
						result += (double)(unitigs[unitig].size() - read.leftClip) / (double)(unitigs[unitig].size());
					}
				}
				if (j > 0 && j < readPaths[i].path.size()-1) result += readPaths[i].reads.size();
				if (j == readPaths[i].path.size()-1)
				{
					for (const auto& read : readPaths[i].reads)
					{
						assert(read.rightClip < unitigs[unitig].size());
						result += (double)(unitigs[unitig].size() - read.rightClip) / (double)(unitigs[unitig].size());
					}
				}
			}
		}
		return result;
	}
	ReadCrosserIteratorHelper iterateCrossingReads(size_t node, const std::vector<PathGroup>& paths) const
	{
		return ReadCrosserIteratorHelper { readsCrossingNode[node], paths };
	}
	const size_t kmerSize;
private:
};

void printReads(const ResolvableUnitigGraph& graph, const std::vector<PathGroup>& paths, size_t node)
{
	std::unordered_set<size_t> closeReads;
	std::unordered_set<size_t> semiCloseReads;
	std::cerr << "around center node " << node << std::endl;
	for (const std::pair<uint32_t, uint32_t> pospair : graph.iterateCrossingReads(node, paths))
	{
		for (const auto read : paths[pospair.first].reads)
		{
			if (closeReads.count(read.readNameIndex) == 1) continue;
			std::cerr << "Read close to assertion: " << graph.readNames[read.readNameIndex] << std::endl;
			closeReads.insert(read.readNameIndex);
		}
	}
	for (auto edge : graph.edges[std::make_pair(node, true)])
	{
		std::cerr << "around edge node " << edge.first << std::endl;
		for (const std::pair<uint32_t, uint32_t> pospair : graph.iterateCrossingReads(edge.first, paths))
		{
			for (const auto read : paths[pospair.first].reads)
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
		for (const std::pair<uint32_t, uint32_t> pospair : graph.iterateCrossingReads(edge.first, paths))
		{
			for (auto read : paths[pospair.first].reads)
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
		decltype(resolvableGraph.readsCrossingNode) tmp;
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
		for (size_t i = 0; i < resolvableGraph.precalcedUnitigLengths.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			newPrecalcedUnitigLengths[newIndex] = resolvableGraph.precalcedUnitigLengths[i];
		}
		std::swap(newPrecalcedUnitigLengths, resolvableGraph.precalcedUnitigLengths);
	}
	{
		std::vector<double> newPrecalcedUnitigCoverages;
		newPrecalcedUnitigCoverages.resize(newSize, 0);
		for (size_t i = 0; i < resolvableGraph.precalcedUnitigCoverages.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			newPrecalcedUnitigCoverages[newIndex] = resolvableGraph.precalcedUnitigCoverages[i];
		}
		std::swap(newPrecalcedUnitigCoverages, resolvableGraph.precalcedUnitigCoverages);
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
		for (size_t j = 0; j < paths[i].path.size(); j++)
		{
			resolvableGraph.readsCrossingNode[paths[i].path[j].first].emplace_back(i, j);
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
		for (auto pair : initial.edgeCov.getValues(fw))
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
		for (auto pair : initial.edgeCov.getValues(bw))
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
			for (auto pair : initial.edgeCov.getValues(tip))
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
	assert(longkmerDivisor > 0 || longkmerSum == 0);
	if (longkmerDivisor > 0)
	{
		result.averageCoverage = longkmerSum / longkmerDivisor;
	}
	else
	{
		assert(allkmerDivisor > 0 || allkmerSum == 0);
		result.averageCoverage = 0;
		if (allkmerDivisor > 0)
		{
			result.averageCoverage = allkmerSum / allkmerDivisor;
		}
	}
	std::cerr << "estimated average coverage " << result.averageCoverage << std::endl;
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

std::pair<UnitigGraph, std::vector<ReadPath>> resolvableToUnitigs(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const std::vector<ReadPathInfo>& readInfos)
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
					assert(read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex) == pathHashCount);
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
				assert(read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex) == pathHashCount);
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
			assert(read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex) == pathHashCount);
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
			result.edges.addEdge(newfw, newEdge);
			result.edges.addEdge(reverse(newEdge), reverse(newfw));
			result.setEdgeCoverage(newfw, newEdge, 0);
			result.setEdgeOverlap(newfw, newEdge, resolvableGraph.overlaps.at(canon(fw, edge)));
		}
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newbw { newIndex.getRank(i), false };
		assert(newbw.first < newSize);
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			assert(newEdge.first < newSize);
			result.edges.addEdge(newbw, newEdge);
			result.edges.addEdge(reverse(newEdge), reverse(newbw));
			result.setEdgeCoverage(newbw, newEdge, 0);
			result.setEdgeOverlap(newbw, newEdge, resolvableGraph.overlaps.at(canon(bw, edge)));
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
			result.setEdgeCoverage(fixPath[j-1], fixPath[j], result.edgeCoverage(fixPath[j-1], fixPath[j]) + path.reads.size());
		}
		for (const auto& read : path.reads)
		{
			resultReads.emplace_back();
			resultReads.back().path = std::vector<Node> { fixPath.begin(), fixPath.end() };
			resultReads.back().readName = resolvableGraph.readNames[read.readNameIndex];
			resultReads.back().readPoses = { readInfos[read.readInfoIndex].readPoses.begin() + read.readPosStartIndex, readInfos[read.readInfoIndex].readPoses.begin() + read.readPosEndIndex };
			resultReads.back().expandedReadPosStart = readInfos[read.readInfoIndex].expandedReadPosStart;
			resultReads.back().expandedReadPosEnd = readInfos[read.readInfoIndex].expandedReadPosEnd;
			resultReads.back().leftClip = read.leftClip;
			resultReads.back().rightClip = read.rightClip;
			resultReads.back().readLength = readInfos[read.readInfoIndex].readLength;
			resultReads.back().readLengthHPC = readInfos[read.readInfoIndex].readLengthHPC;
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
	for (size_t j = 0; j < readPaths[i].path.size(); j++)
	{
		// resolvableGraph.readsCrossingNode[readPaths[i].path[j].first].erase(std::make_pair(i, j));
		if (resolvableGraph.precalcedUnitigCoverages.size() > readPaths[i].path[j].first) resolvableGraph.precalcedUnitigCoverages[readPaths[i].path[j].first] = 0;
	}
	readPaths[i].path.clear();
	readPaths[i].reads.clear();
}

void addPath(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, PathGroup&& newPath)
{
	assert(newPath.path.size() > 0);
	assert(newPath.reads.size() > 0);
	size_t leftStart = 0;
	size_t rightEnd = 0;
	if (newPath.path.size() >= 2)
	{
		leftStart = resolvableGraph.unitigs[newPath.path[0].first].size() - resolvableGraph.overlaps.at(canon(newPath.path[0], newPath.path[1]));
		rightEnd = resolvableGraph.unitigs[newPath.path.back().first].size() - resolvableGraph.overlaps.at(canon(newPath.path[newPath.path.size()-2], newPath.path[newPath.path.size()-1]));
		if (newPath.path[1].second)
		{
			if (resolvableGraph.unitigLeftClipBp[newPath.path[1].first] > 0)
			{
				leftStart += 1;
			}
		}
		else
		{
			if (resolvableGraph.unitigRightClipBp[newPath.path[1].first] > 0)
			{
				leftStart += 1;
			}
		}
		if (newPath.path[newPath.path.size()-2].second)
		{
			if (resolvableGraph.unitigRightClipBp[newPath.path[newPath.path.size()-2].first] > 0)
			{
				rightEnd += 1;
			}
		}
		else
		{
			if (resolvableGraph.unitigLeftClipBp[newPath.path[newPath.path.size()-2].first] > 0)
			{
				rightEnd += 1;
			}
		}
	}
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
			assert(read.leftClip < leftStart);
			assert(read.rightClip < rightEnd);
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
	for (size_t j = 0; j < newPath.path.size(); j++)
	{
		resolvableGraph.readsCrossingNode[newPath.path[j].first].emplace_back(newIndex, j);
		if (resolvableGraph.precalcedUnitigCoverages.size() > newPath.path[j].first) resolvableGraph.precalcedUnitigCoverages[newPath.path[j].first] = 0;
	}
	readPaths.emplace_back(std::move(newPath));
}

void addPathButFirstMaybeTrim(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, PathGroup&& newPath)
{
	assert(newPath.path.size() > 0);
	assert(newPath.reads.size() > 0);
	if (newPath.path.size() == 1)
	{
		addPath(resolvableGraph, readPaths, std::move(newPath));
		return;
	}
	size_t leftTrim = resolvableGraph.unitigs[newPath.path[0].first].size() - resolvableGraph.overlaps.at(canon(newPath.path[0], newPath.path[1]));
	size_t rightTrim = resolvableGraph.unitigs[newPath.path.back().first].size() - resolvableGraph.overlaps.at(canon(newPath.path[newPath.path.size()-2], newPath.path[newPath.path.size()-1]));
	size_t leftCompare = leftTrim;
	size_t rightCompare = rightTrim;
	if (newPath.path[1].second)
	{
		if (resolvableGraph.unitigLeftClipBp[newPath.path[1].first] > 0)
		{
			leftCompare += 1;
		}
	}
	else
	{
		if (resolvableGraph.unitigRightClipBp[newPath.path[1].first] > 0)
		{
			leftCompare += 1;
		}
	}
	if (newPath.path[newPath.path.size()-2].second)
	{
		if (resolvableGraph.unitigRightClipBp[newPath.path[newPath.path.size()-2].first] > 0)
		{
			rightCompare += 1;
		}
	}
	else
	{
		if (resolvableGraph.unitigLeftClipBp[newPath.path[newPath.path.size()-2].first] > 0)
		{
			rightCompare += 1;
		}
	}
	assert(leftCompare > 0);
	assert(rightCompare > 0);
	PathGroup leftTrimPath;
	PathGroup rightTrimPath;
	PathGroup doubleTrimPath;
	for (size_t i = newPath.reads.size()-1; i < newPath.reads.size(); i--)
	{
		const auto& read = newPath.reads[i];
		assert(read.leftClip < leftCompare || read.rightClip < rightCompare || newPath.path.size() >= 3);
		if (read.leftClip >= leftCompare && read.rightClip >= rightCompare)
		{
			if (doubleTrimPath.path.size() == 0)
			{
				assert(newPath.path.size() >= 3);
				doubleTrimPath.path.insert(doubleTrimPath.path.end(), newPath.path.begin()+1, newPath.path.end()-1);
				assert(doubleTrimPath.path.size() > 0);
			}
			doubleTrimPath.reads.push_back(read);
			assert(doubleTrimPath.reads.back().leftClip >= leftTrim);
			doubleTrimPath.reads.back().leftClip -= leftTrim;
			assert(doubleTrimPath.reads.back().rightClip >= rightTrim);
			doubleTrimPath.reads.back().rightClip -= rightTrim;
			std::swap(newPath.reads[i], newPath.reads.back());
			newPath.reads.pop_back();
			continue;
		}
		if (read.leftClip >= leftCompare)
		{
			if (leftTrimPath.path.size() == 0)
			{
				leftTrimPath.path.insert(leftTrimPath.path.end(), newPath.path.begin()+1, newPath.path.end());
				assert(leftTrimPath.path.size() > 0);
			}
			leftTrimPath.reads.push_back(read);
			assert(leftTrimPath.reads.back().leftClip >= leftTrim);
			leftTrimPath.reads.back().leftClip -= leftTrim;
			std::swap(newPath.reads[i], newPath.reads.back());
			newPath.reads.pop_back();
			continue;
		}
		if (read.rightClip >= rightCompare)
		{
			if (rightTrimPath.path.size() == 0)
			{
				rightTrimPath.path.insert(rightTrimPath.path.end(), newPath.path.begin(), newPath.path.end()-1);
				assert(rightTrimPath.path.size() > 0);
			}
			rightTrimPath.reads.push_back(read);
			assert(rightTrimPath.reads.back().rightClip >= rightTrim);
			rightTrimPath.reads.back().rightClip -= rightTrim;
			std::swap(newPath.reads[i], newPath.reads.back());
			newPath.reads.pop_back();
			continue;
		}
	}
	assert((leftTrimPath.path.size() == 0) == (leftTrimPath.reads.size() == 0));
	assert((rightTrimPath.path.size() == 0) == (rightTrimPath.reads.size() == 0));
	assert((doubleTrimPath.path.size() == 0) == (doubleTrimPath.reads.size() == 0));
	if (newPath.reads.size() > 0) addPath(resolvableGraph, readPaths, std::move(newPath));
	if (leftTrimPath.reads.size() > 0) addPath(resolvableGraph, readPaths, std::move(leftTrimPath));
	if (rightTrimPath.reads.size() > 0) addPath(resolvableGraph, readPaths, std::move(rightTrimPath));
	if (doubleTrimPath.reads.size() > 0) addPath(resolvableGraph, readPaths, std::move(doubleTrimPath));
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
		assert((resolvableGraph.unitigLength(newUnitigs[i].first[0].first) < resolvableGraph.unitigLength(newUnitigs[i].second)) || (resolvableGraph.unitigLength(newUnitigs[i].first.back().first) < resolvableGraph.unitigLength(newUnitigs[i].second)));
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
				for (auto pospair : resolvableGraph.iterateCrossingReads(pos.first, readPaths)) relevantReads.insert(pospair.first);
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
		addPathButFirstMaybeTrim(resolvableGraph, readPaths, std::move(newPath));
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
			pos += graph.unitigs[readPaths[i].path[j].id()].size();
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

void createFakeEdgeNode(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes, const phmap::flat_hash_set<size_t>& resolvables, const phmap::flat_hash_set<size_t>& unresolvables, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
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
	resolvableGraph.unitigRightClipBp.push_back(rightClipBp);
	resolvableGraph.unitigLeftClipBp.push_back(leftClipBp);
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	assert(resolvables.count(to.first) == 0 || unresolvables.count(to.first) == 1);
	resolvableGraph.edges[std::make_pair(newIndex, true)].emplace(to);
	resolvableGraph.edges[reverse(to)].emplace(std::make_pair(newIndex, false));
	size_t overlap = resolvableGraph.overlaps.at(canon(from, to));
	assert(resolvableGraph.unitigs[newIndex].size() >= overlap);
	assert(resolvableGraph.unitigs[to.first].size() >= overlap);
	resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = overlap;
	assert(resolvableGraph.getBpOverlap(std::make_pair(newIndex, true), to) < resolvableGraph.unitigLength(to.first));
	assert(resolvableGraph.unitigLength(newIndex) == resolvableGraph.unitigLength(from.first));
	if (resolvableGraph.unitigs[newIndex].size() >= 2)
	{
		assert(resolvableGraph.unitigLeftClipBp[newIndex] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][0], resolvableGraph.unitigs[newIndex][1]));
		assert(resolvableGraph.unitigRightClipBp[newIndex] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-2], resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-1]));
	}
}

void createEdgeNode(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes, const phmap::flat_hash_set<size_t>& resolvables, const phmap::flat_hash_set<size_t>& unresolvables, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
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
			rightClipBp = resolvableGraph.kmerSize - hashlist.getOverlap(resolvableGraph.unitigs.back().back(), add[start]) - 1;
			assert(rightClipBp < resolvableGraph.kmerSize);
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
			assert(rightClipBp < resolvableGraph.kmerSize);
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
		assert(resolvableGraph.unitigLeftClipBp[newIndex] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][0], resolvableGraph.unitigs[newIndex][1]));
		assert(resolvableGraph.unitigRightClipBp[newIndex] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-2], resolvableGraph.unitigs[newIndex][resolvableGraph.unitigs[newIndex].size()-1]));
	}
}

size_t getTrimAmountToCheck(const ResolvableUnitigGraph& resolvableGraph, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	if (from.second && resolvableGraph.unitigRightClipBp[from.first] == 1)
	{
		return resolvableGraph.unitigs[to.first].size() - resolvableGraph.overlaps.at(canon(from, to));
	}
	if (!from.second && resolvableGraph.unitigLeftClipBp[from.first] == 1)
	{
		return resolvableGraph.unitigs[to.first].size() - resolvableGraph.overlaps.at(canon(from, to));
	}
	return 0;
}

std::vector<ResolveTriplet> getRawTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool partTriplets)
{
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> tripletCoverage;
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(node, readPaths))
	{
		const size_t i = pospair.first;
		const size_t j = pospair.second;
		if (readPaths[i].path.size() == 0) continue;
		if (readPaths[i].path.size() == 1) continue;
		if (readPaths[i].reads.size() == 0) continue;
		assert(readPaths[i].path[j].first == node);
		if (j == 0 && resolvableGraph.edges[reverse(readPaths[i].path[0])].size() == 0)
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
			size_t shouldCheckTrimAmount = 0;
			if (readPaths[i].path.size() == 2)
			{
				shouldCheckTrimAmount = getTrimAmountToCheck(resolvableGraph, readPaths[i].path[0], readPaths[i].path[1]);
			}
			size_t coverageHere;
			if (shouldCheckTrimAmount == 0)
			{
				coverageHere = readPaths[i].reads.size();
			}
			else
			{
				coverageHere = 0;
				for (const auto& read : readPaths[i].reads)
				{
					if (read.rightClip < shouldCheckTrimAmount) coverageHere += 1;
				}
			}
			if (coverageHere > 0)
			{
				tripletCoverage[std::make_pair(left, right)] += coverageHere;
			}
		}
		if (j == readPaths[i].path.size()-1 && resolvableGraph.edges[readPaths[i].path.back()].size() == 0)
		{
			std::pair<size_t, bool> left { std::numeric_limits<size_t>::max(), true };
			std::pair<size_t, bool> right { std::numeric_limits<size_t>::max(), true };
			if (readPaths[i].path.back().second)
			{
				left = readPaths[i].path[readPaths[i].path.size()-2];
			}
			else
			{
				right = reverse(readPaths[i].path[readPaths[i].path.size()-2]);
			}
			size_t shouldCheckTrimAmount = 0;
			if (readPaths[i].path.size() == 2)
			{
				shouldCheckTrimAmount = getTrimAmountToCheck(resolvableGraph, reverse(readPaths[i].path[1]), reverse(readPaths[i].path[0]));
			}
			size_t coverageHere;
			if (shouldCheckTrimAmount == 0)
			{
				coverageHere = readPaths[i].reads.size();
			}
			else
			{
				coverageHere = 0;
				for (const auto& read : readPaths[i].reads)
				{
					if (read.leftClip < shouldCheckTrimAmount) coverageHere += 1;
				}
			}
			if (coverageHere > 0)
			{
				tripletCoverage[std::make_pair(left, right)] += coverageHere;
			}
		}
		if (j > 0 && j < readPaths[i].path.size()-1)
		{
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
			assert(!resolvableGraph.unitigRemoved[node]);
			assert(!resolvableGraph.unitigRemoved[left.first]);
			assert(!resolvableGraph.unitigRemoved[right.first]);
			assert(resolvableGraph.edges[left].count(std::make_pair(node, true)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, false)].count(reverse(left)) == 1);
			assert(resolvableGraph.edges[reverse(right)].count(std::make_pair(node, false)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, true)].count(right) == 1);
			size_t shouldCheckLeftTrimAmount = 0;
			size_t shouldCheckRightTrimAmount = 0;
			if (j == 1)
			{
				shouldCheckLeftTrimAmount = getTrimAmountToCheck(resolvableGraph, reverse(readPaths[i].path[1]), reverse(readPaths[i].path[0]));
			}
			if (j == readPaths[i].path.size()-2)
			{
				shouldCheckRightTrimAmount = getTrimAmountToCheck(resolvableGraph, readPaths[i].path[j], readPaths[i].path[j+1]);
			}
			size_t coverageHere;
			if (shouldCheckLeftTrimAmount == 0 && shouldCheckRightTrimAmount == 0)
			{
				coverageHere = readPaths[i].reads.size();
			}
			else
			{
				coverageHere = 0;
				for (const auto& read : readPaths[i].reads)
				{
					if ((shouldCheckLeftTrimAmount == 0 || read.leftClip < shouldCheckLeftTrimAmount) && (shouldCheckRightTrimAmount == 0 || read.rightClip < shouldCheckRightTrimAmount)) coverageHere += 1;
				}
			}
			if (coverageHere > 0)
			{
				tripletCoverage[std::make_pair(left, right)] += coverageHere;
			}
		}
	}
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> partTripletCoverage = tripletCoverage;
	if (partTriplets && resolvableGraph.edges[std::make_pair(node, true)].size() >= 1 && resolvableGraph.edges[std::make_pair(node, false)].size() >= 1)
	{
		std::unordered_map<size_t, std::pair<size_t, std::pair<size_t, size_t>>> fwReads;
		std::unordered_map<size_t, std::pair<size_t, std::pair<size_t, size_t>>> bwReads;
		for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(node, readPaths))
		{
			const size_t i = pospair.first;
			const size_t j = pospair.second;
			if (readPaths[i].path.size() < 2) continue;
			if (readPaths[i].path[0].first == node && readPaths[i].path.back().first == node) continue;
			if (j == 0)
			{
				for (const auto& read : readPaths[i].reads)
				{
					if (fwReads.count(read.readNameIndex) == 1)
					{
						fwReads[read.readNameIndex] = std::make_pair(std::numeric_limits<size_t>::max(), std::make_pair(0, 0));
					}
					else
					{
						fwReads[read.readNameIndex] = std::make_pair(i, std::make_pair(read.readPosZeroOffset, read.readPosStartIndex));
					}
				}
			}
			if (j == readPaths[i].path.size()-1)
			{
				for (const auto& read : readPaths[i].reads)
				{
					if (bwReads.count(read.readNameIndex) == 1)
					{
						bwReads[read.readNameIndex] = std::make_pair(std::numeric_limits<size_t>::max(), std::make_pair(0, 0));
					}
					else
					{
						bwReads[read.readNameIndex] = std::make_pair(i, std::make_pair(read.readPosZeroOffset, read.readPosEndIndex));
					}
				}
			}
		}
		std::unordered_set<std::pair<size_t, size_t>> pairs;
		for (auto pair : fwReads)
		{
			if (pair.second.first == std::numeric_limits<size_t>::max()) continue;
			if (bwReads.count(pair.first) == 0) continue;
			auto other = bwReads.at(pair.first);
			if (other.first == std::numeric_limits<size_t>::max()) continue;
			if (pair.second.first == other.first) continue;
			if (other.second.first > pair.second.second.first) continue;
			if (other.second.first == pair.second.second.first && other.second.second > pair.second.second.second) continue;
			if (readPaths[pair.second.first].path[0] != readPaths[other.first].path.back()) continue;
			pairs.emplace(other.first, pair.second.first);
		}
		for (auto pair : pairs)
		{
			assert(readPaths[pair.first].path.size() >= 2);
			assert(readPaths[pair.second].path.size() >= 2);
			assert(readPaths[pair.first].path.back() == readPaths[pair.second].path[0]);
			size_t minLeftClip = resolvableGraph.unitigs[node].size();
			size_t minRightClip = resolvableGraph.unitigs[node].size();
			size_t coverage = 0;
			for (const auto& read : readPaths[pair.first].reads)
			{
				minLeftClip = std::min(minLeftClip, read.rightClip);
			}
			for (const auto& read : readPaths[pair.second].reads)
			{
				minRightClip = std::min(minRightClip, read.leftClip);
			}
			std::unordered_set<size_t> rightReadNames;
			for (const auto& read : readPaths[pair.second].reads)
			{
				rightReadNames.insert(read.readNameIndex);
			}
			for (const auto& read : readPaths[pair.first].reads)
			{
				if (rightReadNames.count(read.readNameIndex) == 1)
				{
					coverage += 1;
				}
			}
			assert(coverage >= 1);
			if (minLeftClip + minRightClip < resolvableGraph.unitigs[node].size())
			{
				if (readPaths[pair.first].path.back().second)
				{
					partTripletCoverage[std::make_pair(readPaths[pair.first].path[readPaths[pair.first].path.size()-2], readPaths[pair.second].path[1])] += coverage;
				}
				else
				{
					partTripletCoverage[std::make_pair(reverse(readPaths[pair.second].path[1]), reverse(readPaths[pair.first].path[readPaths[pair.first].path.size()-2]))] += coverage;
				}
			}
		}
	}
	std::vector<ResolveTriplet> coveredTriplets;
	for (auto pair : tripletCoverage)
	{
		if (pair.second < minCoverage) continue;
		coveredTriplets.emplace_back(pair.first.first, pair.first.second, pair.second);
	}
	if (partTriplets && resolvableGraph.edges[std::make_pair(node, true)].size() >= 1 && resolvableGraph.edges[std::make_pair(node, false)].size() >= 1)
	{
		bool canAddPartTriplets = true;
		std::unordered_set<std::pair<size_t, bool>> uniqueBwMatches;
		std::unordered_set<std::pair<size_t, bool>> uniqueFwMatches;
		for (auto pair : partTripletCoverage)
		{
			if (pair.second < minCoverage) continue;
			if (uniqueBwMatches.count(pair.first.first) == 1)
			{
				canAddPartTriplets = false;
				break;
			}
			if (uniqueFwMatches.count(pair.first.second) == 1)
			{
				canAddPartTriplets = false;
				break;
			}
			if (resolvableGraph.edges[pair.first.first].size() != 1)
			{
				canAddPartTriplets = false;
				break;
			}
			if (resolvableGraph.edges[reverse(pair.first.second)].size() != 1)
			{
				canAddPartTriplets = false;
				break;
			}
			uniqueBwMatches.insert(pair.first.first);
			uniqueFwMatches.insert(pair.first.second);
		}
		if (canAddPartTriplets)
		{
			coveredTriplets.clear();
			for (auto pair : partTripletCoverage)
			{
				if (pair.second < minCoverage) continue;
				coveredTriplets.emplace_back(pair.first.first, pair.first.second, pair.second);
			}
		}
	}
	return coveredTriplets;
}

std::vector<ResolveTriplet> getReadSupportedTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional, bool guesswork)
{
	auto coveredTriplets = getRawTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, guesswork);
	if (unconditional && coveredTriplets.size() >= 2) return coveredTriplets;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredInNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		if (pair.left.first != std::numeric_limits<size_t>::max()) coveredInNeighbors.emplace(reverse(pair.left));
		if (pair.right.first != std::numeric_limits<size_t>::max()) coveredOutNeighbors.emplace(pair.right);
	}
	assert(coveredInNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, false)].size());
	assert(coveredOutNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, true)].size());
	std::vector<ResolveTriplet> empty;
	if (coveredInNeighbors.size() < resolvableGraph.edges[std::make_pair(node, false)].size()) return empty;
	if (coveredOutNeighbors.size() < resolvableGraph.edges[std::make_pair(node, true)].size()) return empty;
	return coveredTriplets;
}

size_t getAnchorSize(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	size_t result = 0;
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(from.first, readPaths))
	{
		const size_t readi = pospair.first;
		const size_t i = pospair.second;
		if (i < readPaths[readi].path.size()-1 && readPaths[readi].path[i] == from && readPaths[readi].path[i+1] == to)
		{
			if (i < readPaths[readi].path.size()-2)
			{
				return resolvableGraph.unitigs[to.first].size();
			}
			for (const auto& path : readPaths[readi].reads)
			{
				assert(path.rightClip < resolvableGraph.unitigs[to.first].size());
				result = std::max(result, resolvableGraph.unitigs[to.first].size() - path.rightClip);
			}
		}
		if (i > 0 && readPaths[readi].path[i-1] == reverse(to) && readPaths[readi].path[i] == reverse(from))
		{
			if (i > 1)
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
	return result;
}

std::vector<ResolveTriplet> getGuessworkTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional)
{
	size_t fwDegree = resolvableGraph.edges[std::make_pair(node, true)].size();
	size_t bwDegree = resolvableGraph.edges[std::make_pair(node, false)].size();
	std::vector<ResolveTriplet> empty;
	if (fwDegree == 1 || bwDegree == 1) return empty;
	double nodeCoverage = resolvableGraph.getCoverage(readPaths, node);
	size_t nodeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
	if (nodeCopyCount < 2) return empty;
	size_t outneighborCopyCountSum = 0;
	size_t inneighborCopyCountSum = 0;
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> outneighborCopyCounts;
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> inneighborCopyCounts;
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		double nodeCoverage = resolvableGraph.getCoverage(readPaths, edge.first);
		size_t edgeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
		outneighborCopyCountSum += edgeCopyCount;
		outneighborCopyCounts.emplace_back(edge, edgeCopyCount);
		if (resolvables.count(edge.first) == 1) return empty;
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		double nodeCoverage = resolvableGraph.getCoverage(readPaths, edge.first);
		size_t edgeCopyCount = (size_t)((nodeCoverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage);
		inneighborCopyCountSum += edgeCopyCount;
		inneighborCopyCounts.emplace_back(edge, edgeCopyCount);
		if (resolvables.count(edge.first) == 1) return empty;
	}
	if (outneighborCopyCountSum != nodeCopyCount) return empty;
	if (inneighborCopyCountSum != nodeCopyCount) return empty;
	auto coveredTriplets = getRawTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, true);
	if (coveredTriplets.size() == 0) return empty;
	if (coveredTriplets.size() == 1 && minCoverage == 1)
	{
		auto test = getRawTriplets(resolvableGraph, resolvables, readPaths, node, 2, true);
		if (test.size() == 0) return empty;
	}
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredInNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		if (pair.left.first != std::numeric_limits<size_t>::max()) coveredInNeighbors.emplace(reverse(pair.left));
		if (pair.right.first != std::numeric_limits<size_t>::max()) coveredOutNeighbors.emplace(pair.right);
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
		if (coveredInNeighbors.count(pair.first) == 1) continue;
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
			if (inpair.first.first == outpair.first.first) return empty;
			coveredTriplets.emplace_back(reverse(inpair.first), outpair.first, 0);
			addedGuesses += 1;
		}
	}
	assert(addedGuesses > 0);
	return coveredTriplets;
}

void filterCopyCountTriplets(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const size_t node, std::vector<ResolveTriplet>& originals)
{
	double coverage = resolvableGraph.getCoverage(readPaths, node);
	double normalizedCoverage = coverage / resolvableGraph.averageCoverage;
	size_t copyCount = normalizedCoverage + 0.5;
	if (originals.size() <= copyCount) return;
	std::sort(originals.begin(), originals.end(), [](ResolveTriplet left, ResolveTriplet right) {
		return left.coverage > right.coverage;
	});
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredInNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> coveredOutNeighbors;
	phmap::flat_hash_set<std::pair<size_t, bool>> inneighborNeedsCovering;
	phmap::flat_hash_set<std::pair<size_t, bool>> outneighborNeedsCovering;
	phmap::flat_hash_set<std::pair<size_t, bool>> inneighborDoesntNeedCovering;
	phmap::flat_hash_set<std::pair<size_t, bool>> outneighborDoesntNeedCovering;
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		if (resolvableGraph.getCoverage(readPaths, edge.first) > resolvableGraph.averageCoverage * 0.25)
		{
			outneighborNeedsCovering.insert(edge);
		}
		else
		{
			outneighborDoesntNeedCovering.insert(edge);
		}
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		if (resolvableGraph.getCoverage(readPaths, edge.first) > resolvableGraph.averageCoverage * 0.25)
		{
			inneighborNeedsCovering.insert(reverse(edge));
		}
		else
		{
			inneighborDoesntNeedCovering.insert(reverse(edge));
		}
	}
	bool maybeDontRemoveLast = (inneighborDoesntNeedCovering.size() == 1) && (outneighborDoesntNeedCovering.size() == 1);
	size_t threshold = 0;
	for (size_t i = 0; i < copyCount; i++)
	{
		if (inneighborNeedsCovering.count(originals[i].left) == 1)
		{
			inneighborNeedsCovering.erase(originals[i].left);
		}
		if (outneighborNeedsCovering.count(originals[i].right) == 1)
		{
			outneighborNeedsCovering.erase(originals[i].right);
		}
		if (inneighborDoesntNeedCovering.count(originals[i].left) == 1)
		{
			inneighborDoesntNeedCovering.erase(originals[i].left);
		}
		if (outneighborDoesntNeedCovering.count(originals[i].right) == 1)
		{
			outneighborDoesntNeedCovering.erase(originals[i].right);
		}
		if (inneighborNeedsCovering.size() > 0) continue;
		if (outneighborNeedsCovering.size() > 0) continue;
		if (originals[i].coverage < originals[i+1].coverage * 4) continue;
		if (originals[i+1].coverage >= resolvableGraph.averageCoverage * 0.25) continue;
		if (originals[i].coverage < resolvableGraph.averageCoverage * 0.5) break;
		threshold = i+1;
		break;
	}
	if (threshold == 0) return;
	if (threshold == originals.size()) return;
	if (maybeDontRemoveLast && threshold == originals.size()-1)
	{
		if (inneighborDoesntNeedCovering.size() == 1 && outneighborDoesntNeedCovering.size() == 1)
		{
			if (originals[threshold].left == *inneighborDoesntNeedCovering.begin() && originals[threshold].right == *outneighborDoesntNeedCovering.begin())
			{
				// special case, the one extra triplet connects the one missing in- and outneighbors so it might be valid, keep it
				return;
			}
		}
	}
	originals.erase(originals.begin() + threshold, originals.end());
}

std::vector<ResolveTriplet> getValidTriplets(const ResolvableUnitigGraph& resolvableGraph, const phmap::flat_hash_set<size_t>& resolvables, const std::vector<PathGroup>& readPaths, size_t node, size_t minCoverage, bool unconditional, bool guesswork, const bool copycountFilterHeuristic)
{
	auto triplets = getReadSupportedTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork);
	if (triplets.size() > 0 && guesswork && unconditional && copycountFilterHeuristic)
	{
		filterCopyCountTriplets(resolvableGraph, readPaths, node, triplets);
	}
	else if (triplets.size() == 0 && guesswork)
	{
		triplets = getGuessworkTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional);
	}
	std::sort(triplets.begin(), triplets.end(), [](ResolveTriplet left, ResolveTriplet right) {
		assert(left.left != right.left || left.right != right.right);
		if (left.left.first < right.left.first) return true;
		if (left.left.first > right.left.first) return false;
		if (!left.left.second && right.left.second) return true;
		if (left.left.second && !right.left.second) return false;
		if (left.right.first < right.right.first) return true;
		if (left.right.first > right.right.first) return false;
		if (!left.right.second && right.right.second) return true;
		if (left.right.second && !right.right.second) return false;
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
		for (const auto pospair : resolvableGraph.iterateCrossingReads(node, readPaths)) relevantReads.insert(pospair.first);
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
					read.readPosStartIndex = read.readPosStartIndex + nodePosStarts[0] - read.leftClip;
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
					assert(extraClip < read.readPosEndIndex - read.readPosStartIndex);
					read.readPosEndIndex -= extraClip;
					read.rightClip = 0;
				}
				else
				{
					assert(read.rightClip >= (kmerPathLength - nodePosEnds.back()));
					read.rightClip -= (kmerPathLength - nodePosEnds.back());
				}
				assert(nodePosEnds.back() - nodePosStarts[0] == (read.readPosEndIndex - read.readPosStartIndex) + read.leftClip + read.rightClip);
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
					if (posesEnd < read.leftClip + (read.readPosEndIndex - read.readPosStartIndex))
					{
						rightClipRemove = read.rightClip;
						posesEnd -= read.leftClip;
					}
					else if (posesEnd < read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex))
					{
						rightClipRemove = (read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex)) - posesEnd;
						posesEnd = (read.readPosEndIndex - read.readPosStartIndex);
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
					if (posesEnd-posesStart > read.readPosEndIndex - read.readPosStartIndex)
					{
						posesEnd = read.readPosEndIndex - read.readPosStartIndex + posesStart;
					}
					assert(posesStart < posesEnd);
					assert(posesEnd <= (read.readPosEndIndex - read.readPosStartIndex));
					path.reads.emplace_back();
					path.reads.back().readInfoIndex = read.readInfoIndex;
					path.reads.back().readPosZeroOffset = read.readPosZeroOffset;
					path.reads.back().readPosStartIndex = read.readPosStartIndex + posesStart;
					path.reads.back().readPosEndIndex = read.readPosStartIndex + posesEnd;
					path.reads.back().leftClip = read.leftClip;
					path.reads.back().rightClip = read.rightClip;
					assert(path.reads.back().leftClip >= leftClipRemove);
					path.reads.back().leftClip -= leftClipRemove;
					assert(path.reads.back().rightClip >= rightClipRemove);
					path.reads.back().rightClip -= rightClipRemove;
					path.reads.back().readNameIndex = read.readNameIndex;
				}
				addPathButFirstMaybeTrim(resolvableGraph, readPaths, std::move(path));
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
			assert(posesEnd <= read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex));
			if (posesEnd < read.leftClip + (read.readPosEndIndex - read.readPosStartIndex))
			{
				rightClipRemove = read.rightClip;
				posesEnd -= read.leftClip;
			}
			else
			{
				rightClipRemove = (read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex)) - posesEnd;
				posesEnd = (read.readPosEndIndex - read.readPosStartIndex);
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
			assert(posesEnd <= (read.readPosEndIndex - read.readPosStartIndex));
			path.reads.emplace_back();
			path.reads.back().readInfoIndex = read.readInfoIndex;
			path.reads.back().readPosZeroOffset = read.readPosZeroOffset;
			path.reads.back().readPosStartIndex = read.readPosStartIndex + posesStart;
			path.reads.back().readPosEndIndex = read.readPosStartIndex + posesEnd;
			path.reads.back().leftClip = read.leftClip;
			path.reads.back().rightClip = read.rightClip;
			assert(path.reads.back().leftClip >= leftClipRemove);
			path.reads.back().leftClip -= leftClipRemove;
			assert(path.reads.back().rightClip >= rightClipRemove);
			path.reads.back().rightClip -= rightClipRemove;
			path.reads.back().readNameIndex = read.readNameIndex;
		}
		addPathButFirstMaybeTrim(resolvableGraph, readPaths, std::move(path));
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

ResolutionResult resolve(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, std::vector<PathGroup>& readPaths, const phmap::flat_hash_set<size_t>& resolvables, const size_t minCoverage, const bool unconditional, const bool guesswork, const bool copycountFilterHeuristic)
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
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork, copycountFilterHeuristic);
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
			phmap::flat_hash_map<std::pair<size_t, bool>, size_t> fakeEdgeCount;
			auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork, copycountFilterHeuristic);
			bool unresolve = false;
			for (auto triplet : triplets)
			{
				bool bwFake = false;
				bool fwFake = false;
				if (triplet.left.first != std::numeric_limits<size_t>::max() && (resolvables.count(triplet.left.first) == 0 || unresolvables.count(triplet.left.first) == 1))
				{
					if (resolvableGraph.unitigLength(triplet.left.first) == resolvableGraph.getBpOverlap(std::make_pair(node, false), reverse(triplet.left))+1)
					{
						bwFake = true;
					}
				}
				if (triplet.right.first != std::numeric_limits<size_t>::max() && (resolvables.count(triplet.right.first) == 0 || unresolvables.count(triplet.right.first) == 1))
				{
					if (resolvableGraph.unitigLength(triplet.right.first) == resolvableGraph.getBpOverlap(std::make_pair(node, true), triplet.right)+1)
					{
						fwFake = true;
					}
				}
				if (bwFake)
				{
					if (resolvableGraph.edges[triplet.left].size() >= 2) unresolve = true;
				}
				if (fwFake)
				{
					if (resolvableGraph.edges[reverse(triplet.right)].size() >= 2) unresolve = true;
				}
				if (bwFake && fwFake) unresolve = true;
			}
			if (unresolve)
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
	std::unordered_map<size_t, std::vector<ResolveTriplet>> tripletsPerNode;
	for (auto node : actuallyResolvables)
	{
		tripletsPerNode[node] = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage, unconditional, guesswork, copycountFilterHeuristic);
	}
	assert(actuallyResolvables.activeSize() == resolvables.size() - unresolvables.size());
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newEdgeNodes;
	for (auto node : actuallyResolvables)
	{
		assert(tripletsPerNode[node].size() > 0);
		std::unordered_set<std::pair<size_t, bool>> fwCovered;
		std::unordered_set<std::pair<size_t, bool>> bwCovered;
		for (auto triplet : tripletsPerNode[node])
		{
			fwCovered.insert(triplet.right);
			bwCovered.insert(reverse(triplet.left));
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
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1)
			{
				createFakeEdgeNode(resolvableGraph, hashlist, newEdgeNodes, resolvables, unresolvables, pos, edge);
			}
			else
			{
				assert(actuallyResolvables.get(edge.first) || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
				createEdgeNode(resolvableGraph, hashlist, newEdgeNodes, resolvables, unresolvables, pos, edge);
			}
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
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1)
			{
				createFakeEdgeNode(resolvableGraph, hashlist, newEdgeNodes, resolvables, unresolvables, pos, edge);
			}
			else
			{
				assert(actuallyResolvables.get(edge.first) || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
				createEdgeNode(resolvableGraph, hashlist, newEdgeNodes, resolvables, unresolvables, pos, edge);
			}
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
			for (auto triplet : tripletsPerNode[fromnode.first])
			{
				if (fromnode.second && triplet.right == tonode)
				{
					fromHasTo = true;
				}
				if (!fromnode.second && reverse(triplet.left) == tonode)
				{
					fromHasTo = true;
				}
			}
			for (auto triplet : tripletsPerNode[tonode.first])
			{
				if (tonode.second && triplet.right == reverse(fromnode))
				{
					toHasFrom = true;
				}
				if (!tonode.second && reverse(triplet.left) == reverse(fromnode))
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
		auto triplets = tripletsPerNode[node];
		assert(triplets.size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto triplet : triplets)
		{
			const std::pair<size_t, bool> before = triplet.left;
			const std::pair<size_t, bool> after = triplet.right;
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
			assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 0 || newEdgeNodes.count(std::make_pair(before, pos)) == 0);
			assert(newEdgeNodes.count(std::make_pair(pos, after)) == 0 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(before, pos)) == 0);
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(pos), reverse(before))), false);
			}
			else
			{
				assert(newEdgeNodes.count(std::make_pair(before, pos)) == 1);
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(before, pos)), true);
			}
			if (newEdgeNodes.count(std::make_pair(pos, after)) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(pos, after)), true);
			}
			else
			{
				assert(newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1);
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(after), reverse(pos))), false);
			}
			result.maybeUnitigifiable.insert(leftNode.first);
			result.maybeUnitigifiable.insert(rightNode.first);
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
			assert(((resolvables.count(before.first) == 0 || unresolvables.count(before.first) == 1) && resolvableGraph.unitigLength(before.first) == resolvableGraph.getBpOverlap(reverse(pos), reverse(before)) + 1) || resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(leftNode.first));
			assert(((resolvables.count(after.first) == 0 || unresolvables.count(after.first) == 1) && resolvableGraph.unitigLength(after.first) == resolvableGraph.getBpOverlap(pos, after) + 1) || resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(rightNode.first));
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

size_t getEdgeCoverage(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	size_t result = 0;
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(from.first, readPaths))
	{
		const size_t i = pospair.first;
		const size_t j = pospair.second;
		if (j < readPaths[i].path.size()-1 && readPaths[i].path[j] == from && readPaths[i].path[j+1] == to)
		{
			result += readPaths[i].reads.size();
		}
		if (j > 0 && readPaths[i].path[j] == reverse(from) && readPaths[i].path[j-1] == reverse(to))
		{
			result += readPaths[i].reads.size();
		}
	}
	return result;
}

void checkValidity(const ResolvableUnitigGraph& graph, const std::vector<PathGroup>& readPaths)
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
			assert(graph.unitigLeftClipBp[i] < graph.kmerSize - graph.hashlist.getOverlap(graph.unitigs[i][0], graph.unitigs[i][1]));
			assert(graph.unitigRightClipBp[i] < graph.kmerSize - graph.hashlist.getOverlap(graph.unitigs[i][graph.unitigs[i].size()-2], graph.unitigs[i][graph.unitigs[i].size()-1]));
		}
		assert(graph.unitigLeftClipBp[i] < graph.kmerSize);
		assert(graph.unitigRightClipBp[i] < graph.kmerSize);
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges[fw])
		{
			assert(getEdgeCoverage(graph, readPaths, fw, edge) > 0);
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(fw)) == 1);
			assert(graph.edges[fw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2 || edge.first == i);
			assert(graph.getBpOverlap(fw, edge) < graph.unitigLength(i));
			assert(graph.getBpOverlap(fw, edge) < graph.unitigLength(edge.first));
		}
		for (auto edge : graph.edges[bw])
		{
			assert(getEdgeCoverage(graph, readPaths, bw, edge) > 0);
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
		// for (const auto& read : path.reads)
		// {
		// 	if (read.readPoses.size() > 0)
		// 	{
		// 		for (size_t i = 1; i < read.readPoses.size(); i++)
		// 		{
		// 			assert(read.readPoses[i-1] < read.readPoses[i]);
		// 		}
		// 		assert(read.readPoses.back() + graph.kmerSize <= read.readLengthHPC);
		// 	}
		// }
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
				assert(read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex) == pathHashCount);
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
			assert(read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex) == pathHashCount);
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
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		if (i == graph.precalcedUnitigCoverages.size()) break;
		if (graph.unitigRemoved[i]) continue;
		if (graph.precalcedUnitigCoverages[i] == 0) continue;
		double recalcedCoverage = graph.calculateCoverage(readPaths, i);
		if (!(graph.precalcedUnitigCoverages[i] > recalcedCoverage - 0.01 && graph.precalcedUnitigCoverages[i] < recalcedCoverage + 0.01))
		{
			std::cerr << graph.precalcedUnitigCoverages[i] << " " << recalcedCoverage << std::endl;
		}
		assert(graph.precalcedUnitigCoverages[i] > recalcedCoverage - 0.01 && graph.precalcedUnitigCoverages[i] < recalcedCoverage + 0.01);
	}
}

void removeEdgesAndNodes(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::unordered_set<size_t>& removeNodes, const std::unordered_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>& removeEdges)
{
	std::unordered_set<size_t> relevantReads;
	for (auto edge : removeEdges)
	{
		assert(edge == canon(edge.first, edge.second));
		assert(resolvableGraph.edges[edge.first].count(edge.second) == 1);
		assert(resolvableGraph.edges[reverse(edge.second)].count(reverse(edge.first)) == 1);
		resolvableGraph.edges[edge.first].erase(edge.second);
		if (edge.first != reverse(edge.second))
		{
			assert(resolvableGraph.edges[reverse(edge.second)].count(reverse(edge.first)) == 1);
			resolvableGraph.edges[reverse(edge.second)].erase(reverse(edge.first));
		}
		else
		{
			assert(resolvableGraph.edges[reverse(edge.second)].count(reverse(edge.first)) == 0);
		}
		for (const auto pospair : resolvableGraph.iterateCrossingReads(edge.first.first, readPaths)) relevantReads.insert(pospair.first);
		for (const auto pospair : resolvableGraph.iterateCrossingReads(edge.second.first, readPaths)) relevantReads.insert(pospair.first);
	}
	for (auto node : removeNodes)
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
		for (const auto pospair : resolvableGraph.iterateCrossingReads(node, readPaths)) relevantReads.insert(pospair.first);
	}
	std::vector<size_t> relevantReadsDeterministic { relevantReads.begin(), relevantReads.end() };
	std::sort(relevantReadsDeterministic.begin(), relevantReadsDeterministic.end());
	for (const size_t i : relevantReadsDeterministic)
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
			assert(pathKmerLength == (read.readPosEndIndex - read.readPosStartIndex) + read.leftClip + read.rightClip);
		}
		size_t lastStart = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (removeNodes.count(readPaths[i].path[j].first) == 1 || (j > 0 && removeEdges.count(canon(readPaths[i].path[j-1], readPaths[i].path[j])) == 1))
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
					if (posesEnd < read.leftClip + (read.readPosEndIndex - read.readPosStartIndex))
					{
						rightClipRemove = read.rightClip;
						posesEnd -= read.leftClip;
					}
					else
					{
						assert(posesEnd <= read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex));
						rightClipRemove = (read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex)) - posesEnd;
						posesEnd = (read.readPosEndIndex - read.readPosStartIndex);
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
					assert(posesEnd <= (read.readPosEndIndex - read.readPosStartIndex));
					path.reads.emplace_back();
					path.reads.back().readInfoIndex = read.readInfoIndex;
					path.reads.back().readPosZeroOffset = read.readPosZeroOffset;
					path.reads.back().readPosStartIndex = read.readPosStartIndex + posesStart;
					path.reads.back().readPosEndIndex = read.readPosStartIndex + posesEnd;
					path.reads.back().leftClip = read.leftClip;
					path.reads.back().rightClip = read.rightClip;
					assert(path.reads.back().leftClip >= leftClipRemove);
					path.reads.back().leftClip -= leftClipRemove;
					assert(path.reads.back().rightClip >= rightClipRemove);
					path.reads.back().rightClip -= rightClipRemove;
					path.reads.back().readNameIndex = read.readNameIndex;
				}
				addPathButFirstMaybeTrim(resolvableGraph, readPaths, std::move(path));
				lastStart = j;
				if (removeNodes.count(readPaths[i].path[j].first) == 1) lastStart = j+1;
			}
		}
		if (lastStart > 0 && lastStart < readPaths[i].path.size())
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
				assert(posesEnd <= read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex));
				if (posesEnd < read.leftClip + (read.readPosEndIndex - read.readPosStartIndex))
				{
					rightClipRemove = read.rightClip;
					posesEnd -= read.leftClip;
				}
				else
				{
					assert(posesEnd <= read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex));
					rightClipRemove = (read.leftClip + read.rightClip + (read.readPosEndIndex - read.readPosStartIndex)) - posesEnd;
					posesEnd = (read.readPosEndIndex - read.readPosStartIndex);
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
				assert(posesEnd <= (read.readPosEndIndex - read.readPosStartIndex));
				path.reads.emplace_back();
				path.reads.back().readInfoIndex = read.readInfoIndex;
				path.reads.back().readPosZeroOffset = read.readPosZeroOffset;
				path.reads.back().readPosStartIndex = read.readPosStartIndex + posesStart;
				path.reads.back().readPosEndIndex = read.readPosStartIndex + posesEnd;
				path.reads.back().leftClip = read.leftClip;
				path.reads.back().rightClip = read.rightClip;
				assert(path.reads.back().leftClip >= leftClipRemove);
				path.reads.back().leftClip -= leftClipRemove;
				assert(path.reads.back().rightClip >= rightClipRemove);
				path.reads.back().rightClip -= rightClipRemove;
				path.reads.back().readNameIndex = read.readNameIndex;
			}
			addPathButFirstMaybeTrim(resolvableGraph, readPaths, std::move(path));
		}
		if (lastStart != 0) erasePath(resolvableGraph, readPaths, i);
	}
	for (auto node : removeNodes)
	{
		resolvableGraph.readsCrossingNode[node].clear();
		resolvableGraph.unitigRemoved[node] = true;
	}
}

void removeNode(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, size_t node)
{
	std::unordered_set<size_t> oneNode;
	oneNode.insert(node);
	removeEdgesAndNodes(resolvableGraph, readPaths, oneNode, std::unordered_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>{});
}

struct UntippingResult
{
	UntippingResult() :
		nodesRemoved(0),
		edgesRemoved(0),
		maybeUnitigifiable()
	{
	}
	size_t nodesRemoved;
	size_t edgesRemoved;
	phmap::flat_hash_set<size_t> maybeUnitigifiable;
};

void tryRemoveTip(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const size_t i, UntippingResult& result)
{
	assert(!resolvableGraph.unitigRemoved[i]);
	bool fwHasSafeEdge = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)])
	{
		if (resolvableGraph.getCoverage(readPaths, edge.first) < minSafeCoverage) return;
		if (getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, true), edge) > maxRemovableCoverage) return;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) fwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, true)].size() > 0 && !fwHasSafeEdge) return;
	bool bwHasSafeEdge = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)])
	{
		if (resolvableGraph.getCoverage(readPaths, edge.first) < minSafeCoverage) return;
		if (getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, false), edge) > maxRemovableCoverage) return;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) bwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, false)].size() > 0 && !bwHasSafeEdge) return;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)]) result.maybeUnitigifiable.insert(edge.first);
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)]) result.maybeUnitigifiable.insert(edge.first);
	removeNode(resolvableGraph, readPaths, i);
	result.nodesRemoved += 1;
}

UntippingResult tryRemoveCrosslinks(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const double maxRemovableCoverage, const double minSafeCoverage, const std::pair<size_t, bool> start)
{
	UntippingResult result;
	assert(resolvableGraph.edges[start].size() >= 2);
	bool possibleToRemove = false;
	for (auto edge : resolvableGraph.edges[start])
	{
		if (resolvableGraph.edges[reverse(edge)].size() >= 2)
		{
			possibleToRemove = true;
		}
	}
	if (!possibleToRemove) return result;
	bool hasSafe = false;
	std::vector<std::pair<size_t, bool>> checkThese;
	for (auto edge : resolvableGraph.edges[start])
	{
		size_t coverage = getEdgeCoverage(resolvableGraph, readPaths, start, edge);
		if (coverage >= minSafeCoverage)
		{
			hasSafe = true;
		}
		if (coverage <= maxRemovableCoverage && resolvableGraph.edges[reverse(edge)].size() >= 2)
		{
			checkThese.push_back(edge);
		}
	}
	if (!hasSafe) return result;
	std::vector<std::pair<size_t, bool>> removeThese;
	for (auto edge : checkThese)
	{
		bool otherHasSafe = false;
		assert(resolvableGraph.edges[reverse(edge)].size() >= 2);
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (edge2 == reverse(start)) continue;
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage)
			{
				otherHasSafe = true;
			}
		}
		if (otherHasSafe) removeThese.push_back(edge);
	}
	if (removeThese.size() == 0) return result;
	std::unordered_set<size_t> removedNodes;
	std::unordered_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> removedEdges;
	for (auto edge : removeThese)
	{
		removedEdges.emplace(start, edge);
		result.maybeUnitigifiable.insert(edge.first);
	}
	result.maybeUnitigifiable.insert(start.first);
	result.edgesRemoved = removeThese.size();
	removeEdgesAndNodes(resolvableGraph, readPaths, removedNodes, removedEdges);
	return result;
}

UntippingResult removeLowCoverageCrosslinks(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const double maxRemovableCoverage, const double minSafeCoverage)
{
	UntippingResult result;
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2)
		{
			auto partResult = tryRemoveCrosslinks(resolvableGraph, readPaths, maxRemovableCoverage, minSafeCoverage, std::make_pair(i, true));
			result.nodesRemoved += partResult.nodesRemoved;
			result.edgesRemoved += partResult.edgesRemoved;
			result.maybeUnitigifiable.insert(partResult.maybeUnitigifiable.begin(), partResult.maybeUnitigifiable.end());
		}
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2)
		{
			auto partResult = tryRemoveCrosslinks(resolvableGraph, readPaths, maxRemovableCoverage, minSafeCoverage, std::make_pair(i, false));
			result.nodesRemoved += partResult.nodesRemoved;
			result.edgesRemoved += partResult.edgesRemoved;
			result.maybeUnitigifiable.insert(partResult.maybeUnitigifiable.begin(), partResult.maybeUnitigifiable.end());
		}
	}
	return result;
}

UntippingResult removeLowCoverageTips(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const phmap::flat_hash_set<size_t>& maybeUntippable)
{
	for (size_t i = resolvableGraph.lastTippableChecked; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = resolvableGraph.getCoverage(readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	std::vector<size_t> deterministicMaybeUntippable { maybeUntippable.begin(), maybeUntippable.end() };
	std::sort(deterministicMaybeUntippable.begin(), deterministicMaybeUntippable.end());
	for (auto i : deterministicMaybeUntippable)
	{
		if (i >= resolvableGraph.lastTippableChecked) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = resolvableGraph.getCoverage(readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	resolvableGraph.lastTippableChecked = resolvableGraph.unitigs.size();
	UntippingResult result;
	for (size_t index = resolvableGraph.everTippable.size()-1; index < resolvableGraph.everTippable.size(); index--)
	{
		size_t i = resolvableGraph.everTippable[index];
		if (resolvableGraph.unitigRemoved[i])
		{
			std::swap(resolvableGraph.everTippable[index], resolvableGraph.everTippable.back());
			resolvableGraph.everTippable.pop_back();
			continue;
		}
		tryRemoveTip(resolvableGraph, readPaths, hashlist, maxRemovableCoverage, minSafeCoverage, maxRemovableLength, i, result);
	}
	return result;
}

struct CleaningResult
{
	CleaningResult() : 
		nodesRemoved(0),
		edgesRemoved(0),
		maybeUnitigifiable(),
		checked()
	{}
	size_t nodesRemoved;
	size_t edgesRemoved;
	phmap::flat_hash_set<size_t> maybeUnitigifiable;
	phmap::flat_hash_set<std::pair<size_t, bool>> checked;
};

CleaningResult cleanComponent(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t minLongLength, const size_t minUnresolvableLength, const size_t maxUnresolvableLength, std::pair<size_t, bool> start)
{
	CleaningResult result;
	std::vector<std::pair<size_t, bool>> stack;
	stack.push_back(start);
	std::vector<std::pair<size_t, bool>> componentNodeSides;
	bool allValid = true;
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (result.checked.count(top) == 1) continue;
		result.checked.insert(top);
		size_t length = resolvableGraph.unitigLength(top.first);
		if (length > minUnresolvableLength && length < maxUnresolvableLength) allValid = false;
		componentNodeSides.push_back(top);
		double coverage = resolvableGraph.getCoverage(readPaths, top.first);
		size_t estimatedCopyCount = (coverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage;
		if (resolvableGraph.unitigLength(top.first) < minLongLength || estimatedCopyCount == 0) stack.push_back(reverse(top));
		for (auto edge : resolvableGraph.edges[top]) stack.push_back(reverse(edge));
	}
	if (!allValid) return result;
	for (auto pair : componentNodeSides)
	{
		double coverage = resolvableGraph.getCoverage(readPaths, pair.first);
		size_t estimatedCopyCount = (coverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage;
		if (estimatedCopyCount <= 1 && coverage < resolvableGraph.averageCoverage * ((double)estimatedCopyCount - 0.4))
		{
			allValid = false;
			break;
		}
		if (estimatedCopyCount <= 1 && coverage > resolvableGraph.averageCoverage * ((double)estimatedCopyCount + 0.4))
		{
			allValid = false;
			break;
		}
		size_t edgeCopyCountSum = 0;
		for (auto edge : resolvableGraph.edges[pair])
		{
			int coverage = getEdgeCoverage(resolvableGraph, readPaths, pair, edge);
			size_t edgeCopyCount = (coverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage;
			if (estimatedCopyCount <= 1 && coverage < resolvableGraph.averageCoverage * ((double)edgeCopyCount - 0.4))
			{
				allValid = false;
				break;
			}
			if (estimatedCopyCount <= 1 && coverage > resolvableGraph.averageCoverage * ((double)edgeCopyCount + 0.4))
			{
				allValid = false;
				break;
			}
			edgeCopyCountSum += edgeCopyCount;
		}
		if (!allValid) break;
		if (edgeCopyCountSum != estimatedCopyCount)
		{
			allValid = false;
			break;
		}
	}
	if (!allValid)
	{
		return result;
	}
	std::unordered_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> removedEdges;
	std::unordered_set<size_t> removedNodes;
	for (auto pair : componentNodeSides)
	{
		for (auto edge : resolvableGraph.edges[pair])
		{
			int coverage = getEdgeCoverage(resolvableGraph, readPaths, pair, edge);
			size_t edgeCopyCount = (coverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage;
			if (edgeCopyCount == 0) removedEdges.insert(canon(pair, edge));
		}
		double coverage = resolvableGraph.getCoverage(readPaths, pair.first);
		size_t estimatedCopyCount = (coverage + resolvableGraph.averageCoverage / 2) / resolvableGraph.averageCoverage;
		if (estimatedCopyCount == 0) removedNodes.insert(pair.first);
	}
	for (auto edge : removedEdges)
	{
		if (removedNodes.count(edge.first.first) == 0) result.maybeUnitigifiable.insert(edge.first.first);
		if (removedNodes.count(edge.second.first) == 0) result.maybeUnitigifiable.insert(edge.second.first);
	}
	result.nodesRemoved = removedNodes.size();
	result.edgesRemoved = removedEdges.size();
	removeEdgesAndNodes(resolvableGraph, readPaths, removedNodes, removedEdges);
	return result;
}

CleaningResult cleanComponentsByCopynumber(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t minLongLength, const size_t minUnresolvableLength, const size_t maxUnresolvableLength, const phmap::flat_hash_set<size_t>& checkThese, const size_t minNew)
{
	std::unordered_set<std::pair<size_t, bool>> checked;
	CleaningResult result;
	std::vector<size_t> checkTheseDeterministic { checkThese.begin(), checkThese.end() };
	std::sort(checkTheseDeterministic.begin(), checkTheseDeterministic.end());
	for (auto i : checkTheseDeterministic)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		if (checked.count(std::make_pair(i, true)) == 0)
		{
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, true));
			for (auto node : part.checked)
			{
				if (checkThese.count(node.first) == 1 || node.first >= minNew) checked.insert(node);
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
		if (checked.count(std::make_pair(i, false)) == 0)
		{
			assert(!resolvableGraph.unitigRemoved[i]);
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, false));
			for (auto node : part.checked)
			{
				if (checkThese.count(node.first) == 1 || node.first >= minNew) checked.insert(node);
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
	}
	for (size_t i = minNew; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		if (checked.count(std::make_pair(i, true)) == 0)
		{
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, true));
			for (auto node : part.checked)
			{
				if (checkThese.count(node.first) == 1 || node.first >= minNew) checked.insert(node);
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
		if (checked.count(std::make_pair(i, false)) == 0)
		{
			assert(!resolvableGraph.unitigRemoved[i]);
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, false));
			for (auto node : part.checked)
			{
				if (checkThese.count(node.first) == 1 || node.first >= minNew) checked.insert(node);
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
	}
	return result;
}

CleaningResult cleanComponentsByCopynumber(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const size_t minLongLength, const size_t minUnresolvableLength, const size_t maxUnresolvableLength)
{
	CleaningResult result;
	std::vector<bool> fwChecked;
	std::vector<bool> bwChecked;
	fwChecked.resize(resolvableGraph.unitigs.size(), false);
	bwChecked.resize(resolvableGraph.unitigs.size(), false);
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		if (!fwChecked[i])
		{
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, true));
			for (auto node : part.checked)
			{
				if (node.second) fwChecked[node.first] = true;
				if (!node.second) bwChecked[node.first] = true;
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
		if (!bwChecked[i])
		{
			assert(!resolvableGraph.unitigRemoved[i]);
			auto part = cleanComponent(resolvableGraph, readPaths, minLongLength, minUnresolvableLength, maxUnresolvableLength, std::make_pair(i, false));
			for (auto node : part.checked)
			{
				if (node.second) fwChecked[node.first] = true;
				if (!node.second) bwChecked[node.first] = true;
			}
			for (auto node : part.maybeUnitigifiable)
			{
				result.maybeUnitigifiable.insert(node);
			}
			result.nodesRemoved += part.nodesRemoved;
			result.edgesRemoved += part.edgesRemoved;
		}
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
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(pos.first, readPaths))
	{
		const size_t i = pospair.first;
		const size_t j = pospair.second;
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
		if (maxReadTrim == 0) break;
	}
	if (maxReadTrim < trimAmount) return false;
	return true;
}

void trimEndRecursive(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const std::pair<size_t, bool> pos, const size_t trimAmount)
{
	if (resolvableGraph.precalcedUnitigLengths.size() > pos.first) resolvableGraph.precalcedUnitigLengths[pos.first] = 0;
	assert(trimAmount > 0);
	assert(resolvableGraph.unitigs[pos.first].size() > trimAmount);
	if (resolvableGraph.precalcedUnitigCoverages.size() > pos.first) resolvableGraph.precalcedUnitigCoverages[pos.first] = 0;
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
			trimEndRecursive(resolvableGraph, readPaths, reverse(edge), trimThere);
			assert(resolvableGraph.overlaps.at(canon(reverse(pos), edge)) == resolvableGraph.unitigs[pos.first].size() - trimAmount - 1);
		}
		assert(resolvableGraph.overlaps.at(canon(reverse(pos), edge)) < resolvableGraph.unitigs[pos.first].size() - trimAmount);
	}
	if (pos.second)
	{
		assert(resolvableGraph.unitigRightClipBp[pos.first] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[pos.first][resolvableGraph.unitigs[pos.first].size()-2], resolvableGraph.unitigs[pos.first].back()));
		resolvableGraph.unitigs[pos.first].erase(resolvableGraph.unitigs[pos.first].end() - trimAmount, resolvableGraph.unitigs[pos.first].end());
		resolvableGraph.unitigRightClipBp[pos.first] = 0;
	}
	else
	{
		assert(resolvableGraph.unitigLeftClipBp[pos.first] < resolvableGraph.kmerSize - resolvableGraph.hashlist.getOverlap(resolvableGraph.unitigs[pos.first][0], resolvableGraph.unitigs[pos.first][1]));
		resolvableGraph.unitigs[pos.first].erase(resolvableGraph.unitigs[pos.first].begin(), resolvableGraph.unitigs[pos.first].begin() + trimAmount);
		resolvableGraph.unitigLeftClipBp[pos.first] = 0;
	}
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(pos.first, readPaths))
	{
		const size_t i = pospair.first;
		const size_t j = pospair.second;
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

void maybeTrim(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, std::pair<size_t, bool> pos, size_t maxTrim)
{
	if (maxTrim == 0) return;
	if (resolvableGraph.edges[pos].size() > 0) return;
	assertPrintReads(!resolvableGraph.unitigRemoved[pos.first], resolvableGraph, readPaths, pos.first);
	size_t maxReadTrim = resolvableGraph.unitigs[pos.first].size();
	assertPrintReads(resolvableGraph.readsCrossingNode[pos.first].size() > 0, resolvableGraph, readPaths, pos.first);
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(pos.first, readPaths))
	{
		const size_t i = pospair.first;
		const size_t j = pospair.second;
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
			maybeTrim(resolvableGraph, readPaths, pair.first, pair.second);
		}
		return;
	}
	trimEndRecursive(resolvableGraph, readPaths, pos, maxReadTrim);
}

void trimNodes(std::vector<std::pair<std::pair<size_t, bool>, size_t>> maybeTrimmable, ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths)
{
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
		maybeTrim(resolvableGraph, readPaths, pair.first, pair.second);
	}
}

void addPlusOneComponent(const ResolvableUnitigGraph& resolvableGraph, phmap::flat_hash_set<size_t>& resolvables, size_t startNode, size_t maxLength)
{
	phmap::flat_hash_set<size_t> checked;
	std::vector<size_t> checkQueue;
	checkQueue.push_back(startNode);
	while (checkQueue.size() > 0)
	{
		size_t node = checkQueue.back();
		checkQueue.pop_back();
		if (checked.count(node) == 1) continue;
		checked.emplace(node);
		for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
		{
			if (resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(std::make_pair(node, true), edge)+1)
			{
				if (resolvableGraph.unitigLength(edge.first) > maxLength)
				{
					return;
				}
				checkQueue.push_back(edge.first);
			}
		}
		for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
		{
			if (resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(std::make_pair(node, false), edge)+1)
			{
				if (resolvableGraph.unitigLength(edge.first) > maxLength)
				{
					return;
				}
				checkQueue.push_back(edge.first);
			}
		}
	}
	if (checked.size() == 1) return;
	resolvables.insert(checked.begin(), checked.end());
}

bool isLocallyRepetitive(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const size_t node)
{
	for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(node, readPaths))
	{
		const size_t pathi = pospair.first;
		if (readPaths[pathi].path.size() < 2) continue;
		size_t nodeCount = 0;
		for (std::pair<size_t, bool> pathnode : readPaths[pathi].path)
		{
			if (pathnode.first == node) nodeCount += 1;
		}
		if (nodeCount >= 2)
		{
			return true;
		}
	}
	return false;
}

phmap::flat_hash_set<size_t> filterToOnlyLocallyRepetitives(const ResolvableUnitigGraph& resolvableGraph, const std::vector<PathGroup>& readPaths, const phmap::flat_hash_set<size_t>& unfilteredResolvables)
{
	phmap::flat_hash_set<size_t> result;
	for (size_t node : unfilteredResolvables)
	{
		for (const std::pair<uint32_t, uint32_t> pospair : resolvableGraph.iterateCrossingReads(node, readPaths))
		{
			const size_t pathi = pospair.first;
			if (readPaths[pathi].path.size() < 2) continue;
			size_t nodeCount = 0;
			for (std::pair<size_t, bool> pathnode : readPaths[pathi].path)
			{
				if (pathnode.first == node) nodeCount += 1;
			}
			if (nodeCount >= 2)
			{
				result.insert(node);
				break;
			}
		}
	}
	return result;
}

void resolveRound(ResolvableUnitigGraph& resolvableGraph, std::vector<PathGroup>& readPaths, const HashList& hashlist, const size_t minCoverage, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool guesswork, const bool copycountFilterHeuristic, const bool onlyLocalResolve)
{
	checkValidity(resolvableGraph, readPaths);
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
		phmap::flat_hash_set<size_t> thisLengthNodes;
		while (queue.size() > 0 && resolvableGraph.unitigLength(queue.top()) == topSize)
		{
			if (!resolvableGraph.unitigRemoved[queue.top()])
			{
				addPlusOneComponent(resolvableGraph, resolvables, queue.top(), topSize);
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
			thisLengthNodes.insert(queue.top());
			queue.pop();
		}
		if (onlyLocalResolve)
		{
			resolvables = filterToOnlyLocallyRepetitives(resolvableGraph, readPaths, resolvables);
			auto oldResolvables = resolvables;
			for (auto node : oldResolvables)
			{
				addPlusOneComponent(resolvableGraph, resolvables, node, topSize);
			}
		}
		if (resolvables.size() == 0) continue;
		assert(resolvables.size() > 0);
		size_t oldSize = resolvableGraph.unitigs.size();
		checkValidity(resolvableGraph, readPaths);
		std::cerr << "try resolve k=" << topSize;
		auto resolutionResult = resolve(resolvableGraph, hashlist, readPaths, resolvables, minCoverage, topSize < maxUnconditionalResolveLength, guesswork, copycountFilterHeuristic);
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
			trimNodes(maybeTrimmable, resolvableGraph, readPaths);
		}
		if (guesswork)
		{
			phmap::flat_hash_set<size_t> cleanables = thisLengthNodes;
			cleanables.insert(resolutionResult.maybeUnitigifiable.begin(), resolutionResult.maybeUnitigifiable.end());
			auto removed = cleanComponentsByCopynumber(resolvableGraph, readPaths, 50000, topSize, 0, cleanables, oldSize);
			if (removed.nodesRemoved > 0 || removed.edgesRemoved > 0)
			{
				std::cerr << "removed " << removed.nodesRemoved << " nodes and " << removed.edgesRemoved << " edges" << std::endl;
				std::vector<size_t> maybeUnitigifiable;
				maybeUnitigifiable.insert(maybeUnitigifiable.end(), removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end());
				std::sort(maybeUnitigifiable.begin(), maybeUnitigifiable.end());
				auto unitigifiedHere = unitigifySet(resolvableGraph, readPaths, maybeUnitigifiable);
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
		}
		auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, resolutionResult.maybeUnitigifiable);
		resolutionResult.maybeUnitigifiable.insert(removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end());
		auto removed2 = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 2, 5, 10000, resolutionResult.maybeUnitigifiable);
		auto removed3 = removeLowCoverageCrosslinks(resolvableGraph, readPaths, 1, 5);
		auto removed4 = removeLowCoverageCrosslinks(resolvableGraph, readPaths, 2, 10);
		nodesRemoved += removed.nodesRemoved + removed2.nodesRemoved;
		if (removed.nodesRemoved + removed2.nodesRemoved > 0 || removed.edgesRemoved + removed2.edgesRemoved + removed3.edgesRemoved + removed4.edgesRemoved > 0)
		{
			std::cerr << "removed " << removed.nodesRemoved + removed2.nodesRemoved << " tips and " << removed.edgesRemoved + removed2.edgesRemoved + removed3.edgesRemoved + removed4.edgesRemoved << " edges" << std::endl;
			removed.maybeUnitigifiable.insert(removed2.maybeUnitigifiable.begin(), removed2.maybeUnitigifiable.end());
			std::vector<size_t> maybeUnitigifiable2 { removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end() };
			maybeUnitigifiable2.insert(maybeUnitigifiable2.end(), removed3.maybeUnitigifiable.begin(), removed3.maybeUnitigifiable.end());
			maybeUnitigifiable2.insert(maybeUnitigifiable2.end(), removed4.maybeUnitigifiable.begin(), removed4.maybeUnitigifiable.end());
			std::sort(maybeUnitigifiable2.begin(), maybeUnitigifiable2.end());
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
					if (resolvableGraph.precalcedUnitigCoverages.size() > unitigifiedHere[i].second) resolvableGraph.precalcedUnitigCoverages[unitigifiedHere[i].second] = 0;
				}
			}
		}
		checkValidity(resolvableGraph, readPaths);
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
			checkValidity(resolvableGraph, readPaths);
			for (auto node : queueNodes)
			{
				assert(node < resolvableGraph.unitigs.size());
				assert(!resolvableGraph.unitigRemoved[node]);
				queue.emplace(node);
			}
		}
	}
}

bool operator==(const std::vector<std::pair<size_t, bool>>& left, const std::vector<Node>& right)
{
	if (left.size() != right.size()) return false;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != (std::pair<size_t, bool>)right[i]) return false;
	}
	return true;
}

bool operator!=(const std::vector<std::pair<size_t, bool>>& left, const std::vector<Node>& right)
{
	return !(left == right);
}

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, std::vector<ReadPath>& rawReadPaths, const ReadpartIterator& partIterator, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool keepGaps, const bool guesswork, const bool copycountFilterHeuristic, const bool onlyLocalResolve)
{
	auto resolvableGraph = getUnitigs(initial, minCoverage, hashlist, kmerSize, keepGaps);
	std::cerr << rawReadPaths.size() << " raw read paths" << std::endl;
	// todo maybe fix? or does it matter?
	cutRemovedEdgesFromPaths(resolvableGraph, rawReadPaths);
	std::sort(rawReadPaths.begin(), rawReadPaths.end(), [](const ReadPath& left, const ReadPath& right) { return left.path < right.path; });
	std::vector<ReadPathInfo> readInfos;
	std::vector<PathGroup> readPaths;
	{
		std::unordered_map<std::string, size_t> nameLookup;
		for (size_t i = 0; i < rawReadPaths.size(); i++)
		{
			if (rawReadPaths[i].path.size() == 0) continue;
			assert(rawReadPaths[i].readPoses.size() > 0);
			if (readPaths.size() == 0 || readPaths.back().path != rawReadPaths[i].path)
			{
				readPaths.emplace_back();
				readPaths.back().path.insert(readPaths.back().path.end(), rawReadPaths[i].path.begin(), rawReadPaths[i].path.end());
				for (size_t j = 0; j < rawReadPaths[i].path.size(); j++)
				{
					resolvableGraph.readsCrossingNode[rawReadPaths[i].path[j].id()].emplace_back(readPaths.size()-1, j);
				}
			}
			assert(readPaths.size() > 0);
			assert(readPaths.back().path.size() > 0);
			assert(readPaths.back().path == rawReadPaths[i].path);
			if (nameLookup.count(rawReadPaths[i].readName) == 0)
			{
				size_t n = resolvableGraph.readNames.size();
				resolvableGraph.readNames.push_back(rawReadPaths[i].readName);
				nameLookup[rawReadPaths[i].readName] = n;
			}
			readInfos.emplace_back();
			readInfos.back().readPoses = rawReadPaths[i].readPoses;
			readInfos.back().expandedReadPosStart = rawReadPaths[i].expandedReadPosStart;
			readInfos.back().expandedReadPosEnd = rawReadPaths[i].expandedReadPosEnd;
			readInfos.back().readLength = rawReadPaths[i].readLength;
			readInfos.back().readLengthHPC = rawReadPaths[i].readLengthHPC;
			readPaths.back().reads.emplace_back();
			readPaths.back().reads.back().readInfoIndex = readInfos.size()-1;
			readPaths.back().reads.back().readPosZeroOffset = rawReadPaths[i].readPoses[0];
			readPaths.back().reads.back().readNameIndex = nameLookup[rawReadPaths[i].readName];
			readPaths.back().reads.back().readPosStartIndex = 0;
			readPaths.back().reads.back().readPosEndIndex = rawReadPaths[i].readPoses.size();
			readPaths.back().reads.back().leftClip = rawReadPaths[i].leftClip;
			readPaths.back().reads.back().rightClip = rawReadPaths[i].rightClip;
			assert(getNumberOfHashes(resolvableGraph, 0, 0, readPaths.back().path) == (readPaths.back().reads.back().readPosEndIndex - readPaths.back().reads.back().readPosStartIndex) + readPaths.back().reads.back().leftClip + readPaths.back().reads.back().rightClip);
		}
	}
	std::cerr << rawReadPaths.size() << " raw read paths" << std::endl;
	std::cerr << readPaths.size() << " read paths" << std::endl;
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		assert(resolvableGraph.readsCrossingNode[i].size() >= 1);
	}
	unitigifyAll(resolvableGraph, readPaths);
	auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, phmap::flat_hash_set<size_t> {});
	if (removed.nodesRemoved > 0)
	{
		std::cerr << "removed " << removed.nodesRemoved << " tips" << std::endl;
		unitigifyAll(resolvableGraph, readPaths);
	}
	auto removedEdges = removeLowCoverageCrosslinks(resolvableGraph, readPaths, 2, 10);
	if (removedEdges.edgesRemoved > 0)
	{
		std::cerr << "removed " << removedEdges.edgesRemoved << " crosslinks" << std::endl;
		unitigifyAll(resolvableGraph, readPaths);
	}
	if (guesswork)
	{
		auto removed2 = cleanComponentsByCopynumber(resolvableGraph, readPaths, 50000, 0, maxResolveLength);
		if (removed2.nodesRemoved > 0 || removed2.edgesRemoved > 0)
		{
			std::cerr << "removed " << removed2.nodesRemoved << " nodes and " << removed2.edgesRemoved << " edges" << std::endl;
			unitigifyAll(resolvableGraph, readPaths);
		}
	}
	resolveRound(resolvableGraph, readPaths, hashlist, minCoverage, maxResolveLength, maxUnconditionalResolveLength, guesswork, copycountFilterHeuristic, onlyLocalResolve);
	resolveRound(resolvableGraph, readPaths, hashlist, 1, maxResolveLength, maxUnconditionalResolveLength, guesswork, copycountFilterHeuristic, onlyLocalResolve);
	checkValidity(resolvableGraph, readPaths);
	return resolvableToUnitigs(resolvableGraph, readPaths, readInfos);
}
