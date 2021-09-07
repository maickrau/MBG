#include <iostream>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <phmap.h>
#include "MBGCommon.h"
#include "UnitigResolver.h"

class ResolvableUnitigGraph
{
public:
	std::vector<std::vector<std::pair<size_t, bool>>> unitigs;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> overlaps;
	std::vector<bool> unitigRemoved;
	std::vector<std::unordered_set<size_t>> readsCrossingNode;
private:
};

ResolvableUnitigGraph getUnitigs(const UnitigGraph& initial, size_t minCoverage)
{
	ResolvableUnitigGraph result;
	result.unitigs.resize(initial.unitigs.size());
	result.unitigRemoved.resize(initial.unitigs.size(), false);
	result.edges.resize(initial.unitigs.size());
	result.readsCrossingNode.resize(initial.unitigs.size());
	for (size_t i = 0; i < initial.unitigs.size(); i++)
	{
		result.unitigs[i].insert(result.unitigs[i].end(), initial.unitigs[i].begin(), initial.unitigs[i].end());
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : initial.edgeCov[fw])
		{
			if (pair.second < minCoverage) continue;
			auto canonpair = canon(fw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[fw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(fw));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : initial.edgeCov[bw])
		{
			if (pair.second < minCoverage) continue;
			auto canonpair = canon(bw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[bw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(bw));
		}
	}
	return result;
}

std::vector<std::vector<std::pair<size_t, bool>>> getUnitigPaths(const ResolvableUnitigGraph& graph, const HashList& hashlist, const std::vector<std::vector<HashType>>& kmerPaths)
{
	size_t maxKmer = 0;
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			maxKmer = std::max(maxKmer, graph.unitigs[i][j].first);
		}
	}
	std::vector<std::tuple<size_t, size_t, bool>> kmerLocator;
	kmerLocator.resize(maxKmer+1, std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true));
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			assert(std::get<0>(kmerLocator[graph.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			if (graph.unitigs[i][j].second)
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, j, true);
			}
			else
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, graph.unitigs[i].size()-1-j, false);
			}
		}
	}
	std::vector<std::vector<std::pair<size_t, bool>>> result;
	for (size_t i = 0; i < kmerPaths.size(); i++)
	{
		std::vector<std::pair<size_t, bool>> current;
		std::tuple<size_t, size_t, bool> lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
		for (size_t j = 0; j < kmerPaths[i].size(); j++)
		{
			std::pair<size_t, bool> kmer = hashlist.getNodeOrNull(kmerPaths[i][j]);
			if (kmer.first == std::numeric_limits<size_t>::max()) continue;
			if (kmer.first >= maxKmer || std::get<0>(kmerLocator[kmer.first]) == std::numeric_limits<size_t>::max())
			{
				if (current.size() > 0) result.push_back(current);
				current.clear();
				lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
				continue;
			}
			auto pos = kmerLocator[kmer.first];
			assert(std::get<0>(pos) < graph.unitigs.size());
			if (!kmer.second)
			{
				std::get<2>(pos) = !std::get<2>(pos);
				std::get<1>(pos) = graph.unitigs[std::get<0>(pos)].size()-1-std::get<1>(pos);
			}
			if (std::get<0>(lastPos) == std::numeric_limits<size_t>::max())
			{
				current.emplace_back(std::get<0>(pos), std::get<2>(pos));
				lastPos = pos;
				continue;
			}
			if (std::get<0>(pos) == std::get<0>(lastPos) && std::get<2>(pos) == std::get<2>(lastPos) && std::get<1>(pos) == std::get<1>(lastPos)+1)
			{
				lastPos = pos;
				continue;
			}
			std::pair<size_t, bool> fromEdge { std::get<0>(lastPos), std::get<2>(lastPos) };
			std::pair<size_t, bool> toEdge { std::get<0>(pos), std::get<2>(pos) };
			assert(graph.edges[fromEdge].count(toEdge) == graph.edges[reverse(toEdge)].count(reverse(fromEdge)));
			if (graph.edges[fromEdge].count(toEdge) == 1)
			{
				current.emplace_back(std::get<0>(pos), std::get<2>(pos));
				lastPos = pos;
				continue;
			}
			if (current.size() > 0) result.push_back(current);
			current.clear();
			current.emplace_back(std::get<0>(pos), std::get<2>(pos));
			lastPos = pos;
			continue;
		}
		if (current.size() > 0) result.push_back(current);
	}
	return result;
}

UnitigGraph resolvableToUnitigs(const ResolvableUnitigGraph& resolvableGraph)
{
	UnitigGraph result;
	RankBitvector newIndex { resolvableGraph.unitigs.size() };
	assert(resolvableGraph.unitigs.size() == resolvableGraph.unitigRemoved.size());
	assert(resolvableGraph.unitigs.size() == resolvableGraph.edges.size());
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		newIndex.set(i, !resolvableGraph.unitigRemoved[i]);
	}
	newIndex.buildRanks();
	size_t newSize = newIndex.getRank(newIndex.size()-1) + (newIndex.get(newIndex.size()-1) ? 1 : 0);
	result.unitigs.resize(newSize);
	result.unitigCoverage.resize(newSize);
	result.edges.resize(newSize);
	result.edgeCov.resize(newSize);
	result.edgeOvlp.resize(newSize);
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (!newIndex.get(i)) continue;
		size_t unitig = newIndex.getRank(i);
		result.unitigs[unitig].insert(result.unitigs[unitig].end(), resolvableGraph.unitigs[i].begin(), resolvableGraph.unitigs[i].end());
		for (size_t j = 0; j < resolvableGraph.unitigs[i].size(); j++)
		{
			result.unitigCoverage[unitig].push_back(0);
		}
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> newfw { newIndex.getRank(i), true };
		for (auto edge : resolvableGraph.edges[fw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			result.edges[newfw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newfw));
			result.edgeCoverage(newfw, newEdge) = 0;
			result.edgeOverlap(newfw, newEdge) = resolvableGraph.overlaps.at(canon(fw, edge));
		}
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newbw { newIndex.getRank(i), false };
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			result.edges[newbw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newbw));
			result.edgeCoverage(newbw, newEdge) = 0;
			result.edgeOverlap(newbw, newEdge) = resolvableGraph.overlaps.at(canon(bw, edge));
		}
	}
	return result;
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

void erasePath(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, size_t i)
{
	std::unordered_set<size_t> nodes;
	for (auto pos : readPaths[i])
	{
		nodes.emplace(pos.first);
	}
	for (auto node : nodes)
	{
		assert(resolvableGraph.readsCrossingNode[node].count(i) == 1);
		resolvableGraph.readsCrossingNode[node].erase(i);
	}
	readPaths[i].clear();
}

void addPath(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const std::vector<std::pair<size_t, bool>>& newPath)
{
	for (size_t i = 0; i < newPath.size(); i++)
	{
		assert(!resolvableGraph.unitigRemoved[newPath[i].first]);
		if (i > 0)
		{
			assert(resolvableGraph.edges[newPath[i-1]].count(newPath[i]) == 1);
			assert(resolvableGraph.edges[reverse(newPath[i])].count(reverse(newPath[i-1])) == 1);
		}
	}
	size_t newIndex = readPaths.size();
	readPaths.emplace_back(newPath);
	std::unordered_set<size_t> nodes;
	for (auto pos : newPath)
	{
		nodes.insert(pos.first);
	}
	for (auto node : nodes)
	{
		assert(resolvableGraph.readsCrossingNode[node].count(newIndex) == 0);
		resolvableGraph.readsCrossingNode[node].emplace(newIndex);
	}
}

void replacePathNodes(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const std::vector<std::pair<size_t, bool>>& newUnitig, size_t newUnitigIndex)
{
	std::unordered_set<size_t> relevantReads;
	for (auto pos : newUnitig)
	{
		relevantReads.insert(resolvableGraph.readsCrossingNode[pos.first].begin(), resolvableGraph.readsCrossingNode[pos.first].end());
	}
	std::unordered_map<size_t, bool> nodeForwardInUnitig;
	std::unordered_set<size_t> nodeInUnitig;
	for (auto pos : newUnitig)
	{
		assert(nodeInUnitig.count(pos.first) == 0);
		nodeInUnitig.emplace(pos.first);
		nodeForwardInUnitig[pos.first] = pos.second;
	}
	for (size_t i : relevantReads)
	{
		std::vector<std::pair<size_t, bool>> newPath;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (nodeInUnitig.count(readPaths[i][j].first) == 0)
			{
				newPath.emplace_back(readPaths[i][j]);
				continue;
			}
			if (j == 0)
			{
				bool fw = nodeForwardInUnitig.at(readPaths[i][j].first);
				if (!readPaths[i][j].second) fw = !fw;
				newPath.emplace_back(newUnitigIndex, fw);
			}
			else if (readPaths[i][j] == newUnitig[0])
			{
				newPath.emplace_back(newUnitigIndex, true);
			}
			else if (readPaths[i][j] == reverse(newUnitig.back()))
			{
				newPath.emplace_back(newUnitigIndex, false);
			}
		}
		addPath(resolvableGraph, readPaths, newPath);
		erasePath(resolvableGraph, readPaths, i);
	}
}

size_t unitigifyOne(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const size_t unitig)
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
	if (newUnitig.size() == 1) return 1;
	size_t newIndex = resolvableGraph.unitigs.size();
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
	std::unordered_set<size_t> nodesInUnitig;
	for (auto pos : newUnitig)
	{
		assert(nodesInUnitig.count(pos.first) == 0);
		nodesInUnitig.insert(pos.first);
	}
	std::pair<size_t, bool> bw { newIndex, false };
	std::pair<size_t, bool> fw { newIndex, true };
	for (auto edge : resolvableGraph.edges[reverse(newUnitig[0])])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(bw, bw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
		else if (edge == newUnitig[0])
		{
			resolvableGraph.overlaps[canon(bw, fw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(fw);
			resolvableGraph.edges[reverse(fw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(bw, edge)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
	}
	for (auto edge : resolvableGraph.edges[newUnitig.back()])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(fw, bw)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(newUnitig.back()));
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(fw, edge)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(newUnitig.back()));
		}
	}
	for (size_t i = 0; i < newUnitig.size(); i++)
	{
		resolvableGraph.unitigRemoved[newUnitig[i].first] = true;
		resolvableGraph.edges[newUnitig[i]].clear();
		resolvableGraph.edges[reverse(newUnitig[i])].clear();
	}
	replacePathNodes(resolvableGraph, readPaths, newUnitig, newIndex);
	return newUnitig.size();
}

void unitigifyAll(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths)
{
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		unitigifyOne(resolvableGraph, readPaths, i);
	}
}

void unresolveRecursively(const ResolvableUnitigGraph& resolvableGraph, const std::unordered_set<size_t>& resolvables, std::unordered_set<size_t>& unresolvables, const size_t node)
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

void createEdgeNode(ResolvableUnitigGraph& resolvableGraph, std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes, const std::unordered_set<size_t>& resolvables, const std::unordered_set<size_t>& unresolvables, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	size_t newIndex = resolvableGraph.unitigs.size();
	newEdgeNodes[std::make_pair(from, to)] = newIndex;
	resolvableGraph.unitigs.emplace_back();
	std::vector<std::pair<size_t, bool>> add = resolvableGraph.unitigs[from.first];
	if (!from.second) add = revCompPath(add);
	resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin(), add.end());
	size_t start = resolvableGraph.overlaps.at(canon(from, to));
	assert(start < resolvableGraph.unitigs[to.first].size());
	add = resolvableGraph.unitigs[to.first];
	if (!to.second) add = revCompPath(add);
	for (size_t i = 0; i < start; i++)
	{
		assert(resolvableGraph.unitigs.back()[resolvableGraph.unitigs.back().size() - start+i] == add[i]);
	}
	assert(add.size() > start);
	if (resolvables.count(to.first) == 1 && unresolvables.count(to.first) == 0)
	{
		resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin()+start, add.end());
	}
	else
	{
		resolvableGraph.unitigs.back().emplace_back(add[start]);
	}
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	if (resolvables.count(to.first) == 0 || unresolvables.count(to.first) == 1)
	{
		resolvableGraph.edges[std::make_pair(newIndex, true)].emplace(to);
		resolvableGraph.edges[reverse(to)].emplace(std::make_pair(newIndex, false));
		assert(resolvableGraph.unitigs[newIndex].size() > resolvableGraph.overlaps.at(canon(from, to)) + 1);
		assert(resolvableGraph.unitigs[to.first].size() > resolvableGraph.overlaps.at(canon(from, to)) + 1);
		resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = resolvableGraph.overlaps.at(canon(from, to)) + 1;
		// resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = resolvableGraph.unitigs[from.first].size();
	}
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getValidTriplets(const ResolvableUnitigGraph& resolvableGraph, const std::unordered_set<size_t>& resolvables, const std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, size_t node, size_t minCoverage)
{
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> empty;
	std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> tripletCoverage;
	for (size_t i : resolvableGraph.readsCrossingNode[node])
	{
		for (size_t j = 1; j < readPaths[i].size()-1; j++)
		{
			if (readPaths[i][j].first != node) continue;
			std::pair<size_t, bool> left;
			std::pair<size_t, bool> right;
			if (readPaths[i][j].second)
			{
				left = readPaths[i][j-1];
				right = readPaths[i][j+1];
			}
			else
			{
				left = reverse(readPaths[i][j+1]);
				right = reverse(readPaths[i][j-1]);
			}
			tripletCoverage[std::make_pair(left, right)] += 1;
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
	std::unordered_set<std::pair<size_t, bool>> coveredInNeighbors;
	std::unordered_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		coveredInNeighbors.emplace(reverse(pair.first));
		coveredOutNeighbors.emplace(pair.second);
	}
	assert(coveredInNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, false)].size());
	assert(coveredOutNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, true)].size());
	if (coveredInNeighbors.size() < resolvableGraph.edges[std::make_pair(node, false)].size()) return empty;
	if (coveredOutNeighbors.size() < resolvableGraph.edges[std::make_pair(node, true)].size()) return empty;
	return coveredTriplets;
}

void replacePaths(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const std::unordered_set<size_t>& resolvables, const std::unordered_set<size_t>& unresolvables, const std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes)
{
	std::unordered_set<size_t> relevantReads;
	std::unordered_set<size_t> nodesInNewEdgeNodes;
	for (auto pair : newEdgeNodes)
	{
		assert(nodesInNewEdgeNodes.count(pair.second) == 0);
		assert(resolvables.count(pair.first.first.first) == 1);
		assert(unresolvables.count(pair.first.first.first) == 0);
		nodesInNewEdgeNodes.insert(pair.second);
	}
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		relevantReads.insert(resolvableGraph.readsCrossingNode[node].begin(), resolvableGraph.readsCrossingNode[node].end());
	}
	for (auto i : relevantReads)
	{
		std::vector<std::pair<size_t, bool>> newPath;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (resolvables.count(readPaths[i][j].first) == 0 || unresolvables.count(readPaths[i][j].first) == 1)
			{
				newPath.push_back(readPaths[i][j]);
				continue;
			}
			if (j > 0 && newEdgeNodes.count(std::make_pair(reverse(readPaths[i][j]), reverse(readPaths[i][j-1]))) == 1)
			{
				newPath.emplace_back(newEdgeNodes.at(std::make_pair(reverse(readPaths[i][j]), reverse(readPaths[i][j-1]))), false);
			}
			if (j < readPaths[i].size()-1 && newEdgeNodes.count(std::make_pair(readPaths[i][j], readPaths[i][j+1])) == 1)
			{
				newPath.emplace_back(newEdgeNodes.at(std::make_pair(readPaths[i][j], readPaths[i][j+1])), true);
			}
		}
		size_t lastStart = 0;
		for (size_t j = 1; j < newPath.size(); j++)
		{
			assert(resolvableGraph.edges[newPath[j-1]].count(newPath[j]) == resolvableGraph.edges[reverse(newPath[j])].count(reverse(newPath[j-1])));
			if (resolvableGraph.edges[newPath[j-1]].count(newPath[j]) == 0)
			{
				// assert(nodesInNewEdgeNodes.count(newPath[j-1].first) == 1);
				// assert(nodesInNewEdgeNodes.count(newPath[j].first) == 1);
				addPath(resolvableGraph, readPaths, std::vector<std::pair<size_t, bool>> { newPath.begin() + lastStart, newPath.begin() + j });
				lastStart = j;
			}
		}
		addPath(resolvableGraph, readPaths, std::vector<std::pair<size_t, bool>> { newPath.begin() + lastStart, newPath.end() });
		erasePath(resolvableGraph, readPaths, i);
	}
}

struct ResolutionResult
{
public:
	size_t nodesResolved;
	size_t nodesAdded;
	std::unordered_set<size_t> maybeUnitigifiable;
};

ResolutionResult resolve(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const std::unordered_set<size_t>& resolvables, const size_t minCoverage)
{
	ResolutionResult result;
	result.nodesResolved = 0;
	result.nodesAdded = 0;
	std::unordered_set<size_t> unresolvables;
	for (auto node : resolvables)
	{
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
		if (triplets.size() == 0)
		{
			unresolvables.insert(node);
			// unresolveRecursively(resolvableGraph, resolvables, unresolvables, node);
		}
	}
	while (true)
	{
		bool unmadeAny = false;
		for (auto node : resolvables)
		{
			if (unresolvables.count(node) == 1) continue;
			auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
			for (auto triplet : triplets)
			{
				if (resolvables.count(triplet.first.first) == 0 || unresolvables.count(triplet.first.first) == 1)
				{
					if (resolvables.count(triplet.second.first) == 0 || unresolvables.count(triplet.second.first) == 1)
					{
						if (resolvableGraph.unitigs[triplet.first.first].size() == resolvableGraph.overlaps.at(canon(triplet.first, std::make_pair(node, true)))+1)
						{
							if (resolvableGraph.unitigs[triplet.second.first].size() == resolvableGraph.overlaps.at(canon(std::make_pair(node, true), triplet.second))+1)
							{
								assert(unresolvables.count(node) == 0);
								unmadeAny = true;
								unresolvables.insert(node);
								break;
							}
						}
					}
				}
			}
		}
		if (!unmadeAny) break;
	}
	if (unresolvables.size() == resolvables.size()) return result;
	std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newEdgeNodes;
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		assert(getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage).size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto edge : resolvableGraph.edges[pos])
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigs[edge.first].size() == resolvableGraph.overlaps.at(canon(pos, edge)) + 1) continue;
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
		}
		pos = std::make_pair(node, false);
		for (auto edge : resolvableGraph.edges[pos])
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigs[edge.first].size() == resolvableGraph.overlaps.at(canon(pos, edge)) + 1) continue;
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
		}
	}
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
		assert(triplets.size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto triplet : triplets)
		{
			const std::pair<size_t, bool> before = triplet.first;
			const std::pair<size_t, bool> after = triplet.second;
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
				assert(resolvableGraph.unitigs[before.first].size() == resolvableGraph.overlaps.at(canon(before, pos))+1);
				assert(resolvables.count(before.first) == 0 || unresolvables.count(before.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(pos, after)) == 1 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1);
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
				assert(resolvableGraph.unitigs[after.first].size() == resolvableGraph.overlaps.at(canon(pos, after))+1);
				assert(resolvables.count(after.first) == 0 || unresolvables.count(after.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1 || newEdgeNodes.count(std::make_pair(before, pos)) == 1);
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
		}
	}
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
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
	replacePaths(resolvableGraph, readPaths, resolvables, unresolvables, newEdgeNodes);
	result.nodesResolved = resolvables.size() - unresolvables.size();
	result.nodesAdded = newEdgeNodes.size();
	return result;
}

class KmerLengthCompare
{
public:
	KmerLengthCompare(const ResolvableUnitigGraph& graph) :
	graph(graph)
	{
	}
	bool operator()(size_t left, size_t right) const
	{
		return graph.unitigs[left].size() > graph.unitigs[right].size();
	}
private:
	const ResolvableUnitigGraph& graph;
};

void checkValidity(const ResolvableUnitigGraph& graph, const std::vector<std::vector<std::pair<size_t, bool>>>& readPaths)
{
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		if (graph.unitigRemoved[i])
		{
			assert(graph.edges[std::make_pair(i, true)].size() == 0);
			assert(graph.edges[std::make_pair(i, false)].size() == 0);
			continue;
		}
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges[fw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(fw)) == 1);
			assert(graph.edges[fw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2);
		}
		for (auto edge : graph.edges[bw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(bw)) == 1);
			assert(graph.edges[bw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2);
		}
	}
	for (const auto& path : readPaths)
	{
		for (size_t i = 0; i < path.size(); i++)
		{
			assert(!graph.unitigRemoved[path[i].first]);
			if (i > 0) assert(graph.edges[path[i-1]].count(path[i]) == 1);
		}
	}
}

void resolveRound(ResolvableUnitigGraph& resolvableGraph, std::vector<std::vector<std::pair<size_t, bool>>>& readPaths, const size_t minCoverage)
{
	checkValidity(resolvableGraph, readPaths);
	std::priority_queue<size_t, std::vector<size_t>, KmerLengthCompare> queue { KmerLengthCompare { resolvableGraph } };
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		queue.emplace(i);
	}
	size_t lastTopSize = 0;
	while (queue.size() > 0)
	{
		if (lastTopSize >= 1000) break;
		size_t topSize = resolvableGraph.unitigs[queue.top()].size();
		assert(topSize >= lastTopSize);
		lastTopSize = topSize;
		std::unordered_set<size_t> resolvables;
		while (queue.size() > 0 && resolvableGraph.unitigs[queue.top()].size() == topSize)
		{
			if (!resolvableGraph.unitigRemoved[queue.top()])
			{
				if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 1 && resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 1)
				{
					if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 2 || resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 2)
					{
						resolvables.emplace(queue.top());
					}
					else
					{
						assert(resolvableGraph.edges[std::make_pair(queue.top(), true)].size() == 1);
						assert(resolvableGraph.edges[std::make_pair(queue.top(), false)].size() == 1);
						if (resolvableGraph.unitigs[resolvableGraph.edges[std::make_pair(queue.top(), true)].begin()->first].size() == topSize || resolvableGraph.unitigs[resolvableGraph.edges[std::make_pair(queue.top(), false)].begin()->first].size() == topSize)
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
		checkValidity(resolvableGraph, readPaths);
		auto resolutionResult = resolve(resolvableGraph, readPaths, resolvables, minCoverage);
		size_t newSize = resolvableGraph.unitigs.size();
		std::cerr << "try resolve l=" << topSize << ", replaced " << resolutionResult.nodesResolved << " nodes with " << resolutionResult.nodesAdded << " nodes";
		size_t unitigified = 0;
		size_t unitigifiedTo = 0;
		for (auto i : resolutionResult.maybeUnitigifiable)
		// for (size_t i = oldSize; i < newSize; i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			size_t unitigifiedHere = unitigifyOne(resolvableGraph, readPaths, i);
			if (unitigifiedHere > 1)
			{
				unitigified += unitigifiedHere;
				unitigifiedTo += 1;
			}
		}
		if (unitigified > 0)
		{
			std::cerr << ", unitigified " << unitigified << " nodes to " << unitigifiedTo << " nodes";
		}
		std::cerr << std::endl;
		checkValidity(resolvableGraph, readPaths);
		for (size_t i = oldSize; i < resolvableGraph.unitigs.size(); i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			queue.emplace(i);
		}
	}
}

UnitigGraph resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, const std::vector<std::vector<HashType>>& paths, const size_t minCoverage)
{
	auto resolvableGraph = getUnitigs(initial, minCoverage);
	auto readPaths = getUnitigPaths(resolvableGraph, hashlist, paths);
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (auto pos : readPaths[i])
		{
			resolvableGraph.readsCrossingNode[pos.first].emplace(i);
		}
	}
	unitigifyAll(resolvableGraph, readPaths);
	resolveRound(resolvableGraph, readPaths, minCoverage);
	// resolveRound(resolvableGraph, readPaths, minCoverage);
	return resolvableToUnitigs(resolvableGraph);
}
