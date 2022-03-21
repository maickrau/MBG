#include <unordered_set>
#include <cassert>
#include <phmap.h>
#include <iostream>
#include "MBGCommon.h"
#include "UnitigGraph.h"

size_t UnitigGraph::edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw) const
{
	return edgeCoverage(std::make_pair(from, fromFw), std::make_pair(to, toFw));
}
size_t UnitigGraph::edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return edgeCov[from].at(to);
}
size_t& UnitigGraph::edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw)
{
	return edgeCoverage(std::make_pair(from, fromFw), std::make_pair(to, toFw));
}
size_t& UnitigGraph::edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	std::tie(from, to) = canon(from, to);
	return edgeCov[from][to];
}
size_t UnitigGraph::edgeOverlap(size_t from, bool fromFw, size_t to, bool toFw) const
{
	return edgeOverlap(std::make_pair(from, fromFw), std::make_pair(to, toFw));
}
size_t UnitigGraph::edgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
{
	std::tie(from, to) = canon(from, to);
	return edgeOvlp[from].at(to);
}
size_t& UnitigGraph::edgeOverlap(size_t from, bool fromFw, size_t to, bool toFw)
{
	return edgeOverlap(std::make_pair(from, fromFw), std::make_pair(to, toFw));
}
size_t& UnitigGraph::edgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	std::tie(from, to) = canon(from, to);
	return edgeOvlp[from][to];
}
double UnitigGraph::averageCoverage(size_t i) const
{
	size_t total = 0;
	for (auto cov : unitigCoverage[i]) total += cov;
	assert(unitigCoverage[i].size() > 0);
	return (double)total / (double)unitigCoverage[i].size();
}
UnitigGraph UnitigGraph::filterNodes(const RankBitvector& kept) const
{
	if (kept.size() == 0) return *this;
	assert(kept.size() == unitigs.size());
	UnitigGraph result;
	size_t newSize = kept.getRank(kept.size()-1) + (kept.get(kept.size()-1) ? 1 : 0);
	if (newSize == unitigs.size()) return *this;
	result.unitigs.resize(newSize);
	result.unitigCoverage.resize(newSize);
	result.edges.resize(newSize);
	result.edgeCov.resize(newSize);
	result.edgeOvlp.resize(newSize);
	result.leftClip.resize(newSize);
	result.rightClip.resize(newSize);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!kept.get(i)) continue;
		size_t newIndex = kept.getRank(i);
		result.unitigs[newIndex] = unitigs[i];
		result.leftClip[newIndex] = leftClip[i];
		result.rightClip[newIndex] = rightClip[i];
		result.unitigCoverage[newIndex] = unitigCoverage[i];
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newFw { newIndex, true };
		std::pair<size_t, bool> newBw { newIndex, false };
		for (auto to : edges[fw])
		{
			if (!kept.get(to.first)) continue;
			result.edges[newFw].emplace(kept.getRank(to.first), to.second);
		}
		for (auto to : edges[bw])
		{
			if (!kept.get(to.first)) continue;
			result.edges[newBw].emplace(kept.getRank(to.first), to.second);
		}
		for (auto pair : edgeCov[fw])
		{
			if (!kept.get(pair.first.first)) continue;
			result.edgeCov[newFw][std::make_pair(kept.getRank(pair.first.first), pair.first.second)] = pair.second;
		}
		for (auto pair : edgeCov[bw])
		{
			if (!kept.get(pair.first.first)) continue;
			result.edgeCov[newBw][std::make_pair(kept.getRank(pair.first.first), pair.first.second)] = pair.second;
		}
		for (auto pair : edgeOvlp[fw])
		{
			if (!kept.get(pair.first.first)) continue;
			result.edgeOvlp[newFw][std::make_pair(kept.getRank(pair.first.first), pair.first.second)] = pair.second;
		}
		for (auto pair : edgeOvlp[bw])
		{
			if (!kept.get(pair.first.first)) continue;
			result.edgeOvlp[newBw][std::make_pair(kept.getRank(pair.first.first), pair.first.second)] = pair.second;
		}
	}
	return result;
}
size_t UnitigGraph::numNodes() const
{
	return unitigs.size();
}
size_t UnitigGraph::numEdges() const
{
	size_t result = 0;
	for (size_t i = 0; i < edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : edges[fw])
		{
			auto c = canon(fw, edge);
			if (c.first == fw && c.second == edge) result += 1;
		}
		for (auto edge : edges[bw])
		{
			auto c = canon(bw, edge);
			if (c.first == bw && c.second == edge) result += 1;
		}
	}
	return result;
}
std::unordered_set<std::pair<size_t, bool>> UnitigGraph::findReachableNewTips(const RankBitvector& kept, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start) const
{
	assert(kept.get(start.first));
	std::vector<std::pair<size_t, bool>> stack;
	stack.push_back(start);
	std::unordered_set<std::pair<size_t, bool>> visited;
	std::unordered_set<std::pair<size_t, bool>> result;
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		if (kept.get(top.first) && top != start)
		{
			if (newlyTip[reverse(top)])
			{
				result.insert(top);
			}
			continue;
		}
		visited.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	return result;
}
void UnitigGraph::keepReachableNewTips(const RankBitvector& kept, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start, const std::unordered_set<std::pair<size_t, bool>>& reachableTips, std::unordered_set<size_t>& newlyKept) const
{
	std::unordered_set<std::pair<size_t, bool>> reachableFw;
	std::unordered_set<std::pair<size_t, bool>> reachableBw;
	std::vector<std::pair<size_t, bool>> stack;
	stack.push_back(start);
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (reachableFw.count(top) == 1) continue;
		if (kept.get(top.first) && top != start) continue;
		reachableFw.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	for (auto node : reachableTips)
	{
		stack.push_back(reverse(node));
	}
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (reachableBw.count(top) == 1) continue;
		if (kept.get(top.first) && reachableTips.count(reverse(top)) == 0) continue;
		reachableBw.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	for (auto node : reachableFw)
	{
		if (reachableBw.count(reverse(node)) == 1)
		{
			newlyKept.insert(node.first);
		}
	}
}
bool UnitigGraph::isTipGap(const RankBitvector& kept, const std::pair<size_t, bool> start) const
{
	if (!kept.get(start.first)) return false;
	if (edges[start].size() == 0) return false;
	for (auto edge : edges[start])
	{
		if (kept.get(edge.first)) return false;
	}
	return true;
}
void UnitigGraph::keepTipGaps(RankBitvector& kept) const
{
	VectorWithDirection<bool> newlyTip;
	newlyTip.resize(unitigs.size(), false);
	size_t newlyTipped = 0;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!kept.get(i)) continue;
		std::pair<size_t, bool> fw { i, true };
		if (isTipGap(kept, fw))
		{
			newlyTip[fw] = true;
			newlyTipped += 1;
		}
		std::pair<size_t, bool> bw { i, false };
		if (isTipGap(kept, bw))
		{
			newlyTip[bw] = true;
			newlyTipped += 1;
		}
	}
	std::cerr << newlyTipped << " tips created" << std::endl;
	std::unordered_set<size_t> newlyKept;
	size_t keptCheck = 0;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!kept.get(i)) continue;
		std::pair<size_t, bool> fw { i, true };
		if (edges[fw].size() > 0 && newlyTip[fw])
		{
			auto reachableTips = findReachableNewTips(kept, newlyTip, fw);
			if (reachableTips.size() > 0)
			{
				keptCheck += 1;
				keepReachableNewTips(kept, newlyTip, fw, reachableTips, newlyKept);
			}
		}
		std::pair<size_t, bool> bw { i, false };
		if (edges[bw].size() > 0 && newlyTip[bw])
		{
			auto reachableTips = findReachableNewTips(kept, newlyTip, bw);
			if (reachableTips.size() > 0)
			{
				keptCheck += 1;
				keepReachableNewTips(kept, newlyTip, bw, reachableTips, newlyKept);
			}
		}
	}
	std::cerr << keptCheck << " tips checked for keeping" << std::endl;
	std::cerr << "re-inserted " << newlyKept.size() << " nodes to prevent gaps" << std::endl;
	for (auto node : newlyKept)
	{
		assert(!kept.get(node));
		kept.set(node, true);
	}
}
UnitigGraph UnitigGraph::filterUnitigsByCoverage(const double filter, const bool keepGaps)
{
	RankBitvector kept { unitigs.size() };
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		kept.set(i, averageCoverage(i) >= filter);
	}
	if (keepGaps) keepTipGaps(kept);
	kept.buildRanks();
	UnitigGraph filtered = filterNodes(kept);
	return filtered;
}
void UnitigGraph::sort(const std::vector<size_t>& kmerMapping)
{
	std::vector<bool> swapOrientation;
	swapOrientation.resize(unitigs.size(), false);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		for (size_t j = 0; j < unitigs[i].size(); j++)
		{
			unitigs[i][j].first = kmerMapping[unitigs[i][j].first];
		}
		if (unitigs[i].size() == 1)
		{
			swapOrientation[i] = !unitigs[i][0].second;
		}
		else
		{
			assert(unitigs[i].size() >= 2);
			swapOrientation[i] = (unitigs[i].back().first < unitigs[i][0].first);
		}
	}
	std::vector<size_t> unitigOrder;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		unitigOrder.push_back(i);
	}
	std::sort(unitigOrder.begin(), unitigOrder.end(), [this](size_t left, size_t right) { return std::min(unitigs[left][0].first, unitigs[left].back().first) < std::min(unitigs[right][0].first, unitigs[right].back().first); });
	std::vector<size_t> unitigMapping;
	unitigMapping.resize(unitigOrder.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		assert(unitigMapping[unitigOrder[i]] == std::numeric_limits<size_t>::max());
		unitigMapping[unitigOrder[i]] = i;
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		assert(unitigMapping[i] != std::numeric_limits<size_t>::max());
		assert(unitigMapping[i] < unitigs.size());
	}
	{
		std::vector<std::vector<std::pair<NodeType, bool>>> newUnitigs;
		newUnitigs.resize(unitigs.size());
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			if (swapOrientation[i])
			{
				std::reverse(unitigs[i].begin(), unitigs[i].end());
				for (size_t j = 0; j < unitigs[i].size(); j++)
				{
					unitigs[i][j].second = !unitigs[i][j].second;
				}
			}
			assert(unitigs[i][0].first < unitigs[i].back().first || unitigs[i].size() == 1);
			std::swap(unitigs[i], newUnitigs[unitigMapping[i]]);
		}
		std::swap(unitigs, newUnitigs);
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (swapOrientation[i])
		{
			std::swap(leftClip[i], rightClip[i]);
		}
	}
	{
		std::vector<size_t> newLeftClip;
		newLeftClip.resize(leftClip.size());
		for (size_t i = 0; i < leftClip.size(); i++)
		{
			newLeftClip[unitigMapping[i]] = leftClip[i];
		}
		std::swap(leftClip, newLeftClip);
	}
	{
		std::vector<size_t> newRightClip;
		newRightClip.resize(rightClip.size());
		for (size_t i = 0; i < rightClip.size(); i++)
		{
			newRightClip[unitigMapping[i]] = rightClip[i];
		}
		std::swap(rightClip, newRightClip);
	}
	{
		std::vector<std::vector<size_t>> newUnitigCoverage;
		newUnitigCoverage.resize(unitigCoverage.size());
		for (size_t i = 0; i < unitigCoverage.size(); i++)
		{
			if (swapOrientation[i])
			{
				std::reverse(unitigCoverage[i].begin(), unitigCoverage[i].end());
			}
			std::swap(newUnitigCoverage[unitigMapping[i]], unitigCoverage[i]);
		}
		std::swap(unitigCoverage, newUnitigCoverage);
	}
	{
		VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> newEdges;
		newEdges.resize(edges.size());
		for (size_t i = 0; i < newEdges.size(); i++)
		{
			for (auto to : edges[std::make_pair(i, true)])
			{
				newEdges[std::make_pair(unitigMapping[i], true ^ swapOrientation[i])].emplace(unitigMapping[to.first], to.second ^ swapOrientation[to.first]);
			}
			for (auto to : edges[std::make_pair(i, false)])
			{
				newEdges[std::make_pair(unitigMapping[i], false ^ swapOrientation[i])].emplace(unitigMapping[to.first], to.second ^ swapOrientation[to.first]);
			}
		}
		std::swap(newEdges, edges);
	}
	{
		VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> newEdgeCov;
		newEdgeCov.resize(edgeCov.size());
		for (size_t i = 0; i < newEdgeCov.size(); i++)
		{
			std::pair<size_t, bool> from { unitigMapping[i], true ^ swapOrientation[i] };
			for (auto pair : edgeCov[std::make_pair(i, true)])
			{
				std::pair<size_t, bool> to { std::make_pair(unitigMapping[pair.first.first], pair.first.second ^ swapOrientation[pair.first.first]) };
				auto key = canon(from, to);
				newEdgeCov[key.first].emplace(key.second, pair.second);
			}
			from.second = !from.second;
			for (auto pair : edgeCov[std::make_pair(i, false)])
			{
				std::pair<size_t, bool> to { std::make_pair(unitigMapping[pair.first.first], pair.first.second ^ swapOrientation[pair.first.first]) };
				auto key = canon(from, to);
				newEdgeCov[key.first].emplace(key.second, pair.second);
			}
		}
		std::swap(newEdgeCov, edgeCov);
	}
	{
		VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> newEdgeOvlp;
		newEdgeOvlp.resize(edgeOvlp.size());
		for (size_t i = 0; i < newEdgeOvlp.size(); i++)
		{
			std::pair<size_t, bool> from { unitigMapping[i], true ^ swapOrientation[i] };
			for (auto pair : edgeOvlp[std::make_pair(i, true)])
			{
				std::pair<size_t, bool> to { std::make_pair(unitigMapping[pair.first.first], pair.first.second ^ swapOrientation[pair.first.first]) };
				auto key = canon(from, to);
				newEdgeOvlp[key.first].emplace(key.second, pair.second);
			}
			from.second = !from.second;
			for (auto pair : edgeOvlp[std::make_pair(i, false)])
			{
				std::pair<size_t, bool> to { std::make_pair(unitigMapping[pair.first.first], pair.first.second ^ swapOrientation[pair.first.first]) };
				auto key = canon(from, to);
				newEdgeOvlp[key.first].emplace(key.second, pair.second);
			}
		}
		std::swap(newEdgeOvlp, edgeOvlp);
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		assert(unitigs[i].size() >= 1);
		assert(unitigs[i][0].first < unitigs[i].back().first || unitigs[i].size() == 1);
		assert(i == 0 || (unitigs[i][0].first > unitigs[i-1][0].first));
	}
}