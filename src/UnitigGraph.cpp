#include <unordered_set>
#include <cassert>
#include <phmap.h>
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
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!kept.get(i)) continue;
		size_t newIndex = kept.getRank(i);
		result.unitigs[newIndex] = unitigs[i];
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
UnitigGraph UnitigGraph::filterUnitigsByCoverage(const double filter)
{
	RankBitvector kept { unitigs.size() };
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		kept.set(i, averageCoverage(i) >= filter);
	}
	kept.buildRanks();
	UnitigGraph filtered = filterNodes(kept);
	return filtered;
}
