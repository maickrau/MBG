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
UnitigGraph UnitigGraph::filterNodes(const std::vector<bool>& kept) const
{
	assert(kept.size() == unitigs.size());
	UnitigGraph result;
	std::vector<size_t> newIndex;
	newIndex.resize(unitigs.size(), std::numeric_limits<size_t>::max());
	size_t nextIndex = 0;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!kept[i]) continue;
		newIndex[i] = nextIndex;
		nextIndex += 1;
	}
	result.unitigs.resize(nextIndex);
	result.unitigCoverage.resize(nextIndex);
	result.edges.resize(nextIndex);
	result.edgeCov.resize(nextIndex);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (newIndex[i] == std::numeric_limits<size_t>::max()) continue;
		result.unitigs[newIndex[i]] = unitigs[i];
		result.unitigCoverage[newIndex[i]] = unitigCoverage[i];
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newFw { newIndex[i], true };
		std::pair<size_t, bool> newBw { newIndex[i], false };
		for (auto to : edges[fw])
		{
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			result.edges[newFw].emplace(newIndex[to.first], to.second);
		}
		for (auto to : edges[bw])
		{
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			result.edges[newBw].emplace(newIndex[to.first], to.second);
		}
		for (auto pair : edgeCov[fw])
		{
			if (newIndex[pair.first.first] == std::numeric_limits<size_t>::max()) continue;
			result.edgeCov[newFw][std::make_pair(newIndex[pair.first.first], pair.first.second)] = pair.second;
		}
		for (auto pair : edgeCov[bw])
		{
			if (newIndex[pair.first.first] == std::numeric_limits<size_t>::max()) continue;
			result.edgeCov[newBw][std::make_pair(newIndex[pair.first.first], pair.first.second)] = pair.second;
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
	std::vector<bool> kept;
	kept.resize(unitigs.size(), true);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (averageCoverage(i) < filter) kept[i] = false;
	}
	UnitigGraph filtered = filterNodes(kept);
	return filtered;
}

std::pair<std::string, std::vector<uint16_t>> UnitigGraph::getSequenceAndLength(size_t unitig, const HashList& hashlist) const
{
	std::string sequenceRLE;
	std::vector<uint16_t> sequenceCharacterLength;
	for (size_t j = 0; j < unitigs[unitig].size(); j++)
	{
		auto to = unitigs[unitig][j];
		std::string sequenceRLEHere;
		std::vector<uint16_t> sequenceCharacterLengthHere = hashlist.getHashCharacterLength(to.first);
		if (to.second)
		{
			sequenceRLEHere = hashlist.getHashSequenceRLE(to.first).toString();
		}
		else
		{
			sequenceRLEHere = revCompRLE(hashlist.getHashSequenceRLE(to.first)).toString();
			std::reverse(sequenceCharacterLengthHere.begin(), sequenceCharacterLengthHere.end());
		}
		if (j > 0)
		{
			auto from = unitigs[unitig][j-1];
			size_t overlap = hashlist.getOverlap(from, to);
			assert(overlap < sequenceRLEHere.size());
			sequenceRLEHere = sequenceRLEHere.substr(overlap);
			sequenceCharacterLengthHere.erase(sequenceCharacterLengthHere.begin(), sequenceCharacterLengthHere.begin() + overlap);
		}
		sequenceRLE.insert(sequenceRLE.end(), sequenceRLEHere.begin(), sequenceRLEHere.end());
		sequenceCharacterLength.insert(sequenceCharacterLength.end(), sequenceCharacterLengthHere.begin(), sequenceCharacterLengthHere.end());
	}
	return std::make_pair(sequenceRLE, sequenceCharacterLength);
}