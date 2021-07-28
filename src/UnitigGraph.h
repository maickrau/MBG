#ifndef UnitigGraph_h
#define UnitigGraph_h

#include <unordered_set>
#include <vector>
#include <tuple>
#include "VectorWithDirection.h"
#include "HashList.h"
#include "RankBitvector.h"

class UnitigGraph
{
public:
	std::vector<std::vector<std::pair<NodeType, bool>>> unitigs;
	std::vector<std::vector<size_t>> unitigCoverage;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCov;
	size_t edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw) const;
	size_t edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t& edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw);
	size_t& edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	double averageCoverage(size_t i) const;
	UnitigGraph filterNodes(const RankBitvector& kept) const;
	size_t numNodes() const;
	size_t numEdges() const;
	UnitigGraph filterUnitigsByCoverage(const double filter);
};

#endif
