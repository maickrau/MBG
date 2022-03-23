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
	std::vector<size_t> leftClip;
	std::vector<size_t> rightClip;
	std::vector<std::vector<size_t>> unitigCoverage;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCov;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeOvlp;
	size_t edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw) const;
	size_t edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t& edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw);
	size_t& edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	size_t edgeOverlap(size_t from, bool fromFw, size_t to, bool toFw) const;
	size_t edgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t& edgeOverlap(size_t from, bool fromFw, size_t to, bool toFw);
	size_t& edgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	double averageCoverage(size_t i) const;
	UnitigGraph filterNodes(const RankBitvector& kept) const;
	size_t numNodes() const;
	size_t numEdges() const;
	UnitigGraph filterUnitigsByCoverage(const double filter, const bool keepGaps);
	void sort(const std::vector<size_t>& kmerMapping);
private:
	bool isTipGap(const RankBitvector& kept, const std::pair<size_t, bool> start) const;
	void keepReachableNewTips(const RankBitvector& kept, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start, const std::unordered_set<std::pair<size_t, bool>>& reachableTips, std::unordered_set<size_t>& newlyKept) const;
	std::unordered_set<std::pair<size_t, bool>> findReachableNewTips(const RankBitvector& kept, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start) const;
	void keepTipGaps(RankBitvector& kept) const;
};

#endif
