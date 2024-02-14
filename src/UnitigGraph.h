#ifndef MBGUnitigGraph_h
#define MBGUnitigGraph_h

#include <unordered_set>
#include <vector>
#include <tuple>
#include "VectorWithDirection.h"
#include "HashList.h"
#include "RankBitvector.h"
#include "SparseEdgeContainer.h"
#include "MostlySparse2DHashmap.h"

namespace MBG
{

class UnitigGraph
{
public:
	std::vector<std::vector<std::pair<NodeType, bool>>> unitigs;
	std::vector<uint32_t> leftClip;
	std::vector<uint32_t> rightClip;
	std::vector<std::vector<size_t>> unitigCoverage;
	SparseEdgeContainer edges;
	MostlySparse2DHashmap<uint8_t, size_t> edgeCov;
	MostlySparse2DHashmap<uint16_t, size_t> edgeOvlp;
	size_t edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw) const;
	size_t edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void setEdgeCoverage(size_t from, bool fromFw, size_t to, bool toFw, size_t val);
	void setEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t val);
	size_t edgeOverlap(size_t from, bool fromFw, size_t to, bool toFw) const;
	size_t edgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void setEdgeOverlap(size_t from, bool fromFw, size_t to, bool toFw, size_t val);
	void setEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t val);
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

}

#endif
