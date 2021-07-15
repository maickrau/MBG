#ifndef BluntGraph_h
#define BluntGraph_h

#include <vector>
#include <string>
#include <tuple>
#include "HashList.h"
#include "UnitigGraph.h"

// DBG bluntifying algorithm from Hassan Nikaein (personal communication)
// generalized to arbitrary overlaps
class BluntGraph
{
public:
	BluntGraph(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex);
	std::vector<std::string> nodes;
	std::vector<std::tuple<size_t, bool, size_t, bool, size_t>> edges;
	std::vector<float> nodeAvgCoverage;
private:
	void initializeNodes(const std::vector<CompressedSequenceType>& unbluntSequences, const StringIndex& stringIndex, const VectorWithDirection<size_t>& maxOverlap, const UnitigGraph& unitigs);
	void initializeEdgesAndEdgenodes(const HashList& hashlist, const UnitigGraph& unitigs, const VectorWithDirection<size_t>& maxOverlap, const std::vector<CompressedSequenceType>& unbluntSequences, const StringIndex& stringIndex);
	VectorWithDirection<size_t> getMaxOverlaps(const HashList& hashlist, const UnitigGraph& unitigs) const;
};

#endif
