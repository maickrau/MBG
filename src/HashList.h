#ifndef AdjacentMinimizerList_h
#define AdjacentMinimizerList_h

#include <vector>
#include <string>
#include <mutex>
#include <phmap.h>
#include "MBGCommon.h"
#include "TwobitString.h"
#include "VectorWithDirection.h"

class HashList
{
public:
	HashList(size_t kmerSize);
	std::vector<size_t> coverage;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> sequenceOverlap;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCoverage;
	phmap::flat_hash_map<HashType, std::pair<size_t, bool>> hashToNode;
	size_t numSequenceOverlaps() const;
	size_t getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap);
	void addEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	size_t size() const;
	std::pair<size_t, bool> getNodeOrNull(std::string_view sequence) const;
	std::pair<std::pair<size_t, bool>, HashType> addNode(std::string_view sequence, std::string_view reverse, HashType previousHash, size_t overlap, uint64_t bucketHash);
private:
	std::mutex indexMutex;
	const size_t kmerSize;
};

#endif
