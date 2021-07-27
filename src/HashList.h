#ifndef AdjacentMinimizerList_h
#define AdjacentMinimizerList_h

#include <memory>
#include <vector>
#include <string>
#include <mutex>
#include <phmap.h>
#include "MBGCommon.h"
#include "VectorWithDirection.h"
#include "LittleBigVector.h"
#include "MostlySparse2DHashmap.h"
#include "RankBitvector.h"

class HashList
{
public:
	HashList(size_t kmerSize);
	LittleBigVector<uint8_t, size_t> coverage;
	phmap::flat_hash_map<HashType, std::pair<size_t, bool>> hashToNode;
	size_t numSequenceOverlaps() const;
	size_t getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void setEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t coverage);
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> getEdgeCoverages(std::pair<size_t, bool> from) const;
	size_t getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	std::vector<std::pair<std::pair<size_t, bool>, size_t>> getSequenceOverlaps(std::pair<size_t, bool> from) const;
	bool hasSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap);
	void addEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	size_t size() const;
	void resize(size_t size);
	std::pair<size_t, bool> getNodeOrNull(VectorView<CharType> sequence) const;
	std::pair<std::pair<size_t, bool>, HashType> addNode(VectorView<CharType> sequence, VectorView<CharType> reverse, HashType previousHash, size_t overlap, uint64_t bucketHash);
	void filter(const RankBitvector& kept);
private:
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	MostlySparse2DHashmap<uint16_t, size_t> sequenceOverlap;
	std::shared_ptr<std::mutex> indexMutex;
	size_t kmerSize;
};

#endif
