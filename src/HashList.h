#ifndef AdjacentMinimizerList_h
#define AdjacentMinimizerList_h

#include <vector>
#include <string>
#include <phmap.h>
#include "MBGCommon.h"
#include "TwobitString.h"
#include "VectorWithDirection.h"

class AdjacentMinimizerList
{
public:
	AdjacentMinimizerList();
	TwobitView getView(size_t coord1, size_t coord2, size_t size) const;
	std::pair<size_t, size_t> addString(std::string_view str, HashType currentHash, HashType previousHash, size_t overlap);
	AdjacentMinimizerList getReverseComplementStorage() const;
	std::pair<size_t, size_t> getRevCompLocation(size_t coord1, size_t coord2, size_t size) const;
private:
	std::vector<TwobitString> data;
	HashType lastHash;
};

class AdjacentLengthList
{
public:
	AdjacentLengthList();
	std::vector<uint16_t> getData(size_t coord1, size_t coord2, size_t size) const;
	std::pair<size_t, size_t> addData(const std::vector<uint16_t>& lens, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap);
	void addCounts(const std::vector<uint16_t>& lens, bool fw, size_t start, size_t end, size_t coord1, size_t coord2);
	size_t size() const;
private:
	std::vector<std::vector<uint16_t>> sums;
	std::vector<std::vector<uint8_t>> counts;
	HashType lastHash;
};

class HashList
{
public:
	HashList(size_t kmerSize, bool collapseRunLengths);
	std::vector<size_t> coverage;
	std::vector<uint64_t> fakeFwHashes;
	std::vector<uint64_t> fakeBwHashes;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> sequenceOverlap;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCoverage;
	phmap::flat_hash_map<HashType, std::pair<size_t, bool>> hashToNode;
	size_t numSequenceOverlaps() const;
	size_t getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	void addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap);
	size_t size() const;
	std::vector<uint16_t> getHashCharacterLength(size_t index) const;
	void addHashCharacterLength(const std::vector<uint16_t>& data, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap);
	void addHashCharacterLength(const std::vector<uint16_t>& data, bool fw, size_t start, size_t end, size_t node);
	TwobitView getHashSequenceRLE(size_t index) const;
	TwobitView getRevCompHashSequenceRLE(size_t index) const;
	void addHashSequenceRLE(std::string_view seq, HashType currentHash, HashType previousHash, size_t overlap);
	void buildReverseCompHashSequences();
	std::pair<size_t, bool> getNodeOrNull(std::string_view sequence) const;
	std::pair<std::pair<size_t, bool>, HashType> addNode( std::string_view sequence, std::string_view reverse, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap, uint64_t fakeFwHash, uint64_t fakeBwHash);
private:
	AdjacentLengthList hashCharacterLengths;
	std::vector<std::pair<size_t, size_t>> hashCharacterLengthPtr;
	std::vector<size_t> hashCharacterCounts;
	AdjacentMinimizerList hashSequences;
	std::vector<std::pair<size_t, size_t>> hashSeqPtr;
	AdjacentMinimizerList hashSequencesRevComp;
	const size_t kmerSize;
	const bool collapseRunLengths;
};

#endif
