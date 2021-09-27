#ifndef UnitigResolver_h
#define UnitigResolver_h

#include <tuple>
#include <vector>
#include "UnitigGraph.h"
#include "HashList.h"
#include "CumulativeVector.h"

class HashPath
{
public:
	std::vector<HashType> hashes;
	std::string readName;
	size_t readLength;
	CumulativeVector<uint8_t> hashPoses;
	CumulativeVector<uint8_t> hashPosesExpandedStart;
	CumulativeVector<uint8_t> hashPosesExpandedEnd;
};

class ReadPath
{
public:
	std::string readName;
	std::vector<std::pair<size_t, bool>> path;
	std::vector<size_t> readPoses;
	std::vector<size_t> readPosesExpandedStart;
	std::vector<size_t> readPosesExpandedEnd;
	size_t leftClip;
	size_t rightClip;
	size_t readLength;
private:
};

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, const std::vector<HashPath>& paths, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength);

#endif
