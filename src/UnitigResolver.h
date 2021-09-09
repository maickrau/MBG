#ifndef UnitigResolver_h
#define UnitigResolver_h

#include <tuple>
#include <vector>
#include "UnitigGraph.h"
#include "HashList.h"

class HashPath
{
public:
	std::vector<HashType> hashes;
	std::string readName;
	std::vector<size_t> hashPoses;
};

class ReadPath
{
public:
	std::string readName;
	std::vector<std::pair<size_t, bool>> path;
	std::vector<size_t> readPoses;
	size_t leftClip;
	size_t rightClip;
private:
};

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, const std::vector<HashPath>& paths, const size_t minCoverage);

#endif
