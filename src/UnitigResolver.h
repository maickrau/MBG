#ifndef UnitigResolver_h
#define UnitigResolver_h

#include <tuple>
#include <string>
#include <vector>
#include "UnitigGraph.h"
#include "HashList.h"
#include "CumulativeVector.h"
#include "ReadHelper.h"

class ReadPath
{
public:
	std::string readName;
	std::vector<std::pair<size_t, bool>> path;
	std::vector<size_t> readPoses;
	size_t expandedReadPosStart;
	size_t expandedReadPosEnd;
	size_t leftClip;
	size_t rightClip;
	size_t readLength;
	size_t readLengthHPC;
private:
};

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, std::vector<ReadPath>& readPaths, const ReadpartIterator& partIterator, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool keepGaps, const bool guesswork);

#endif
