#ifndef MBG_h
#define MBG_h

#include <vector>
#include <string>
#include "ReadHelper.h"

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const ErrorMasking errorMasking, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths, const size_t maxResolveLength, const bool blunt, const size_t maxUnconditionalResolveLength, const std::string& nodeNamePrefix, const std::string& sequenceCacheFile, const bool keepGaps, const double hpcVariantOnecopyCoverage, const bool guesswork, const bool copycountFilterHeuristic);

#endif
