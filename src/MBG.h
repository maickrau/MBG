#ifndef MBG_h
#define MBG_h

#include <vector>
#include <string>

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool blunt, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths);

#endif
