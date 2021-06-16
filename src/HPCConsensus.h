#ifndef HPCConsensus_h
#define HPCConsensus_h

#include <vector>
#include <tuple>
#include <string>
#include "HashList.h"
#include "UnitigGraph.h"
#include "ReadHelper.h"

std::vector<std::pair<std::vector<uint16_t>, std::vector<uint8_t>>> getHPCUnitigSequences(const HashList& hash, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads);

#endif
