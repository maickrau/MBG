#ifndef HPCConsensus_h
#define HPCConsensus_h

#include <vector>
#include <tuple>
#include <string>
#include "HashList.h"
#include "UnitigGraph.h"
#include "ReadHelper.h"

std::vector<std::pair<std::string, std::vector<uint8_t>>> getHPCUnitigSequences(const HashList& hash, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads);
std::string getExpandedSequence(const std::string& seq, const std::vector<uint8_t>& length);
std::string revCompRLE(const std::string& str);
size_t getOverlapFromRLE(const std::vector<std::pair<std::string, std::vector<uint8_t>>>& unitigSequences, std::pair<size_t, bool> fromUnitig, size_t rleOverlap);

#endif
