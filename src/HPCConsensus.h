#ifndef HPCConsensus_h
#define HPCConsensus_h

#include <vector>
#include <tuple>
#include <string>
#include "HashList.h"
#include "UnitigGraph.h"
#include "ReadHelper.h"
#include "MBGCommon.h"
#include "StringIndex.h"
#include "UnitigResolver.h"

std::pair<std::vector<CompressedSequenceType>, StringIndex> getHPCUnitigSequences(const HashList& hash, const UnitigGraph& unitigs, std::vector<ReadPath>& readPaths, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads);
void getHpcVariantsAndReadPaths(const HashList& hash, const UnitigGraph& unitigs, const size_t kmerSize, ReadpartIterator& partIterator, const size_t numThreads, const double minUnitigCoverage, const size_t minVariantCoverage);

#endif
