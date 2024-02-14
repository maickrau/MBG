#ifndef UnitigHelper_h
#define UnitigHelper_h

#include <tuple>
#include "HashList.h"
#include "UnitigGraph.h"

namespace MBG
{

size_t getUnitigOverlap(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to);
size_t getUnitigSize(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const size_t unitig);

}

#endif
