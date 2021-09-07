#ifndef UnitigResolver_h
#define UnitigResolver_h

#include <tuple>
#include <vector>
#include "UnitigGraph.h"
#include "HashList.h"

UnitigGraph resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, const std::vector<std::vector<HashType>>& paths, const size_t minCoverage);

#endif
