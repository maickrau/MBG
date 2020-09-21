#ifndef TransitiveCleaner_h
#define TransitiveCleaner_h

#include <vector>
#include <unordered_set>
#include "VectorWithDirection.h"
#include "HashList.h"
#include "LazyString.h"

class TransitiveCleaner
{
public:
	TransitiveCleaner(size_t kmerSize, const HashList& hashlist);
	std::vector<std::pair<size_t, bool>> insertMiddles(std::vector<std::pair<size_t, bool>> raw) const;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> newSequenceOverlaps;
private:
	void addMiddles(size_t kmerSize, std::pair<size_t, bool> start, std::pair<size_t, bool> end, LazyString& seq, const HashList& list, const std::unordered_set<size_t>& minimizerPrefixes);
	std::unordered_set<size_t> getMinimizerPrefixes(size_t kmerSize, const HashList& hashlist);
	void getMiddles(size_t kmerSize, const HashList& hashlist);
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, std::vector<std::pair<size_t, bool>>>> transitiveMiddle;
};

#endif
