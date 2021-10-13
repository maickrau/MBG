#include "UnitigHelper.h"

size_t getUnitigOverlap(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	auto kmerOverlap = unitigs.edgeOverlap(from, to);
	if (kmerOverlap == 0)
	{
		std::pair<size_t, bool> fromKmer = unitigs.unitigs[from.first].back();
		if (!from.second) fromKmer = reverse(unitigs.unitigs[from.first][0]);
		std::pair<size_t, bool> toKmer = unitigs.unitigs[to.first][0];
		if (!to.second) toKmer = reverse(unitigs.unitigs[to.first].back());
		return hashlist.getOverlap(fromKmer, toKmer);
	}
	assert(unitigs.unitigs[from.first].size() >= kmerOverlap);
	size_t result = kmerSize;
	for (size_t i = 1; i < kmerOverlap; i++)
	{
		std::pair<size_t, bool> prev = unitigs.unitigs[from.first][unitigs.unitigs[from.first].size() - kmerOverlap + i - 1];
		std::pair<size_t, bool> curr = unitigs.unitigs[from.first][unitigs.unitigs[from.first].size() - kmerOverlap + i];
		if (!from.second)
		{
			prev = unitigs.unitigs[from.first][i-1];
			curr = unitigs.unitigs[from.first][i];
		}
		result += kmerSize - hashlist.getOverlap(prev, curr);
	}
	return result;
}
