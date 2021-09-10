#ifndef ConsensusMaker_h
#define ConsensusMaker_h

#include <vector>
#include <string>
#include <mutex>
#include <cassert>
#include <phmap.h>
#include "MBGCommon.h"
#include "StringIndex.h"
#include "ErrorMaskHelper.h"

class ConsensusMaker
{
private:
public:
	void init(const std::vector<size_t>& unitigLens);
	std::pair<std::vector<CompressedSequenceType>, StringIndex> getSequences();
	template <typename F>
	void addStrings(size_t unitig, size_t unitigStart, size_t unitigEnd, F sequenceGetter)
	{
		assert(unitig < simpleCounts.size());
		assert(unitigEnd > unitigStart);
		assert(unitigEnd <= simpleCounts[unitig].size());
		for (size_t i = 0; i < unitigEnd - unitigStart; i++)
		{
			uint16_t compressed;
			std::string expanded;
			std::tie(compressed, expanded) = sequenceGetter(i);
			size_t expandedIndex = stringIndex.getIndex(compressed, expanded);
			size_t off = unitigStart + i;
			auto found = find(unitig, off);
			size_t realUnitig = std::get<0>(found);
			size_t realOff = std::get<1>(found);
			if (!std::get<2>(found))
			{
				expandedIndex = stringIndex.getReverseIndex(compressed, expandedIndex);
				compressed = complement(compressed);
			}
			assert(compressedSequences[realUnitig].get(realOff) == 0 || compressedSequences[realUnitig].get(realOff) == compressed);
			compressedSequences[realUnitig].set(realOff, compressed);
			bool didSimple = false;
			if (expandedIndex < 256)
			{
				if (simpleCounts[realUnitig][realOff].second == 0 || (expandedIndex == simpleCounts[realUnitig][realOff].first && simpleCounts[realUnitig][realOff].second < 255))
				{
					simpleCounts[realUnitig][realOff].first = expandedIndex;
					simpleCounts[realUnitig][realOff].second += 1;
					didSimple = true;
				}
			}
			if (!didSimple)
			{
				complexCounts[realUnitig][realOff][expandedIndex] += 1;
			}
		}
	}
	void addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap);
private:
	std::tuple<size_t, size_t, bool> find(size_t unitig, size_t index) const;
	void merge(size_t leftUnitig, size_t leftIndex, size_t rightUnitig, size_t rightIndex, bool fw);
	StringIndex stringIndex;
	std::vector<std::vector<std::pair<uint8_t, uint8_t>>> simpleCounts;
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint32_t, uint32_t>>> complexCounts;
	std::vector<TwobitLittleBigVector<uint16_t>> compressedSequences;
	mutable std::vector<std::vector<std::tuple<size_t, size_t, bool>>> parent;
};

#endif
