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
	~ConsensusMaker();
	void init(const std::vector<size_t>& unitigLens);
	std::pair<std::vector<CompressedSequenceType>, StringIndex> getSequences();
	void findParentLinks();
	template <typename F>
	void addStrings(size_t unitig, size_t unitigStart, size_t unitigEnd, F sequenceGetter)
	{
		assert(unitig < simpleCounts.size());
		assert(unitigEnd > unitigStart);
		assert(unitigEnd <= simpleCounts[unitig].size());
		std::vector<std::pair<uint16_t, uint32_t>> sequences;
		sequences.resize(unitigEnd - unitigStart);
		{
			std::lock_guard guard { stringIndexMutex };
			for (size_t i = 0; i < unitigEnd - unitigStart; i++)
			{
				uint16_t compressed;
				std::variant<size_t, std::string> expanded;
				std::tie(compressed, expanded) = sequenceGetter(i);
				uint32_t expandedIndex = stringIndex.getIndex(compressed, expanded);
				sequences[i].first = compressed;
				sequences[i].second = expandedIndex;
			}
		}
		std::vector<std::tuple<size_t, size_t, size_t>> complexes;
		size_t currentUnitig = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < unitigEnd - unitigStart; i++)
		{
			uint16_t compressed = sequences[i].first;
			size_t expandedIndex = sequences[i].second;
			size_t off = unitigStart + i;
			// no find because it might mutate parent, instead rely on parent being correct already
			auto found = parent[unitig][off];
			assert(std::get<0>(parent[std::get<0>(found)][std::get<1>(found)]) == std::get<0>(found));
			assert(std::get<1>(parent[std::get<0>(found)][std::get<1>(found)]) == std::get<1>(found));
			size_t realUnitig = std::get<0>(found);
			size_t realOff = std::get<1>(found);
			if (!std::get<2>(found))
			{
				expandedIndex = stringIndex.getReverseIndex(compressed, expandedIndex); // thread safe so no lock
				compressed = complement(compressed);
			}
			if (realUnitig != currentUnitig)
			{
				if (currentUnitig != std::numeric_limits<size_t>::max()) simpleSequenceMutexes[currentUnitig]->unlock();
				simpleSequenceMutexes[realUnitig]->lock();
				currentUnitig = realUnitig;
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
				complexes.emplace_back(realUnitig, realOff, expandedIndex);
			}
		}
		if (currentUnitig != std::numeric_limits<size_t>::max()) simpleSequenceMutexes[currentUnitig]->unlock();
		if (complexes.size() > 0)
		{
			std::lock_guard<std::mutex> guard { complexCountMutex };
			for (auto tuple : complexes)
			{
				complexCounts[std::get<0>(tuple)][std::get<1>(tuple)][std::get<2>(tuple)] += 1;
			}
		}
	}
	void addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap);
private:
	std::tuple<size_t, size_t, bool> find(size_t unitig, size_t index);
	void merge(size_t leftUnitig, size_t leftIndex, size_t rightUnitig, size_t rightIndex, bool fw);
	StringIndex stringIndex;
	std::vector<std::vector<std::pair<uint8_t, uint8_t>>> simpleCounts;
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint32_t, uint32_t>>> complexCounts;
	std::vector<std::mutex*> simpleSequenceMutexes;
	std::mutex complexCountMutex;
	std::mutex stringIndexMutex;
	std::vector<TwobitLittleBigVector<uint16_t>> compressedSequences;
	std::vector<std::pair<size_t, size_t>> needsComplementVerification;
	mutable std::vector<std::vector<std::tuple<size_t, size_t, bool>>> parent;
};

#endif
