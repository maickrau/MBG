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
	// allow multiple threads to update the same contig sequence but in different regions
	// each mutex covers MutexLength bp in one contig
	// ------
	//       ------
	//             ------
	// etc
	// size arbitrarily 1Mbp, so ~(3000 + num_contigs) mutexes in a human genome, hopefully not too many
	// and the chance of two random hifis falling in the same 1Mbp bucket is ~.03% so hopefully not too much waiting
	static constexpr size_t MutexLength = 1000000;
public:
	~ConsensusMaker();
	void init(const std::vector<size_t>& unitigLens);
	std::pair<std::vector<CompressedSequenceType>, StringIndex> getSequences();
	template <typename F>
	void addStrings(size_t unitig, size_t unitigStart, size_t unitigEnd, F sequenceGetter)
	{
		assert(unitig < simpleCounts.size());
		assert(unitigEnd > unitigStart);
		assert(unitigEnd <= simpleCounts[unitig].size());
		std::vector<std::pair<uint16_t, uint32_t>> processedChars;
		processedChars.resize(unitigEnd - unitigStart);
		{
			std::lock_guard<std::mutex> lock { stringIndexMutex };
			for (size_t i = 0; i < unitigEnd - unitigStart; i++)
			{
				uint16_t compressed;
				std::string expanded;
				std::tie(compressed, expanded) = sequenceGetter(i);
				processedChars[i].first = compressed;
				processedChars[i].second = stringIndex.getIndex(compressed, expanded);
			}
		}
		// size_t lowMutexIndex = unitigStart / MutexLength;
		// if (unitigStart > 64) lowMutexIndex = (unitigStart - 64) / MutexLength;
		// size_t highMutexIndex = (unitigEnd + 64 + MutexLength - 1) / MutexLength;
		// if (highMutexIndex >= seqMutexes[unitig].size()) highMutexIndex = seqMutexes[unitig].size();
		// std::vector<std::lock_guard<std::mutex>*> guards;
		// assert(unitig < seqMutexes.size());
		// assert(lowMutexIndex < seqMutexes[unitig].size());
		// assert(highMutexIndex <= seqMutexes[unitig].size());
		// assert(highMutexIndex > lowMutexIndex);
		// for (size_t i = lowMutexIndex; i < highMutexIndex; i++)
		// {
		// 	guards.emplace_back(new std::lock_guard<std::mutex>{*seqMutexes[unitig][i]});
		// }
		std::vector<std::pair<uint32_t, uint32_t>> addComplexes;
		for (size_t i = 0; i < processedChars.size(); i++)
		{
			size_t off = unitigStart + i;
			auto found = find(unitig, off);
			size_t realUnitig = std::get<0>(found);
			size_t realOff = std::get<1>(found);
			uint16_t compressed = processedChars[i].first;
			uint32_t expandedIndex = processedChars[i].second;
			if (!std::get<2>(found))
			{
				expandedIndex = stringIndex.getReverseIndex(compressed, expandedIndex);
				compressed = complement(compressed);
			}
			size_t compressIndex = realOff / MutexLength;
			assert(compressedSequences[realUnitig][compressIndex].get(realOff % MutexLength) == 0 || compressedSequences[realUnitig][compressIndex].get(realOff % MutexLength) == compressed);
			compressedSequences[realUnitig][compressIndex].set(realOff % MutexLength, compressed);
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
		// for (size_t i = 0; i < guards.size(); i++)
		// {
		// 	delete guards[i];
		// }
	}
	void addEdgeOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, size_t overlap);
private:
	std::tuple<size_t, size_t, bool> find(size_t unitig, size_t index) const;
	void merge(size_t leftUnitig, size_t leftIndex, size_t rightUnitig, size_t rightIndex, bool fw);
	StringIndex stringIndex;
	std::vector<std::vector<std::pair<uint8_t, uint8_t>>> simpleCounts;
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint32_t, uint32_t>>> complexCounts;
	std::vector<std::mutex*> complexCountMutexes;
	std::vector<std::vector<std::mutex*>> seqMutexes;
	std::vector<std::vector<TwobitLittleBigVector<uint16_t>>> compressedSequences;
	std::mutex stringIndexMutex;
	mutable std::vector<std::vector<std::tuple<size_t, size_t, bool>>> parent;
};

#endif
