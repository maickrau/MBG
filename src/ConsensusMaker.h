#ifndef ConsensusMaker_h
#define ConsensusMaker_h

#include <vector>
#include <string>
#include <mutex>
#include <cassert>
#include <phmap.h>
#include "MBGCommon.h"
#include "StringIndex.h"

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
		size_t lowMutexIndex = unitigStart / MutexLength;
		if (unitigStart > 64) lowMutexIndex = (unitigStart - 64) / MutexLength;
		size_t highMutexIndex = (unitigEnd + 64 + MutexLength - 1) / MutexLength;
		if (highMutexIndex >= seqMutexes[unitig].size()) highMutexIndex = seqMutexes[unitig].size();
		std::vector<std::lock_guard<std::mutex>*> guards;
		for (size_t i = lowMutexIndex; i < highMutexIndex; i++)
		{
			guards.emplace_back(new std::lock_guard<std::mutex>{*seqMutexes[unitig][i]});
		}
		std::vector<std::pair<uint32_t, uint32_t>> addComplexes;
		for (size_t i = 0; i < processedChars.size(); i++)
		{
			size_t off = unitigStart + i;
			uint16_t compressed = processedChars[i].first;
			uint32_t expandedIndex = processedChars[i].second;
			assert(compressedSequences[unitig][off] == 0 || compressedSequences[unitig][off] == compressed);
			compressedSequences[unitig][off] = compressed;
			bool didSimple = false;
			if (expandedIndex < 256)
			{
				if (simpleCounts[unitig][off].second == 0 || (expandedIndex == simpleCounts[unitig][off].first && simpleCounts[unitig][off].second < 256))
				{
					simpleCounts[unitig][off].first = expandedIndex;
					simpleCounts[unitig][off].second += 1;
					didSimple = true;
				}
			}
			if (!didSimple)
			{
				addComplexes.emplace_back(off, expandedIndex);
			}
		}
		for (size_t i = 0; i < guards.size(); i++)
		{
			delete guards[i];
		}
		if (addComplexes.size() > 0)
		{
			std::lock_guard<std::mutex> lock { *complexCountMutexes[unitig] };
			for (auto pair : addComplexes)
			{
				complexCounts[unitig][pair] += 1;
			}
		}
	}
private:
	StringIndex stringIndex;
	std::vector<std::vector<std::pair<uint8_t, uint8_t>>> simpleCounts;
	std::vector<phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint32_t>> complexCounts;
	std::vector<std::mutex*> complexCountMutexes;
	std::vector<std::vector<std::mutex*>> seqMutexes;
	std::vector<std::vector<uint16_t>> compressedSequences;
	std::mutex stringIndexMutex;
};

#endif
