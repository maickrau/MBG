#ifndef ConsensusMaker_h
#define ConsensusMaker_h

#include <vector>
#include <string>
#include <mutex>
#include <cassert>
#include <phmap.h>
#include "MBGCommon.h"

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
	std::vector<CompressedSequenceType> getSequences();
	template <typename F>
	void addStrings(size_t unitig, size_t unitigStart, size_t unitigEnd, F sequenceGetter)
	{
		assert(unitig < expandedCounts.size());
		assert(unitigEnd > unitigStart);
		assert(unitigEnd <= expandedCounts[unitig].size());
		size_t lowMutexIndex = unitigStart / MutexLength;
		if (unitigStart > 64) lowMutexIndex = (unitigStart - 64) / MutexLength;
		size_t highMutexIndex = (unitigEnd + 64 + MutexLength - 1) / MutexLength;
		if (highMutexIndex >= seqMutexes[unitig].size()) highMutexIndex = seqMutexes[unitig].size();
		std::vector<std::lock_guard<std::mutex>*> guards;
		for (size_t i = lowMutexIndex; i < highMutexIndex; i++)
		{
			guards.emplace_back(new std::lock_guard<std::mutex>{*seqMutexes[unitig][i]});
		}
		for (size_t i = 0; i < unitigEnd - unitigStart; i++)
		{
			size_t off = unitigStart + i;
			uint16_t compressed;
			std::string expanded;
			std::tie(compressed, expanded) = sequenceGetter(i);
			assert(compressedSequences[unitig][off] == 0 || compressedSequences[unitig][off] == compressed);
			compressedSequences[unitig][off] = compressed;
			expandedCounts[unitig][off][expanded] += 1;
		}
		for (size_t i = 0; i < guards.size(); i++)
		{
			delete guards[i];
		}
	}
private:
	std::vector<std::vector<phmap::flat_hash_map<std::string, size_t>>> expandedCounts;
	std::vector<std::vector<std::mutex*>> seqMutexes;
	std::vector<std::vector<uint16_t>> compressedSequences;
};

#endif
