#ifndef ContextMerIterator_h
#define ContextMerIterator_h

#include <vector>
#include <tuple>
#include <cstdint>
#include "MBGCommon.h"
#include "FastHasher.h"
#include "ReadHelper.h"

class MinimizerIterator
{
public:
	MinimizerIterator() = default;
	MinimizerIterator(size_t k);
	void init(const SequenceCharType& start, size_t posOffset);
	void moveChar(uint16_t added, uint16_t removed);
	size_t minimizerPosition() const;
	uint64_t minimizerHash() const;
private:
	size_t k;
	size_t windowSize;
	std::vector<std::pair<size_t, uint64_t>> windowKmers;
	FastHasher lastHash;
	size_t pos;
};

template <typename F>
void iterateContextWindowMerlists(const SequenceCharType& seq, size_t k, size_t numWindows, size_t windowSize, F callback)
{
	if (seq.size() < numWindows * windowSize + k) return;
	std::vector<MinimizerIterator> windowIterators;
	for (size_t i = 0; i < numWindows; i++)
	{
		windowIterators.emplace_back(k);
		SequenceCharType startWindow { seq.begin() + i * windowSize, seq.begin() + (i+1) * windowSize + k };
		windowIterators[i].init(startWindow, i * windowSize);
	}
	std::vector<uint64_t> hashesHere;
	hashesHere.resize(numWindows);
	for (size_t i = 0; i < numWindows; i++)
	{
		hashesHere[i] = windowIterators[i].minimizerHash();
	}
	callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
	for (size_t i = 0; i + numWindows * windowSize + k < seq.size(); i++)
	{
		bool changed = false;
		for (size_t j = 0; j < numWindows; j++)
		{
			windowIterators[j].moveChar(seq[(j+1)*windowSize + k + i], seq[(j+1)*windowSize+i]);
			if (windowIterators[j].minimizerHash() != hashesHere[j])
			{
				changed = true;
				hashesHere[j] = windowIterators[j].minimizerHash();
			}
		}
		if (changed)
		{
			callback(hashesHere, windowIterators[0].minimizerPosition(), windowIterators.back().minimizerPosition()+k-1);
		}
	}
}

template <typename F>
void iterateContextHashes(const std::string& readSequence, const ErrorMasking errorMasking, const size_t k, const size_t w, const size_t numWindows, F callback)
{
	std::vector<std::string> readFiles { };
	ReadpartIterator partIterator { 31, 1, errorMasking, 1, readFiles, false, "" };
	partIterator.iteratePartsOfRead("", readSequence, [callback, k, numWindows, w](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
	{
		if (seq.size() < numWindows * w + k) return;
		iterateContextWindowMerlists(seq, k, numWindows, w, [callback, &poses](const std::vector<uint64_t>& hashes, const size_t startPos, const size_t endPos)
		{
			uint64_t totalhash = 0;
			for (auto hash : hashes)
			{
				totalhash *= 3;
				totalhash += hash;
			}
			assert(endPos > startPos);
			assert(startPos < poses.size());
			assert(endPos < poses.size());
			assert(endPos+1 < poses.size());
			size_t realStartPos = poses[startPos];
			size_t realEndPos = poses[endPos+1]-1; // inclusive!
			callback(totalhash, realStartPos, realEndPos);
		});
	});
}
#endif