#ifndef ContextMerStorage_h
#define ContextMerStorage_h

#include <iostream>
#include <vector>
#include <phmap.h>
#include "MBGCommon.h"
#include "CommonUtils.h"

class ReadpartIterator;

class ContextMerStorage
{
public:
	ContextMerStorage();
	ContextMerStorage(size_t contextk, size_t contextw, size_t rawk, size_t raww, size_t numWindows, const ErrorMasking masking, const size_t minCoverage);
	void setParams(size_t contextk, size_t contextw, size_t rawk, size_t raww, size_t numWindows, const ErrorMasking masking, const size_t minCoverage);
	void constructStorage(const std::vector<std::string>& readFiles, const ReadpartIterator* const partIterator, const size_t numThreads);
	template <typename F>
	void iterateHashesInContext(const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, std::vector<HashType>& hashes, F callback) const
	{
		auto fwContextMers = getCoveredContextMers(rawSeq);
		auto bwContextMers = getCoveredContextMers(CommonUtils::ReverseComplement(rawSeq));
		std::reverse(bwContextMers.begin(), bwContextMers.end());
		for (size_t i = 0; i < bwContextMers.size(); i++)
		{
			std::swap(std::get<1>(bwContextMers[i]), std::get<2>(bwContextMers[i]));
			std::get<1>(bwContextMers[i]) = rawSeq.size()-1-std::get<1>(bwContextMers[i]);
			std::get<2>(bwContextMers[i]) = rawSeq.size()-1-std::get<2>(bwContextMers[i]);
		}
		for (size_t i = 1; i < fwContextMers.size(); i++)
		{
			assert(std::get<1>(fwContextMers[i]) >= std::get<1>(fwContextMers[i-1]));
			assert(std::get<2>(fwContextMers[i]) > std::get<1>(fwContextMers[i]));
		}
		for (size_t i = 1; i < bwContextMers.size(); i++)
		{
			assert(std::get<2>(bwContextMers[i]) > std::get<1>(bwContextMers[i]));
			assert(std::get<1>(bwContextMers[i]) >= std::get<1>(bwContextMers[i-1]));
		}
		size_t fwContextmerPos = 0;
		size_t bwContextmerPos = 0;
		std::vector<bool> deletePositions;
		deletePositions.resize(positions.size(), false);
		for (size_t i = 0; i < positions.size(); i++)
		{
			size_t startPos = poses[positions[i]];
			size_t endPos = poses[positions[i]+rawk]-1; // inclusive!
			while (fwContextmerPos < fwContextMers.size() && std::get<2>(fwContextMers[fwContextmerPos]) < endPos) fwContextmerPos += 1;
			while (bwContextmerPos < bwContextMers.size() && std::get<2>(bwContextMers[bwContextmerPos]) < endPos) bwContextmerPos += 1;
//			uint64_t fwHashModifier = 0;
//			uint64_t bwHashModifier = 0;
			HashType fwHash = hashes[i];
			HashType bwHash = (fwHash << 64) + (fwHash >> 64);
			HashType key = std::min(fwHash, bwHash);
			bool hasAnyHash = false;
			uint64_t anyHash = 0;
			if (fwContextmerPos < fwContextMers.size() && std::get<1>(fwContextMers[fwContextmerPos]) <= startPos)
			{
				anyHash = std::get<0>(fwContextMers[fwContextmerPos]);
				hasAnyHash = true;
//				fwHashModifier = findNoMutate(fwHash, std::get<0>(fwContextMers[fwContextmerPos]));
			}
			if (bwContextmerPos < bwContextMers.size() && std::get<1>(bwContextMers[bwContextmerPos]) <= startPos)
			{
				if (hasAnyHash) assert(findNoMutate(key, anyHash) == findNoMutate(key, std::get<0>(bwContextMers[bwContextmerPos])));
				anyHash = std::get<0>(bwContextMers[bwContextmerPos]);
				hasAnyHash = true;
//				bwHashModifier = findNoMutate(bwHash, std::get<0>(bwContextMers[bwContextmerPos]));
			}
			if (!hasAnyHash)
			{
				deletePositions[i] = true;
				continue;
			}
			uint64_t hashKey = findNoMutate(key, anyHash);
			for (size_t j = fwContextmerPos; j < fwContextMers.size(); j++)
			{
				if (std::get<1>(fwContextMers[j]) > startPos) break;
				assert(std::get<2>(fwContextMers[j]) >= endPos);
				assert(findNoMutate(key, std::get<0>(fwContextMers[j])) == hashKey);
			}
			for (size_t j = bwContextmerPos; j < bwContextMers.size(); j++)
			{
				if (std::get<1>(bwContextMers[j]) > startPos) break;
				assert(std::get<2>(bwContextMers[j]) >= endPos);
				assert(findNoMutate(key, std::get<0>(bwContextMers[j])) == hashKey);
			}
//			hashes[i] ^= (HashType)fwHashModifier << (HashType)64;
//			hashes[i] ^= (HashType)bwHashModifier;
			hashes[i] ^= (HashType)hashKey << (HashType)64;
			hashes[i] ^= (HashType)hashKey;
		}
		size_t lastNoDelete = 0;
		for (size_t i = 0; i < positions.size(); i++)
		{
			if (deletePositions[i])
			{
				if (i > lastNoDelete)
				{
					std::vector<size_t> newPositions { positions.begin()+lastNoDelete, positions.begin()+i };
					std::vector<HashType> newHashes { hashes.begin() + lastNoDelete, hashes.begin() + i };
					callback(read, seq, poses, rawSeq, newPositions, newHashes);
				}
				lastNoDelete = i+1;
			}
		}
		if (lastNoDelete < positions.size())
		{
			std::vector<size_t> newPositions { positions.begin()+lastNoDelete, positions.end() };
			std::vector<HashType> newHashes { hashes.begin() + lastNoDelete, hashes.end() };
			callback(read, seq, poses, rawSeq, newPositions, newHashes);
		}
//		callback(read, seq, poses, rawSeq, positions, hashes);
	}
private:
	std::vector<std::tuple<uint64_t, size_t, size_t>> getCoveredContextMers(const std::string& rawSeq) const;
	std::vector<std::tuple<uint64_t, size_t, size_t>> getRawContextMers(const std::string& rawSeq) const;
	uint64_t findNoMutate(HashType kmerHash, uint64_t contextHash) const;
	uint64_t find(HashType kmerHash, uint64_t contextHash);
	void merge(HashType kmerHash, uint64_t leftContextHash, uint64_t rightContextHash);
	phmap::flat_hash_map<HashType, phmap::flat_hash_map<uint64_t, uint64_t>> contextParentPerKmer;
	phmap::flat_hash_map<uint64_t, uint8_t> contextMerCoverage;
	ErrorMasking masking;
	size_t minCoverage;
	size_t contextk;
	size_t contextw;
	size_t numWindows;
	size_t rawk;
	size_t raww;
};

#endif
