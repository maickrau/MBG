#include "ContextMerStorage.h"
#include "ContextMerIterator.h"

ContextMerStorage::ContextMerStorage()
{
}

ContextMerStorage::ContextMerStorage(size_t contextk, size_t contextw, size_t rawk, size_t raww, size_t numWindows, const ErrorMasking masking, const size_t minCoverage)
{
	setParams(contextk, contextw, rawk, raww, numWindows, masking, minCoverage);
}

void ContextMerStorage::setParams(size_t contextk, size_t contextw, size_t rawk, size_t raww, const size_t numWindows, const ErrorMasking masking, const size_t minCoverage)
{
	assert(contextParentPerKmer.size() == 0);
	assert(contextMerCoverage.size() == 0);
	this->contextk = contextk;
	this->contextw = contextw;
	this->rawk = rawk;
	this->raww = raww;
	this->numWindows = numWindows;
	this->masking = masking;
	this->minCoverage = minCoverage;
}

void ContextMerStorage::constructStorage(const std::vector<std::string>& readFiles, const ReadpartIterator* const partIterator, const size_t numThreads)
{
	assert(contextParentPerKmer.size() == 0);
	assert(contextMerCoverage.size() == 0);
	std::mutex indexMutex;
	iterateReadsMultithreaded(readFiles, numThreads, [this, &indexMutex](ReadInfo readInfo, const std::string& sequence)
	{
		auto fwHashes = getRawContextMers(sequence);
		auto bwHashes = getRawContextMers(CommonUtils::ReverseComplement(sequence));
		std::lock_guard<std::mutex> lock { indexMutex };
		for (auto hash : fwHashes)
		{
			if (contextMerCoverage[std::get<0>(hash)] < minCoverage) contextMerCoverage[std::get<0>(hash)] += 1;
		}
		for (auto hash : bwHashes)
		{
			if (contextMerCoverage[std::get<0>(hash)] < minCoverage) contextMerCoverage[std::get<0>(hash)] += 1;
		}
	});
	//ReadpartIterator partIterator { rawk, raww, masking, numThreads, readFiles, false, "" };
	partIterator->iterateHashes([this, &indexMutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		auto fwHashes = getCoveredContextMers(rawSeq);
		auto bwHashes = getCoveredContextMers(CommonUtils::ReverseComplement(rawSeq));
		std::reverse(bwHashes.begin(), bwHashes.end());
		for (size_t i = 0; i < bwHashes.size(); i++)
		{
			std::swap(std::get<1>(bwHashes[i]), std::get<2>(bwHashes[i]));
			std::get<1>(bwHashes[i]) = rawSeq.size()-1-std::get<1>(bwHashes[i]);
			std::get<2>(bwHashes[i]) = rawSeq.size()-1-std::get<2>(bwHashes[i]);
		}
		for (size_t i = 1; i < fwHashes.size(); i++)
		{
			assert(std::get<1>(fwHashes[i]) >= std::get<1>(fwHashes[i-1]));
			assert(std::get<2>(fwHashes[i]) > std::get<1>(fwHashes[i]));
		}
		for (size_t i = 1; i < bwHashes.size(); i++)
		{
			assert(std::get<2>(bwHashes[i]) > std::get<1>(bwHashes[i]));
			assert(std::get<1>(bwHashes[i]) >= std::get<1>(bwHashes[i-1]));
		}
		size_t fwStart = 0;
		size_t fwEnd = 0;
		size_t bwStart = 0;
		size_t bwEnd = 0;
		for (size_t i = 0; i < positions.size(); i++)
		{
			size_t startPos = poses[positions[i]];
			size_t endPos = poses[positions[i]+rawk]-1; // inclusive!
			while (fwStart < fwHashes.size() && std::get<2>(fwHashes[fwStart]) < endPos) fwStart += 1;
			fwEnd = std::max(fwEnd, fwStart);
			while (fwEnd < fwHashes.size() && std::get<1>(fwHashes[fwEnd]) <= startPos) fwEnd += 1;
			while (bwStart < bwHashes.size() && std::get<2>(bwHashes[bwStart]) < endPos) bwStart += 1;
			bwEnd = std::max(bwEnd, bwStart);
			while (bwEnd < bwHashes.size() && std::get<1>(bwHashes[bwEnd]) <= startPos) bwEnd += 1;
			HashType fwHash = hashes[i];
			HashType bwHash = (fwHash << 64) + (fwHash >> 64);
			HashType key = std::min(fwHash, bwHash);
			std::lock_guard<std::mutex> guard { indexMutex };
			uint64_t anyHash = 0;
			bool hasAnyHash = false;
			if (fwEnd > fwStart)
			{
				anyHash = std::get<0>(fwHashes[fwStart]);
				hasAnyHash = true;
			}
			else if (bwEnd > bwStart)
			{
				anyHash = std::get<0>(bwHashes[bwStart]);
				hasAnyHash = true;
			}
			if (hasAnyHash) contextParentPerKmer[key];
			for (size_t j = fwStart; j < fwEnd; j++)
			{
				merge(key, anyHash, std::get<0>(fwHashes[j]));
			}
			for (size_t j = bwStart; j < bwEnd; j++)
			{
				merge(key, anyHash, std::get<0>(bwHashes[j]));
			}
/*			if (fwEnd == fwStart+1) contextParentPerKmer[fwHash]; // make sure the k-mer level map exists if there's just one context hash
			if (bwEnd == bwStart+1) contextParentPerKmer[bwHash]; // make sure the k-mer level map exists if there's just one context hash
			for (size_t i = fwStart+1; i < fwEnd; i++)
			{
				merge(fwHash, std::get<0>(fwHashes[fwStart]), std::get<0>(fwHashes[i]));
			}
			for (size_t i = bwStart+1; i < bwEnd; i++)
			{
				merge(bwHash, std::get<0>(bwHashes[fwStart]), std::get<0>(bwHashes[i]));
			}*/
		}
	});
	for (auto& pair : contextParentPerKmer)
	{
		phmap::flat_hash_set<uint64_t> keys;
		for (auto& pair2 : pair.second)
		{
			keys.insert(pair2.first);
		}
		for (auto key : keys)
		{
			find(pair.first, key);
		}
	}
}

// end pos is inclusive
std::vector<std::tuple<uint64_t, size_t, size_t>> ContextMerStorage::getRawContextMers(const std::string& rawSeq) const
{
	std::vector<std::tuple<uint64_t, size_t, size_t>> result;
	iterateContextHashes(rawSeq, masking, contextk, contextw, numWindows, [&result](uint64_t hash, size_t start, size_t end)
	{
		result.emplace_back(hash, start, end);
	});
	return result;
}

// end pos is inclusive
std::vector<std::tuple<uint64_t, size_t, size_t>> ContextMerStorage::getCoveredContextMers(const std::string& rawSeq) const
{
	std::vector<std::tuple<uint64_t, size_t, size_t>> result;
	iterateContextHashes(rawSeq, masking, contextk, contextw, numWindows, [this, &result](uint64_t hash, size_t start, size_t end)
	{
		if (contextMerCoverage.at(hash) < minCoverage) return;
		result.emplace_back(hash, start, end);
	});
	return result;
}

uint64_t ContextMerStorage::findNoMutate(HashType kmerHash, uint64_t contextHash) const
{
	const phmap::flat_hash_map<uint64_t, uint64_t>& map = contextParentPerKmer.at(kmerHash);
	if (map.count(contextHash) == 0) return contextHash;
	return map.at(contextHash); // assume the parent links are correct at this point
}

uint64_t ContextMerStorage::find(HashType kmerHash, uint64_t contextHash)
{
	uint64_t original = contextHash;
	phmap::flat_hash_map<uint64_t, uint64_t>& map = contextParentPerKmer[kmerHash];
	if (map.count(contextHash) == 0) return contextHash;
	while (map.count(contextHash) == 1 && map.at(contextHash) != contextHash)
	{
		contextHash = map.at(contextHash);
	}
	map[original] = contextHash;
	return contextHash;
}

void ContextMerStorage::merge(HashType kmerHash, uint64_t leftContextHash, uint64_t rightContextHash)
{
	phmap::flat_hash_map<uint64_t, uint64_t>& map = contextParentPerKmer[kmerHash];
	auto left = find(kmerHash, leftContextHash);
	auto right = find(kmerHash, rightContextHash);
	map[right] = left;
}
