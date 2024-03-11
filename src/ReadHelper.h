#ifndef ReadHelper_h
#define ReadHelper_h

#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <atomic>
#include <thread>
#include <iostream>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include <unordered_set>
#include "fastqloader.h"
#include "FastHasher.h"
#include "ErrorMaskHelper.h"
#include "Serializer.h"

namespace MBG
{

extern thread_local std::vector<size_t> memoryIterables;

enum ErrorMasking
{
	No,
	Hpc,
	Collapse,
	Dinuc,
	Microsatellite,
	CollapseDinuc,
	CollapseMicrosatellite,
};

template <typename F, typename EdgeCheckFunction>
void findSyncmerPositions(const SequenceCharType& sequence, size_t kmerSize, size_t smerSize, std::vector<std::tuple<size_t, uint64_t>>& smerOrder, EdgeCheckFunction endSmer, F callback)
{
	if (sequence.size() < kmerSize) return;
	assert(smerSize <= kmerSize);
	size_t windowSize = kmerSize - smerSize + 1;
	assert(windowSize >= 1);
	FastHasher fwkmerHasher { smerSize };
	for (size_t i = 0; i < smerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	auto thisHash = fwkmerHasher.hash();
	if (endSmer(thisHash)) thisHash = 0;
	smerOrder.emplace_back(0, thisHash);
	for (size_t i = 1; i < windowSize; i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		uint64_t hash = fwkmerHasher.hash();
		// palindromic k-mer prevention: those s-mers which, if picked, might lead to picking a palindromic k-mer will arbitrarily have max hash value so they won't get picked
		if (smerSize*2 > kmerSize)
		{
			if (sequence[seqPos-kmerSize/2] == complement(sequence[seqPos-kmerSize/2]) || sequence[seqPos-smerSize+1+kmerSize/2] == complement(sequence[seqPos-smerSize+1+kmerSize/2]))
			{
				hash = std::numeric_limits<uint64_t>::max();
			}
		}
		if (endSmer(hash)) hash = 0;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == windowSize-1))
	{
		callback(0);
	}
	for (size_t i = windowSize; smerSize+i-1 < sequence.size(); i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		uint64_t hash = fwkmerHasher.hash();
		// palindromic k-mer prevention: those s-mers which, if picked, might lead to picking a palindromic k-mer will arbitrarily have max hash value so they won't get picked
		if (smerSize*2 > kmerSize)
		{
			if (sequence[seqPos-kmerSize/2] == complement(sequence[seqPos-kmerSize/2]) || sequence[seqPos-smerSize+1+kmerSize/2] == complement(sequence[seqPos-smerSize+1+kmerSize/2]))
			{
				hash = std::numeric_limits<uint64_t>::max();
			}
		}
		if (endSmer(hash)) hash = 0;
		// even though pop_front is used it turns out std::vector is faster than std::deque ?!
		// because pop_front is O(w), but it is only called in O(1/w) fraction of loops
		// so the performace penalty of pop_front does not scale with w!
		// and std::vector's speed in normal, non-popfront operation outweighs the slow pop_front
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i - windowSize) smerOrder.erase(smerOrder.begin());
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
		if ((std::get<0>(smerOrder.front()) == i-windowSize+1) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i))
		{
			callback(i-windowSize+1);
		}
	}
}

class ReadInfo
{
public:
	ReadInfo() = default;
	ReadName readName;
	size_t readLength;
	size_t readLengthHpc;
};

class ReadBundle
{
public:
	ReadBundle() = default;
	ReadInfo readInfo;
	SequenceCharType seq;
	SequenceLengthType poses;
	std::string rawSeq;
	std::vector<size_t> positions;
	std::vector<HashType> hashes;
};

template <typename F>
void iterateReadsMultithreaded(const std::vector<std::string>& files, const size_t numThreads, F readCallback)
{
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> sequenceQueue;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &sequenceQueue, readCallback, i]()
		{
			while (true)
			{
				std::shared_ptr<FastQ> read;
				if (!sequenceQueue.try_dequeue(read))
				{
					bool tryBreaking = readDone;
					if (!sequenceQueue.try_dequeue(read))
					{
						if (tryBreaking) return;
						std::this_thread::sleep_for(std::chrono::milliseconds(10));
						continue;
					}
				}
				assert(read != nullptr);
				ReadInfo info;
				info.readName.first = read->seq_id;
				info.readName.second = 0;
				info.readLength = read->sequence.size();
				readCallback(info, read->sequence);
			}
		});
	}
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&sequenceQueue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			bool queued = sequenceQueue.try_enqueue(ptr);
			if (queued) return;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(ptr);
				if (queued) return;
				triedSleeping += 1;
			}
			sequenceQueue.enqueue(ptr);
		});
	}
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

class ReadpartIterator
{
public:
	ReadpartIterator(const size_t kmerSize, const size_t windowSize, const ErrorMasking errorMasking, const size_t numThreads, const std::vector<std::string>& readFiles, const bool includeEndSmers, const std::string& cacheFileName);
	~ReadpartIterator();
	template <typename F>
	void iterateHashes(F callback) const
	{
		if (memoryReads.size() > 0)
		{
			iterateHashesFromMemory(callback);
		}
		else if (cacheFileName.size() > 0)
		{
			iterateHashesFromCache(callback);
		}
		else
		{
			iterateHashesFromFiles(callback);
		}
	}
	template <typename F>
	void iterateOnlyHashes(F callback) const
	{
		if (memoryReads.size() > 0)
		{
			iterateOnlyHashesFromMemory(callback);
		}
		else if (cacheFileName.size() > 0)
		{
			iterateOnlyHashesFromCache(callback);
		}
		else
		{
			iterateOnlyHashesFromFiles(callback);
		}
	}
	template <typename F>
	void iterateParts(F callback) const
	{
		if (memoryReads.size() > 0)
		{
			iteratePartsFromMemory(callback);
		}
		else if (cacheFileName.size() > 0)
		{
			iteratePartsFromCache(callback);
		}
		else
		{
			iteratePartsFromFiles(callback);
		}
	}
	template <typename F>
	void iterateHashesOfRead(const std::string& name, const std::string& seq, F callback) const
	{
		ReadInfo info;
		info.readName.first = name;
		info.readName.second = 0;
		info.readLength = seq.size();
		errorMask(info, seq, [this, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq)
		{
			iterateNonpalindromeHashes(read, seq, poses, rawSeq, callback);
		});
	}
	template <typename F>
	void iteratePartsOfRead(const std::string& name, const std::string& seq, F callback) const
	{
		ReadInfo info;
		info.readName.first = name;
		info.readName.second = 0;
		info.readLength = seq.size();
		errorMask(info, seq, callback);
	}
	void addHpcVariants(const HashType hash, const size_t offset, const std::vector<size_t>& variants);
	void clearCache();
	void clearCacheHashes();
	void setMemoryReads(const std::vector<std::pair<std::string, std::string>>& rawSeqs);
	void addMemoryRead(const std::pair<std::string, std::string>& seq);
	void setMemoryReadIterables(const std::vector<size_t>& iterables);
private:
	const size_t kmerSize;
	const size_t windowSize;
	const ErrorMasking errorMasking;
	std::vector<bool> endSmers;
	const size_t numThreads;
	const std::vector<std::string> readFiles;
	const std::string cacheFileName;
	phmap::flat_hash_map<HashType, std::vector<std::pair<size_t, std::vector<size_t>>>> hpcVariants;
	mutable size_t cacheItems;
	mutable bool cacheBuilt;
	mutable bool cache2Built;
	std::vector<ReadBundle> memoryReads;
	void collectEndSmers();
	template <typename F>
	void iteratePartsFromMemory(F callback) const
	{
		if (memoryIterables.size() > 0)
		{
			for (size_t i : memoryIterables)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].seq, memoryReads[i].poses, memoryReads[i].rawSeq);
			}
		}
		else
		{
			for (size_t i = 0; i < memoryReads.size(); i++)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].seq, memoryReads[i].poses, memoryReads[i].rawSeq);
			}
		}
	}
	template <typename F>
	void iterateHashesFromMemory(F callback) const
	{
		if (memoryIterables.size() > 0)
		{
			for (size_t i : memoryIterables)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].seq, memoryReads[i].poses, memoryReads[i].rawSeq, memoryReads[i].positions, memoryReads[i].hashes);
			}
		}
		else
		{
			for (size_t i = 0; i < memoryReads.size(); i++)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].seq, memoryReads[i].poses, memoryReads[i].rawSeq, memoryReads[i].positions, memoryReads[i].hashes);
			}
		}
	}
	template <typename F>
	void iterateOnlyHashesFromMemory(F callback) const
	{
		if (memoryIterables.size() > 0)
		{
			for (size_t i : memoryIterables)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].positions, memoryReads[i].hashes);
			}
		}
		else
		{
			for (size_t i = 0; i < memoryReads.size(); i++)
			{
				callback(memoryReads[i].readInfo, memoryReads[i].positions, memoryReads[i].hashes);
			}
		}
	}
	template <typename F>
	void iteratePartsFromCache(F callback) const
	{
		iterateHashesFromCache([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) {
			callback(read, seq, poses, rawSeq);
		});
	}
	template <typename F>
	void buildSecondCacheAndIterateHashes(F callback) const
	{
		assert(cacheFileName.size() > 0);
		std::atomic<bool> readDone;
		readDone = false;
		std::vector<std::thread> threads;
		moodycamel::ConcurrentQueue<std::shared_ptr<ReadBundle>> sequenceQueue;
		std::ofstream cachePart2 { cacheFileName + "2", std::ios::binary };
		if (!cachePart2.good())
		{
			std::cerr << "Could not build cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		for (size_t i = 0; i < numThreads; i++)
		{
			threads.emplace_back([this, &readDone, &sequenceQueue, &cachePart2, callback]()
			{
				while (true)
				{
					std::shared_ptr<ReadBundle> read;
					if (!sequenceQueue.try_dequeue(read))
					{
						bool tryBreaking = readDone;
						if (!sequenceQueue.try_dequeue(read))
						{
							if (tryBreaking) return;
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
							continue;
						}
					}
					assert(read != nullptr);
					callback(read->readInfo, read->seq, read->poses, read->rawSeq, read->positions, read->hashes);
				}
			});
		}
		std::ifstream cache { cacheFileName, std::ios::binary };
		if (!cache.good())
		{
			std::cerr << "Could not read cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		size_t itemsRead = 0;
		while (cache.good())
		{
			cache.peek();
			if (!cache.good()) break;
			std::shared_ptr<ReadBundle> readInfo { new ReadBundle };
			Serializer::read(cache, readInfo->readInfo.readName.first);
			Serializer::read(cache, readInfo->readInfo.readName.second);
			Serializer::readMostlyTwobits(cache, readInfo->seq);
			Serializer::readMonotoneIncreasing(cache, readInfo->poses);
			Serializer::readTwobits(cache, readInfo->rawSeq);
			assert(readInfo->poses.size() == 0 || readInfo->poses[0] == 0);
			assert(readInfo->poses.size() == 0 || readInfo->poses.back() == readInfo->rawSeq.size());
			iterateNonpalindromeHashes(readInfo->readInfo, readInfo->seq, readInfo->poses, readInfo->rawSeq, [&readInfo](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, std::vector<HashType>& hashes)
			{
				readInfo->positions = positions;
				readInfo->hashes = hashes;
			});
			Serializer::write(cachePart2, readInfo->readInfo.readName.first);
			Serializer::write(cachePart2, readInfo->readInfo.readName.second);
			Serializer::write(cachePart2, readInfo->rawSeq.size());
			Serializer::write(cachePart2, readInfo->seq.size());
			Serializer::writeMonotoneIncreasing(cachePart2, readInfo->positions);
			Serializer::write(cachePart2, readInfo->hashes);
			itemsRead += 1;
			bool queued = sequenceQueue.try_enqueue(readInfo);
			if (queued) continue;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(readInfo);
				if (queued) break;
				triedSleeping += 1;
			}
			if (queued) continue;
			sequenceQueue.enqueue(readInfo);
		}
		if (itemsRead != cacheItems)
		{
			std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
			std::abort();
		}
		readDone = true;
		for (size_t i = 0; i < threads.size(); i++)
		{
			threads[i].join();
		}
		if (!cachePart2.good())
		{
			std::cerr << "Could not build cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		cache2Built = true;
	}
	template <typename F>
	void buildCacheAndIterateHashes(F callback) const
	{
		std::cerr << "Building sequence cache" << std::endl;
		assert(cacheFileName.size() > 0);
		std::ofstream cache { cacheFileName, std::ios::binary };
		std::ofstream cachePart2 { cacheFileName + "2", std::ios::binary };
		if (!cache.good() || !cachePart2.good())
		{
			std::cerr << "Could not build cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		std::mutex writeMutex;
		cacheItems = 0;
		iterateHashesFromFiles([this, &cache, &cachePart2, &writeMutex, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
		{
			callback(read, seq, poses, rawSeq, positions, hashes);
			std::lock_guard<std::mutex> lock { writeMutex };
			Serializer::write(cache, read.readName.first);
			Serializer::write(cache, read.readName.second);
			Serializer::writeMostlyTwobits(cache, seq);
			Serializer::writeMonotoneIncreasing(cache, poses);
			Serializer::writeTwobits(cache, rawSeq);
			Serializer::write(cachePart2, read.readName.first);
			Serializer::write(cachePart2, read.readName.second);
			Serializer::write(cachePart2, rawSeq.size());
			Serializer::write(cachePart2, seq.size());
			Serializer::writeMonotoneIncreasing(cachePart2, positions);
			Serializer::write(cachePart2, hashes);
			cacheItems += 1;
		});
		if (!cache.good() || !cachePart2.good())
		{
			std::cerr << "Could not build cache. Try running without sequence cache." << std::endl;
			cache.close();
			cachePart2.close();
			remove(cacheFileName.c_str());
			std::string cache2name = cacheFileName;
			cache2name += '2';
			remove(cache2name.c_str());
			std::abort();
		}
		std::cerr << "Stored " << cacheItems << " sequences in the cache" << std::endl;
		cacheBuilt = true;
		cache2Built = true;
	}
	template <typename F>
	void iterateOnlyHashesFromCache(F callback) const
	{
		if (!cacheBuilt)
		{
			buildCacheAndIterateHashes([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) {
				callback(read, positions, hashes);
			});
			assert(cacheBuilt);
			return;
		}
		if (cacheBuilt && !cache2Built)
		{
			buildSecondCacheAndIterateHashes([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) {
				callback(read, positions, hashes);
			});
			assert(cache2Built);
			return;
		}
		std::cerr << "Reading sequences from cache" << std::endl;
		std::atomic<bool> readDone;
		readDone = false;
		std::vector<std::thread> threads;
		moodycamel::ConcurrentQueue<std::shared_ptr<ReadBundle>> sequenceQueue;
		for (size_t i = 0; i < numThreads; i++)
		{
			threads.emplace_back([&readDone, &sequenceQueue, callback]()
			{
				while (true)
				{
					std::shared_ptr<ReadBundle> read;
					if (!sequenceQueue.try_dequeue(read))
					{
						bool tryBreaking = readDone;
						if (!sequenceQueue.try_dequeue(read))
						{
							if (tryBreaking) return;
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
							continue;
						}
					}
					assert(read != nullptr);
					callback(read->readInfo, read->positions, read->hashes);
				}
			});
		}
		std::ifstream cachePart2 { cacheFileName + "2", std::ios::binary };
		if (!cachePart2.good())
		{
			std::cerr << "Could not read cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		size_t itemsRead = 0;
		while (cachePart2.good())
		{
			cachePart2.peek();
			if (!cachePart2.good()) break;
			std::shared_ptr<ReadBundle> readInfo { new ReadBundle };
			Serializer::read(cachePart2, readInfo->readInfo.readName.first);
			Serializer::read(cachePart2, readInfo->readInfo.readName.second);
			Serializer::read(cachePart2, readInfo->readInfo.readLength);
			Serializer::read(cachePart2, readInfo->readInfo.readLengthHpc);
			Serializer::readMonotoneIncreasing(cachePart2, readInfo->positions);
			Serializer::read(cachePart2, readInfo->hashes);
			itemsRead += 1;
			bool queued = sequenceQueue.try_enqueue(readInfo);
			if (queued) continue;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(readInfo);
				if (queued) break;
				triedSleeping += 1;
			}
			if (queued) continue;
			sequenceQueue.enqueue(readInfo);
		}
		if (itemsRead != cacheItems)
		{
			std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
			std::abort();
		}
		readDone = true;
		for (size_t i = 0; i < threads.size(); i++)
		{
			threads[i].join();
		}
	}
	template <typename F>
	void iterateHashesFromCache(F callback) const
	{
		if (!cacheBuilt)
		{
			buildCacheAndIterateHashes(callback);
			assert(cacheBuilt);
			return;
		}
		if (cacheBuilt && !cache2Built)
		{
			buildSecondCacheAndIterateHashes(callback);
			assert(cache2Built);
			return;
		}
		std::cerr << "Reading sequences from cache" << std::endl;
		std::atomic<bool> readDone;
		readDone = false;
		std::vector<std::thread> threads;
		moodycamel::ConcurrentQueue<std::shared_ptr<ReadBundle>> sequenceQueue;
		for (size_t i = 0; i < numThreads; i++)
		{
			threads.emplace_back([&readDone, &sequenceQueue, callback]()
			{
				while (true)
				{
					std::shared_ptr<ReadBundle> read;
					if (!sequenceQueue.try_dequeue(read))
					{
						bool tryBreaking = readDone;
						if (!sequenceQueue.try_dequeue(read))
						{
							if (tryBreaking) return;
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
							continue;
						}
					}
					assert(read != nullptr);
					callback(read->readInfo, read->seq, read->poses, read->rawSeq, read->positions, read->hashes);
				}
			});
		}
		std::ifstream cache { cacheFileName, std::ios::binary };
		std::ifstream cachePart2 { cacheFileName + "2", std::ios::binary };
		if (!cache.good() || !cachePart2.good())
		{
			std::cerr << "Could not read cache. Try running without sequence cache." << std::endl;
			std::abort();
		}
		size_t itemsRead = 0;
		while (cache.good())
		{
			cache.peek();
			if (!cache.good()) break;
			std::shared_ptr<ReadBundle> readInfo { new ReadBundle };
			Serializer::read(cache, readInfo->readInfo.readName.first);
			Serializer::read(cache, readInfo->readInfo.readName.second);
			Serializer::readMostlyTwobits(cache, readInfo->seq);
			Serializer::readMonotoneIncreasing(cache, readInfo->poses);
			Serializer::readTwobits(cache, readInfo->rawSeq);
			std::string tmp;
			Serializer::read(cachePart2, tmp);
			if (tmp != readInfo->readInfo.readName.first)
			{
				std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
				std::abort();
			}
			size_t tmp2;
			Serializer::read(cachePart2, tmp2);
			if (tmp2 != readInfo->readInfo.readName.second)
			{
				std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
				std::abort();
			}
			Serializer::read(cachePart2, tmp2);
			if (tmp2 != readInfo->rawSeq.size())
			{
				std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
				std::abort();
			}
			Serializer::read(cachePart2, tmp2);
			if (tmp2 != readInfo->seq.size())
			{
				std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
				std::abort();
			}
			Serializer::readMonotoneIncreasing(cachePart2, readInfo->positions);
			Serializer::read(cachePart2, readInfo->hashes);
			assert(readInfo->poses.size() == 0 || readInfo->poses[0] == 0);
			assert(readInfo->poses.size() == 0 || readInfo->poses.back() == readInfo->rawSeq.size());
			itemsRead += 1;
			bool queued = sequenceQueue.try_enqueue(readInfo);
			if (queued) continue;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(readInfo);
				if (queued) break;
				triedSleeping += 1;
			}
			if (queued) continue;
			sequenceQueue.enqueue(readInfo);
		}
		if (itemsRead != cacheItems)
		{
			std::cerr << "Sequence cache has been corrupted. Try re-running, or running without sequence cache." << std::endl;
			std::abort();
		}
		readDone = true;
		for (size_t i = 0; i < threads.size(); i++)
		{
			threads[i].join();
		}
	}
	template <typename F>
	void iterateOnlyHashesFromFiles(F callback) const
	{
		iterateHashesFromFiles([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) {
			callback(read, positions, hashes);
		});
	}
	template <typename F>
	void iterateHashesFromFiles(F callback) const
	{
		if (hpcVariants.size() == 0)
		{
			iterateHashesFromFilesInternal(callback);
		}
		else
		{
			iterateHashesFromFilesInternal([this, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, std::vector<HashType>& hashes)
			{
				bool variant = false;
				for (size_t i = 0; i < hashes.size(); i++)
				{
					HashType fwHash = hashes[i];
					HashType bwHash = (fwHash << 64) + (fwHash >> 64);
					if (hpcVariants.count(fwHash) == 0 && hpcVariants.count(bwHash) == 0) continue;
					std::vector<size_t> firstVariantLengths;
					std::vector<size_t> secondVariantLengths;
					if (hpcVariants.count(fwHash) == 1)
					{
						assert(hpcVariants.count(bwHash) == 0);
						for (const auto& pair : hpcVariants.at(fwHash))
						{
							size_t lengthHere = poses[positions[i] + pair.first + 1] - poses[positions[i] + pair.first];
							assert(lengthHere <= pair.second.back());
							size_t lengthCluster = *std::lower_bound(pair.second.begin(), pair.second.end(), lengthHere);
							if (pair.first <= kmerSize/2)
							{
								firstVariantLengths.push_back(lengthCluster);
							}
							if (pair.first >= kmerSize/2)
							{
								secondVariantLengths.push_back(lengthCluster);
							}
						}
						std::reverse(secondVariantLengths.begin(), secondVariantLengths.end());
					}
					else
					{
						assert(hpcVariants.count(fwHash) == 0);
						for (const auto& pair : hpcVariants.at(bwHash))
						{
							size_t lengthHere = poses[positions[i] + kmerSize - pair.first] - poses[positions[i] + kmerSize - pair.first - 1];
							assert(lengthHere <= pair.second.back());
							size_t lengthCluster = *std::lower_bound(pair.second.begin(), pair.second.end(), lengthHere);
							if (pair.first >= kmerSize/2)
							{
								firstVariantLengths.push_back(lengthCluster);
							}
							if (pair.first <= kmerSize/2)
							{
								secondVariantLengths.push_back(lengthCluster);
							}
						}
						std::reverse(firstVariantLengths.begin(), firstVariantLengths.end());
					}
					if (firstVariantLengths.size() > 0 || secondVariantLengths.size() > 0) variant = true;
					if (variant)
					{
						uint64_t firstVariantHash = std::hash<const std::vector<size_t>&>{}(firstVariantLengths);
						uint64_t secondVariantHash = std::hash<const std::vector<size_t>&>{}(secondVariantLengths);
						hashes[i] = hashes[i] ^ (HashType)firstVariantHash ^ (((HashType)secondVariantHash) << 64);
					}
				}
				callback(read, seq, poses, rawSeq, positions, hashes);
			});
		}
	}
	template <typename F>
	void iterateNonpalindromeHashes(const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, F callback) const
	{
		iterateKmers(read, seq, rawSeq, poses, [this, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positionsWithPalindromes) {
			SequenceCharType revSeq = revCompRLE(seq);
			std::vector<size_t> positions;
			std::vector<HashType> hashes;
			for (size_t i = 0; i < positionsWithPalindromes.size(); i++)
			{
				const auto pos = positionsWithPalindromes[i];
				VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				VectorView<CharType> revMinimizerSequence { revSeq, revPos, revPos + kmerSize };
				HashType fwHash = hash(minimizerSequence, revMinimizerSequence);
				HashType bwHash = (fwHash << 64) + (fwHash >> 64);
				if (fwHash == bwHash)
				{
					bool palindrome = true;
					for (size_t j = 0; j < kmerSize/2; j++)
					{
						if (minimizerSequence[j] != complement(minimizerSequence[kmerSize-1-j]))
						{
							palindrome = false;
							break;
						}
					}
					if (palindrome)
					{
						if (i == 0 || i == positionsWithPalindromes.size()-1 || positionsWithPalindromes[i+1] - positionsWithPalindromes[i-1] < kmerSize)
						{
							// palindromic k-mer, but it can be safely dropped without creating gaps
							continue;
						}
						else
						{
							std::cerr << "The genome has a palindromic k-mer. Cannot build a graph. Try running with a different -w" << std::endl;
							std::cerr << "Example read around the palindromic k-mer: " << read.readName.first << std::endl;
							std::abort();
						}
					}
					else
					{
						std::cerr << "Unhashable k-mer around read: " << read.readName.first << std::endl;
						std::abort();
					}
				}
				positions.push_back(positionsWithPalindromes[i]);
				hashes.push_back(fwHash);
			}
			for (size_t i = 1; i < positions.size(); i++)
			{
				if (positions[i] >= positions[i-1]+kmerSize)
				{
					std::cerr << "The genome has too many palindromic k-mers. Cannot build a graph. Try running with a different -w or different -k" << std::endl;
					std::cerr << "Example read with palindromic k-mers: " << read.readName.first << std::endl;
					std::abort();
				}
			}
			callback(read, seq, poses, rawSeq, positions, hashes);
		});
	}
	template <typename F>
	void iterateHashesFromFilesInternal(F callback) const
	{
		iteratePartsFromFiles([this, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq) {
			if (seq.size() < kmerSize) return;
			iterateNonpalindromeHashes(read, seq, poses, rawSeq, callback);
		});
	}
	template <typename F>
	void errorMask(ReadInfo& read, const std::string& rawSeq, F callback) const
	{
		if (errorMasking == ErrorMasking::Hpc)
		{
			iterateRLE(read, rawSeq, callback);
		}
		else if (errorMasking == ErrorMasking::Collapse)
		{
			iterateCollapse(read, rawSeq, callback);
		}
		else if (errorMasking == ErrorMasking::Dinuc)
		{
			iterateDinuc(read, rawSeq, callback);
		}
		else if (errorMasking == ErrorMasking::Microsatellite)
		{
			iterateMicrosatellite(read, rawSeq, callback);
		}
		else if (errorMasking == ErrorMasking::CollapseDinuc)
		{
			iterateCollapseDinuc(read, rawSeq, callback);
		}
		else if (errorMasking == ErrorMasking::CollapseMicrosatellite)
		{
			iterateCollapseMicrosatellite(read, rawSeq, callback);
		}
		else
		{
			assert(errorMasking == ErrorMasking::No);
			iterateNoRLE(read, rawSeq, callback);
		}
	}
	template <typename F>
	void iteratePartsFromFiles(F callback) const
	{
		iterateReadsMultithreaded(readFiles, numThreads, [this, callback](ReadInfo& read, const std::string& rawSeq)
		{
			if (rawSeq.size() < 32) return;
			errorMask(read, rawSeq, callback);
		});
	}
	template <typename F>
	void iterateMicrosatellite(ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateRLE(read, seq, [callback](ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			if (seq.size() < 32) return;
			auto pieces = multiRLECompress(seq, poses, 6);
			read.readLengthHpc = pieces.first.size();
			callback(read, pieces.first, pieces.second, raw);
		});
	}
	template <typename F>
	void iterateCollapseMicrosatellite(ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateCollapse(read, seq, [callback](ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			if (seq.size() < 32) return;
			auto pieces = multiRLECompress(seq, poses, 6);
			read.readLengthHpc = pieces.first.size();
			callback(read, pieces.first, pieces.second, raw);
		});
	}
	template <typename F>
	void iterateCollapse(ReadInfo& read, const std::string& seq, F callback) const
	{
		if (seq.size() == 0) return;
		std::string collapseSeq;
		collapseSeq.reserve(seq.size());
		collapseSeq.push_back(seq[0]);
		for (size_t i = 1; i < seq.size(); i++)
		{
			if (seq[i] != seq[i-1]) collapseSeq.push_back(seq[i]);
		}
		iterateRLE(read, collapseSeq, callback);
	}
	template <typename F>
	void iterateDinuc(ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateRLE(read, seq, [callback](ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			if (seq.size() < 32) return;
			auto pieces = multiRLECompress(seq, poses, 2);
			read.readLengthHpc = pieces.first.size();
			callback(read, pieces.first, pieces.second, raw);
		});
	}
	template <typename F>
	void iterateCollapseDinuc(ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateCollapse(read, seq, [callback](ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			if (seq.size() < 32) return;
			auto pieces = multiRLECompress(seq, poses, 2);
			read.readLengthHpc = pieces.first.size();
			callback(read, pieces.first, pieces.second, raw);
		});
	}
	template <typename F>
	void iterateRLE(ReadInfo& read, const std::string& seq, F callback) const
	{
		SequenceCharType currentSeq;
		SequenceLengthType currentPos;
		currentSeq.reserve(seq.size());
		currentPos.reserve(seq.size());
		size_t i = 0;
		assert(currentSeq.size() == 0);
		assert(currentPos.size() == 0);
		while (i < seq.size() && seq[i] != 'a' && seq[i] != 'A' && seq[i] != 'c' && seq[i] != 'C' && seq[i] != 'g' && seq[i] != 'G' && seq[i] != 't' && seq[i] != 'T') i += 1;
		if (i == seq.size()) return;
		size_t lastStart = i;
		switch(seq[i])
		{
			case 'a':
			case 'A':
				currentSeq.push_back(0);
				currentPos.push_back(i);
				break;
			case 'c':
			case 'C':
				currentSeq.push_back(1);
				currentPos.push_back(i);
				break;
			case 'g':
			case 'G':
				currentSeq.push_back(2);
				currentPos.push_back(i);
				break;
			case 't':
			case 'T':
				currentSeq.push_back(3);
				currentPos.push_back(i);
				break;
			default:
				assert(false);
		}
		i += 1;
		for (; i < seq.size(); i++)
		{
			switch(seq[i])
			{
				case 'a':
				case 'A':
					if (currentSeq.back() != 0)
					{
						currentSeq.push_back(0);
						currentPos.push_back(i);
					}
					break;
				case 'c':
				case 'C':
					if (currentSeq.back() != 1)
					{
						currentSeq.push_back(1);
						currentPos.push_back(i);
					}
					break;
				case 'g':
				case 'G':
					if (currentSeq.back() != 2)
					{
						currentSeq.push_back(2);
						currentPos.push_back(i);
					}
					break;
				case 't':
				case 'T':
					if (currentSeq.back() != 3)
					{
						currentSeq.push_back(3);
						currentPos.push_back(i);
					}
					break;
				default:
					currentPos.push_back(i);
					read.readLengthHpc = currentSeq.size();
					read.readName.second = lastStart;
					callback(read, currentSeq, currentPos, seq);
					currentSeq.clear();
					currentPos.clear();
					while (i < seq.size() && seq[i] != 'a' && seq[i] != 'A' && seq[i] != 'c' && seq[i] != 'C' && seq[i] != 'g' && seq[i] != 'G' && seq[i] != 't' && seq[i] != 'T') i += 1;
					lastStart = i;
					if (i == seq.size()) return;
					switch(seq[i])
					{
						case 'a':
						case 'A':
							currentSeq.push_back(0);
							currentPos.push_back(i);
							break;
						case 'c':
						case 'C':
							currentSeq.push_back(1);
							currentPos.push_back(i);
							break;
						case 'g':
						case 'G':
							currentSeq.push_back(2);
							currentPos.push_back(i);
							break;
						case 't':
						case 'T':
							currentSeq.push_back(3);
							currentPos.push_back(i);
							break;
						default:
							assert(false);
					}
			}
		}
		if (currentSeq.size() > 0)
		{
			currentPos.push_back(seq.size());
			read.readLengthHpc = currentSeq.size();
			read.readName.second = lastStart;
			callback(read, currentSeq, currentPos, seq);
		}
	}
	template <typename F>
	void iterateNoRLE(ReadInfo& read, const std::string& seq, F callback) const
	{
		SequenceCharType currentSeq;
		SequenceLengthType currentPos;
		currentSeq.reserve(seq.size());
		currentPos.reserve(seq.size());
		size_t i = 0;
		size_t lastStart = 0;
		assert(currentSeq.size() == 0);
		assert(currentPos.size() == 0);
		for (; i < seq.size(); i++)
		{
			switch(seq[i])
			{
				case 'a':
				case 'A':
					currentSeq.push_back(0);
					currentPos.push_back(i);
					break;
				case 'c':
				case 'C':
					currentSeq.push_back(1);
					currentPos.push_back(i);
					break;
				case 'g':
				case 'G':
					currentSeq.push_back(2);
					currentPos.push_back(i);
					break;
				case 't':
				case 'T':
					currentSeq.push_back(3);
					currentPos.push_back(i);
					break;
				default:
					if (currentSeq.size() > 0)
					{
						currentPos.push_back(i);
						read.readLengthHpc = currentSeq.size();
						read.readName.second = lastStart;
						callback(read, currentSeq, currentPos, seq);
					}
					currentSeq.clear();
					currentPos.clear();
					lastStart = i;
			}
		}
		if (currentSeq.size() > 0)
		{
			currentPos.push_back(i);
			read.readLengthHpc = currentSeq.size();
			read.readName.second = lastStart;
			callback(read, currentSeq, currentPos, seq);
		}
	}
	template <typename F>
	void iterateKmers(const ReadInfo& read, const SequenceCharType& seq, const std::string& rawSeq, const SequenceLengthType& poses, F callback) const
	{
		if (seq.size() < kmerSize) return;
		// keep the same smerOrder to reduce mallocs which destroy multithreading performance
		thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
		smerOrder.resize(0);
		std::vector<size_t> positions;
		if (endSmers.size() > 0)
		{
			findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [this](uint64_t hash) { return endSmers[hash % endSmers.size()]; }, [this, &positions](size_t pos)
			{
				assert(positions.size() == 0 || pos > positions.back());
				assert(positions.size() == 0 || pos - positions.back() <= windowSize);
				positions.push_back(pos);
			});
		}
		else
		{
			findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [](uint64_t hash) { return false; }, [this, &positions](size_t pos)
			{
				// assert(positions.size() == 0 || pos > positions.back());
				// assert(positions.size() == 0 || pos - positions.back() <= windowSize);
				positions.push_back(pos);
			});
		}
		callback(read, seq, poses, rawSeq, positions);
	}
};

}

#endif
