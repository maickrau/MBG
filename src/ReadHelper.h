#ifndef ReadHelper_h
#define ReadHelper_h

#include <vector>
#include <string>
#include <fstream>
#include <atomic>
#include <thread>
#include <iostream>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include "fastqloader.h"
#include "FastHasher.h"
#include "ErrorMaskHelper.h"
#include "Serializer.h"

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
		size_t hash = fwkmerHasher.hash();
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
		size_t hash = fwkmerHasher.hash();
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
	std::string readName;
};

class ReadBundle
{
public:
	ReadInfo readInfo;
	SequenceCharType seq;
	SequenceLengthType poses;
	std::string rawSeq;
	std::vector<size_t> positions;
	std::vector<HashType> hashes;
};

class ReadpartIterator
{
public:
	ReadpartIterator(const size_t kmerSize, const size_t windowSize, const ErrorMasking errorMasking, const size_t numThreads, const std::vector<std::string>& readFiles, const bool includeEndSmers, const std::string& cacheFileName);
	~ReadpartIterator();
	template <typename F>
	void iterateHashes(F callback) const
	{
		if (cacheFileName.size() > 0)
		{
			iterateHashesFromCache(callback);
		}
		else
		{
			iterateHashesFromFiles(callback);
		}
	}
	template <typename F>
	void iterateParts(F callback) const
	{
		if (cacheFileName.size() > 0)
		{
			iteratePartsFromCache(callback);
		}
		else
		{
			iteratePartsFromFiles(callback);
		}
	}
private:
	const size_t kmerSize;
	const size_t windowSize;
	const ErrorMasking errorMasking;
	std::vector<bool> endSmers;
	const size_t numThreads;
	const std::vector<std::string> readFiles;
	const std::string cacheFileName;
	size_t cacheItems;
	void collectEndSmers();
	void buildCache();
	template <typename F>
	void iteratePartsFromCache(F callback) const
	{
		iterateHashesFromCache([callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes) {
			callback(read, seq, poses, rawSeq);
		});
	}
	template <typename F>
	void iterateHashesFromCache(F callback) const
	{
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
			Serializer::read(cache, readInfo->readInfo.readName);
			Serializer::read(cache, readInfo->seq);
			Serializer::read(cache, readInfo->poses);
			Serializer::read(cache, readInfo->rawSeq);
			Serializer::read(cache, readInfo->positions);
			Serializer::read(cache, readInfo->hashes);
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
	void iterateHashesFromFiles(F callback) const
	{
		iteratePartsFromFiles([this, callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq) {
			if (seq.size() < kmerSize) return;
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
								std::cerr << "Example read around the palindromic k-mer: " << read.readName << std::endl;
								std::abort();
							}
						}
						else
						{
							std::cerr << "Unhashable k-mer around read: " << read.readName << std::endl;
							std::abort();
						}
					}
					positions.push_back(positionsWithPalindromes[i]);
					hashes.push_back(fwHash);
				}
				callback(read, seq, poses, rawSeq, positions, hashes);
			});
		});
	}
	template <typename F>
	void iteratePartsFromFiles(F callback) const
	{
		iterateReadsMultithreaded(readFiles, numThreads, [this, callback](const ReadInfo& read, const std::string& rawSeq)
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
		});
	}
	template <typename F>
	void iterateMicrosatellite(const ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateRLE(read, seq, [callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			auto pieces = multiRLECompress(seq, poses, 6);
			for (const auto& pair : pieces)
			{
				callback(read, pair.first, pair.second, raw);
			}
		});
	}
	template <typename F>
	void iterateCollapseMicrosatellite(const ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateCollapse(read, seq, [callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			auto pieces = multiRLECompress(seq, poses, 6);
			for (const auto& pair : pieces)
			{
				callback(read, pair.first, pair.second, raw);
			}
		});
	}
	template <typename F>
	void iterateCollapse(const ReadInfo& read, const std::string& seq, F callback) const
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
	void iterateDinuc(const ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateRLE(read, seq, [callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			auto pieces = multiRLECompress(seq, poses, 2);
			for (const auto& pair : pieces)
			{
				callback(read, pair.first, pair.second, raw);
			}
		});
	}
	template <typename F>
	void iterateCollapseDinuc(const ReadInfo& read, const std::string& seq, F callback) const
	{
		iterateCollapse(read, seq, [callback](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& raw)
		{
			auto pieces = multiRLECompress(seq, poses, 2);
			for (const auto& pair : pieces)
			{
				callback(read, pair.first, pair.second, raw);
			}
		});
	}
	template <typename F>
	void iterateRLE(const ReadInfo& read, const std::string& seq, F callback) const
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
					callback(read, currentSeq, currentPos, seq);
					currentSeq.clear();
					currentPos.clear();
					while (i < seq.size() && seq[i] != 'a' && seq[i] != 'A' && seq[i] != 'c' && seq[i] != 'C' && seq[i] != 'g' && seq[i] != 'G' && seq[i] != 't' && seq[i] != 'T') i += 1;
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
			callback(read, currentSeq, currentPos, seq);
		}
	}
	template <typename F>
	void iterateNoRLE(const ReadInfo& read, const std::string& seq, F callback) const
	{
		SequenceCharType currentSeq;
		SequenceLengthType currentPos;
		currentSeq.reserve(seq.size());
		currentPos.reserve(seq.size());
		size_t i = 0;
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
						callback(read, currentSeq, currentPos, seq);
					}
					currentSeq.clear();
					currentPos.clear();
			}
		}
		if (currentSeq.size() > 0)
		{
			currentPos.push_back(i);
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
	template <typename F>
	void iterateReadsMultithreaded(const std::vector<std::string>& files, const size_t numThreads, F readCallback) const
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
					info.readName = read->seq_id;
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

};

#endif
