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

enum ErrorMasking
{
	No,
	Hpc,
};

template <typename F, typename EdgeCheckFunction>
uint64_t findSyncmerPositions(const std::string& sequence, size_t kmerSize, size_t smerSize, std::vector<std::tuple<size_t, uint64_t>>& smerOrder, EdgeCheckFunction endSmer, F callback)
{
	if (sequence.size() < kmerSize) return 0;
	assert(smerSize <= kmerSize);
	size_t windowSize = kmerSize - smerSize + 1;
	assert(windowSize >= 1);
	uint64_t minHash = std::numeric_limits<uint64_t>::max();
	FastHasher fwkmerHasher { smerSize };
	for (size_t i = 0; i < smerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	minHash = std::min(minHash, fwkmerHasher.hash());
	auto thisHash = fwkmerHasher.hash();
	if (endSmer(thisHash)) thisHash = 0;
	smerOrder.emplace_back(0, thisHash);
	for (size_t i = 1; i < windowSize; i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		size_t hash = fwkmerHasher.hash();
		minHash = std::min(minHash, hash);
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
		minHash = std::min(minHash, hash);
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
	return minHash;
}

class ReadpartIterator
{
public:
	ReadpartIterator(const size_t kmerSize, const size_t windowSize, const ErrorMasking errorMasking) :
	kmerSize(kmerSize),
	windowSize(windowSize),
	errorMasking(errorMasking)
	{
	}
	template <typename F>
	void iteratePartKmers(const FastQ& read, F callback) const
	{
		if (read.sequence.size() < kmerSize) return;
		iterateParts(read, [this, callback](const std::string& seq, const std::vector<uint16_t>& lens) {
			iterateKmers(seq, lens, callback);
		});
	}
	template <typename F>
	void iterateParts(const FastQ& read, F callback) const
	{
		if (errorMasking == ErrorMasking::Hpc)
		{
			iterateRLE(read.sequence, callback);
		}
		else
		{
			assert(errorMasking == ErrorMasking::No);
			iterateNoRLE(read.sequence, callback);
		}
	}
	const size_t kmerSize;
	const size_t windowSize;
	const ErrorMasking errorMasking;
	std::vector<bool> endSmers;
private:
	template <typename F>
	void iterateRLE(const std::string& seq, F callback) const
	{
		std::string currentSeq;
		std::vector<uint16_t> currentLens;
		currentSeq.reserve(seq.size());
		currentLens.reserve(seq.size());
		size_t i = 0;
		assert(currentSeq.size() == 0);
		assert(currentLens.size() == 0);
		while (i < seq.size() && seq[i] != 'a' && seq[i] != 'A' && seq[i] != 'c' && seq[i] != 'C' && seq[i] != 'g' && seq[i] != 'G' && seq[i] != 't' && seq[i] != 'T') i += 1;
		if (i == seq.size()) return;
		switch(seq[i])
		{
			case 'a':
			case 'A':
				currentSeq.push_back(1);
				currentLens.push_back(1);
				break;
			case 'c':
			case 'C':
				currentSeq.push_back(2);
				currentLens.push_back(1);
				break;
			case 'g':
			case 'G':
				currentSeq.push_back(3);
				currentLens.push_back(1);
				break;
			case 't':
			case 'T':
				currentSeq.push_back(4);
				currentLens.push_back(1);
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
					if (currentSeq.back() != 1)
					{
						currentSeq.push_back(1);
						currentLens.push_back(1);
					}
					else
					{
						currentLens.back() += 1;
					}
					break;
				case 'c':
				case 'C':
					if (currentSeq.back() != 2)
					{
						currentSeq.push_back(2);
						currentLens.push_back(1);
					}
					else
					{
						currentLens.back() += 1;
					}
					break;
				case 'g':
				case 'G':
					if (currentSeq.back() != 3)
					{
						currentSeq.push_back(3);
						currentLens.push_back(1);
					}
					else
					{
						currentLens.back() += 1;
					}
					break;
				case 't':
				case 'T':
					if (currentSeq.back() != 4)
					{
						currentSeq.push_back(4);
						currentLens.push_back(1);
					}
					else
					{
						currentLens.back() += 1;
					}
					break;
				default:
					callback(currentSeq, currentLens);
					currentSeq.clear();
					currentLens.clear();
					while (i < seq.size() && seq[i] != 'a' && seq[i] != 'A' && seq[i] != 'c' && seq[i] != 'C' && seq[i] != 'g' && seq[i] != 'G' && seq[i] != 't' && seq[i] != 'T') i += 1;
					if (i == seq.size()) return;
					switch(seq[i])
					{
						case 'a':
						case 'A':
							currentSeq.push_back(1);
							currentLens.push_back(1);
							break;
						case 'c':
						case 'C':
							currentSeq.push_back(2);
							currentLens.push_back(1);
							break;
						case 'g':
						case 'G':
							currentSeq.push_back(3);
							currentLens.push_back(1);
							break;
						case 't':
						case 'T':
							currentSeq.push_back(4);
							currentLens.push_back(1);
							break;
						default:
							assert(false);
					}
			}
		}
		if (currentSeq.size() > 0)
		{
			callback(currentSeq, currentLens);
		}
	}
	template <typename F>
	void iterateNoRLE(const std::string& seq, F callback) const
	{
		std::string currentSeq;
		std::vector<uint16_t> currentLens;
		currentSeq.reserve(seq.size());
		currentLens.reserve(seq.size());
		size_t i = 0;
		assert(currentSeq.size() == 0);
		assert(currentLens.size() == 0);
		for (; i < seq.size(); i++)
		{
			switch(seq[i])
			{
				case 'a':
				case 'A':
					currentSeq.push_back(1);
					currentLens.push_back(1);
					break;
				case 'c':
				case 'C':
					currentSeq.push_back(2);
					currentLens.push_back(1);
					break;
				case 'g':
				case 'G':
					currentSeq.push_back(3);
					currentLens.push_back(1);
					break;
				case 't':
				case 'T':
					currentSeq.push_back(4);
					currentLens.push_back(1);
					break;
				default:
					if (currentSeq.size() > 0) callback(currentSeq, currentLens);
					currentSeq.clear();
					currentLens.clear();
			}
		}
		if (currentSeq.size() > 0)
		{
			callback(currentSeq, currentLens);
		}
	}
	template <typename F>
	void iterateKmers(const std::string& seq, const std::vector<uint16_t>& lens, F callback) const
	{
		if (seq.size() < kmerSize) return;
		// keep the same smerOrder to reduce mallocs which destroy multithreading performance
		thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
		smerOrder.resize(0);
		std::vector<size_t> positions;
		uint64_t minHash;
		if (endSmers.size() > 0)
		{
			minHash = findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [this](uint64_t hash) { return endSmers[hash % endSmers.size()]; }, [this, &positions](size_t pos)
			{
				assert(positions.size() == 0 || pos > positions.back());
				assert(positions.size() == 0 || pos - positions.back() <= windowSize);
				positions.push_back(pos);
			});
		}
		else
		{
			minHash = findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [](uint64_t hash) { return false; }, [this, &positions](size_t pos)
			{
				// assert(positions.size() == 0 || pos > positions.back());
				// assert(positions.size() == 0 || pos - positions.back() <= windowSize);
				positions.push_back(pos);
			});
		}
		callback(seq, lens, minHash, positions);
	}
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
				readCallback(i, *read);
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

#endif
