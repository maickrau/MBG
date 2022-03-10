#include <mutex>
#include <cstdio>
#include "ReadHelper.h"
#include "FastHasher.h"

ReadpartIterator::ReadpartIterator(const size_t kmerSize, const size_t windowSize, const ErrorMasking errorMasking, const size_t numThreads, const std::vector<std::string>& readFiles, const bool includeEndSmers, const std::string& cacheFileName) :
	kmerSize(kmerSize),
	windowSize(windowSize),
	errorMasking(errorMasking),
	numThreads(numThreads),
	readFiles(readFiles),
	cacheFileName(cacheFileName),
	cacheItems(0),
	cacheBuilt(false)
{
	if (includeEndSmers)
	{
		std::cerr << "Collecting end k-mers" << std::endl;
		// 2^32 arbitrarily, should be big enough for a human genome
		// a 30x coverage error-free hifi dataset will have an edge s-mer about every 250-300bp
		// so a human genome will have about 10000000-12000000 set bits
		// so about ~0.3% filled -> ~0.3% of kmers have a hash collision and are picked even though they are not edge k-mers
		// note that this also limits the effective window size to 250-300 regardless of what parameter is given
		endSmers.resize(4294967296, false);
		collectEndSmers();
	}
}

ReadpartIterator::~ReadpartIterator()
{
	if (cacheFileName.size() > 0)
	{
		remove(cacheFileName.c_str());
	}
}

void ReadpartIterator::collectEndSmers()
{
	size_t smerSize = kmerSize - windowSize + 1;
	size_t addedEndSmers = 0;
	std::mutex vectorMutex;
	iterateParts([this, smerSize, &addedEndSmers, &vectorMutex](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq) {
		if (seq.size() < smerSize) return;
		FastHasher startKmer { smerSize };
		FastHasher endKmer { smerSize };
		for (size_t i = 0; i < smerSize; i++)
		{
			startKmer.addChar(seq[i]);
			endKmer.addChar(complement(seq[seq.size()-1-i]));
		}
		{
			std::lock_guard<std::mutex> guard { vectorMutex };
			if (!endSmers[startKmer.hash() % endSmers.size()]) addedEndSmers += 1;
			if (!endSmers[endKmer.hash() % endSmers.size()]) addedEndSmers += 1;
			endSmers[startKmer.hash() % endSmers.size()] = true;
			endSmers[endKmer.hash() % endSmers.size()] = true;
		}
	});
	std::cerr << addedEndSmers << " end k-mers" << std::endl;
}
