#include <queue>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <string_view>
#include <phmap.h>
#include <thread>
#include "fastqloader.h"
#include "CommonUtils.h"
#include "MBGCommon.h"
#include "MBG.h"
#include "VectorWithDirection.h"
#include "TwobitString.h"
#include "FastHasher.h"
#include "SparseEdgeContainer.h"
#include "HashList.h"
#include "UnitigGraph.h"
#include "BluntGraph.h"
#include "ReadHelper.h"
#include "HPCConsensus.h"

struct AssemblyStats
{
public:
	size_t nodes;
	size_t edges;
	size_t size;
	size_t N50;
	size_t approxKmers;
};

std::vector<std::pair<std::string, std::vector<uint16_t>>> runLengthEncode(const std::string& original)
{
	assert(original.size() > 0);
	std::string resultStr;
	std::vector<uint16_t> lens;
	std::vector<std::pair<std::string, std::vector<uint16_t>>> result;
	resultStr.reserve(original.size());
	lens.reserve(original.size());
	lens.push_back(1);
	switch(original[0])
	{
		case 'a':
		case 'A':
			resultStr.push_back(1);
			break;
		case 'c':
		case 'C':
			resultStr.push_back(2);
			break;
		case 'g':
		case 'G':
			resultStr.push_back(3);
			break;
		case 't':
		case 'T':
			resultStr.push_back(4);
			break;
		default:
			break;
	}
	for (size_t i = 1; i < original.size(); i++)
	{
		switch(original[i])
		{
			case 'a':
			case 'A':
			case 'c':
			case 'C':
			case 'g':
			case 'G':
			case 't':
			case 'T':
			{
				if (original[i] == original[i-1])
				{
					lens.back() += 1;
					continue;
				}
				lens.push_back(1);
				switch(original[i])
				{
					case 'a':
					case 'A':
						resultStr.push_back(1);
						break;
					case 'c':
					case 'C':
						resultStr.push_back(2);
						break;
					case 'g':
					case 'G':
						resultStr.push_back(3);
						break;
					case 't':
					case 'T':
						resultStr.push_back(4);
						break;
					default:
						assert(false);
						break;
				}
				break;
			}
			default:
			{
				if (resultStr.size() == 0) continue;
				result.emplace_back(std::move(resultStr), std::move(lens));
				resultStr.clear();
				lens.clear();
				break;
			}
		}
	}
	if (resultStr.size() > 0)
	{
		result.emplace_back(std::move(resultStr), std::move(lens));
	}
	return result;
}

std::vector<std::pair<std::string, std::vector<uint16_t>>> noRunLengthEncode(const std::string& original)
{
	assert(original.size() > 0);
	std::vector<std::pair<std::string, std::vector<uint16_t>>> result;
	std::string resultStr;
	resultStr.reserve(original.size());
	for (size_t i = 0; i < original.size(); i++)
	{
		switch(original[i])
		{
			case 'a':
			case 'A':
				resultStr.push_back(1);
				break;
			case 'c':
			case 'C':
				resultStr.push_back(2);
				break;
			case 'g':
			case 'G':
				resultStr.push_back(3);
				break;
			case 't':
			case 'T':
				resultStr.push_back(4);
				break;
			default:
				if (resultStr.size() == 0) continue;
				result.emplace_back(std::move(resultStr), std::vector<uint16_t>{});
				resultStr.clear();
		}
	}
	if (resultStr.size() > 0) result.emplace_back(std::move(resultStr), std::vector<uint16_t>{});
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].second.resize(result[i].first.size(), 1);
	}
	return result;
}

char complement(char c)
{
	static std::vector<char> comp { 0, 4, 3, 2, 1 };
	return comp[c];
}

void collectEndSmers(std::vector<bool>& endSmer, const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const ReadpartIterator& partIterator)
{
	size_t smerSize = kmerSize - windowSize + 1;
	size_t addedEndSmers = 0;
	std::cerr << "Collecting end k-mers" << std::endl;
	iterateReadsMultithreaded(files, 1, [&endSmer, smerSize, &addedEndSmers, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iterateParts(read, [&endSmer, smerSize, &addedEndSmers](const std::string& seq, const std::vector<uint16_t>& lens) {
			FastHasher startKmer { smerSize };
			FastHasher endKmer { smerSize };
			for (size_t i = 0; i < smerSize; i++)
			{
				startKmer.addChar(seq[i]);
				endKmer.addChar(complement(seq[seq.size()-1-i]));
			}
			if (!endSmer[startKmer.hash() % endSmer.size()]) addedEndSmers += 1;
			if (!endSmer[endKmer.hash() % endSmer.size()]) addedEndSmers += 1;
			endSmer[startKmer.hash() % endSmer.size()] = true;
			endSmer[endKmer.hash() % endSmer.size()] = true;
		});
	});
	std::cerr << addedEndSmers << " end k-mers" << std::endl;
}

void loadReadsAsHashesMultithread(HashList& result, const std::vector<std::string>& files, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	std::atomic<size_t> totalNodes = 0;
	std::cerr << "Collecting selected k-mers" << std::endl;
	iterateReadsMultithreaded(files, numThreads, [&result, &totalNodes, kmerSize, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &totalNodes, kmerSize](const std::string& seq, const std::vector<uint16_t>& lens, uint64_t minHash, const std::vector<size_t>& positions)
		{
			std::string revSeq = revCompRLE(seq);
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			HashType lastHash = 0;
			for (auto pos : positions)
			{
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				std::string_view minimizerSequence { seq.data() + pos, kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
				std::pair<size_t, bool> current;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				std::tie(current, lastHash) = result.addNode(minimizerSequence, revMinimizerSequence, lens, pos, pos + kmerSize, lastHash, overlap, minHash);
				assert(pos - lastMinimizerPosition < kmerSize);
				if (last.first != std::numeric_limits<size_t>::max())
				{
					assert(lastMinimizerPosition + kmerSize >= pos);
					result.addSequenceOverlap(last, current, overlap);
					auto pair = canon(last, current);
					result.addEdgeCoverage(pair.first, pair.second);
				}
				lastMinimizerPosition = pos;
				last = current;
				totalNodes += 1;
			};
		});
	});
	std::cerr << totalNodes << " total selected k-mers in reads" << std::endl;
	std::cerr << result.size() << " distinct selected k-mers in reads" << std::endl;
}

std::vector<size_t> getRLEExpandedPositions(const std::string& seq, const std::vector<uint16_t>& lens)
{
	std::vector<size_t> result;
	result.resize(seq.size()+1);
	result[0] = 0;
	for (size_t i = 0; i < seq.size(); i++)
	{
		result[i+1] = result[i] + lens[i];
	}
	return result;
}

void outputSequencePathLine(std::ofstream& outPaths, const std::string& readName, const std::vector<size_t>& readExpandedPositions, const size_t currentSeqStart, const size_t currentSeqEnd, const std::vector<std::tuple<size_t, bool, size_t>>& kmerPath, const size_t currentPathStart, const size_t currentPathEnd, const std::vector<size_t>& unitigLength, const std::vector<size_t>& kmerUnitigStart, const std::vector<size_t>& kmerUnitigEnd, const UnitigGraph& unitigs, const HashList& hashlist)
{
	if (kmerPath.size() == 0) return;
	assert(readExpandedPositions.size() >= currentSeqEnd);
	assert(currentSeqEnd > currentSeqStart);
	assert(currentPathEnd > currentPathStart);
	assert(currentPathEnd - currentPathStart == currentSeqEnd - currentSeqStart);
	std::pair<size_t, bool> firstUnitig;
	std::pair<size_t, bool> lastUnitig;
	firstUnitig = std::make_pair(std::get<0>(kmerPath[0]), std::get<1>(kmerPath[0]));
	lastUnitig = firstUnitig;
	size_t endAddition = 0;
	std::string pathStr;
	if (firstUnitig.second)
	{
		pathStr += ">";
	}
	else
	{
		pathStr += "<";
	}
	// +1 because indices start at 0 but node ids at 1
	pathStr += std::to_string(firstUnitig.first+1);
	for (size_t i = 1; i < kmerPath.size(); i++)
	{
		if (std::get<2>(kmerPath[i]) == 0)
		{
			lastUnitig = std::make_pair(std::get<0>(kmerPath[i]), std::get<1>(kmerPath[i]));
			if (lastUnitig.second)
			{
				pathStr += ">";
			}
			else
			{
				pathStr += "<";
			}
			// +1 because indices start at 0 but node ids at 1
			pathStr += std::to_string(lastUnitig.first+1);
			std::pair<size_t, bool> from, to;
			if (std::get<1>(kmerPath[i-1]))
			{
				from = unitigs.unitigs[std::get<0>(kmerPath[i-1])].back();
			}
			else
			{
				from = reverse(unitigs.unitigs[std::get<0>(kmerPath[i-1])][0]);
			}
			if (std::get<1>(kmerPath[i]))
			{
				to = unitigs.unitigs[std::get<0>(kmerPath[i])][0];
			}
			else
			{
				to = reverse(unitigs.unitigs[std::get<0>(kmerPath[i])].back());
			}
			size_t overlap = hashlist.getOverlap(from, to);
			endAddition += unitigLength[std::get<0>(kmerPath[i-1])] - overlap;
		}
	}
	std::pair<size_t, bool> startKmer;
	if (std::get<1>(kmerPath[0]))
	{
		startKmer = unitigs.unitigs[std::get<0>(kmerPath[0])][std::get<2>(kmerPath[0])];
	}
	else
	{
		startKmer = reverse(unitigs.unitigs[std::get<0>(kmerPath[0])][unitigs.unitigs[std::get<0>(kmerPath[0])].size() - std::get<2>(kmerPath[0]) - 1]);
	}
	std::pair<size_t, bool> endKmer;
	if (std::get<1>(kmerPath.back()))
	{
		endKmer = unitigs.unitigs[std::get<0>(kmerPath.back())][std::get<2>(kmerPath.back())];
	}
	else
	{
		endKmer = reverse(unitigs.unitigs[std::get<0>(kmerPath.back())][unitigs.unitigs[std::get<0>(kmerPath.back())].size() - std::get<2>(kmerPath.back()) - 1]);
	}
	size_t RLExpandedPathStart = 0;
	size_t RLExpandedPathEnd = 0;
	if (startKmer.second)
	{
		RLExpandedPathStart = kmerUnitigStart[startKmer.first];
	}
	else
	{
		RLExpandedPathStart = unitigLength[firstUnitig.first] - 1 - kmerUnitigEnd[startKmer.first];
	}
	if (endKmer.second)
	{
		RLExpandedPathEnd = kmerUnitigEnd[endKmer.first] + 1;
	}
	else
	{
		RLExpandedPathEnd = unitigLength[lastUnitig.first] - kmerUnitigStart[endKmer.first];
	}
	RLExpandedPathEnd += endAddition;
	size_t RLExpandedSeqSize = readExpandedPositions.back();
	size_t RLExpandedSeqStart = readExpandedPositions[currentSeqStart];
	size_t RLExpandedSeqEnd = readExpandedPositions[currentSeqEnd];
	size_t RLExpandedPathSize = endAddition + unitigLength[lastUnitig.first];
	assert(RLExpandedSeqStart < RLExpandedSeqEnd);
	assert(RLExpandedSeqEnd <= RLExpandedSeqSize);
	assert(RLExpandedPathStart < RLExpandedPathEnd);
	assert(RLExpandedPathEnd <= RLExpandedPathSize);
	outPaths << readName << "\t" << RLExpandedSeqSize << "\t" << RLExpandedSeqStart << "\t" << RLExpandedSeqEnd << "\t" << "+" << "\t" << pathStr << "\t" << RLExpandedPathSize << "\t" << RLExpandedPathStart << "\t" << RLExpandedPathEnd << "\t" << (RLExpandedPathEnd - RLExpandedPathStart) << "\t" << (RLExpandedPathEnd - RLExpandedPathStart) << "\t" << "60" << std::endl;
}

// std::unordered_set<uint64_t> collectApproxHashes(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs)
// {
// 	std::unordered_set<uint64_t> approxHashes;
// 	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
// 	{
// 		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
// 		{
// 			FastHasher hasher { kmerSize };
// 			auto seq = hashlist.getHashSequenceRLE(unitigs.unitigs[i][j].first);
// 			for (size_t k = 0; k < kmerSize; k++)
// 			{
// 				hasher.addChar(seq[k]);
// 			}
// 			approxHashes.insert(hasher.hash());
// 		}
// 	}
// 	return approxHashes;
// }

// std::unordered_set<HashType> collectExactHashes(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs)
// {
// 	std::unordered_set<HashType> exactHashes;
// 	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
// 	{
// 		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
// 		{
// 			auto seq = hashlist.getHashSequenceRLE(unitigs.unitigs[i][j].first);
// 			exactHashes.insert(hash(seq.toString()));
// 			auto revseq = revCompRLE(seq);
// 			exactHashes.insert(hash(revseq.toString()));
// 		}
// 	}
// 	return exactHashes;
// }

// template <typename F>
// void findCollectedKmers(const std::string& seq, const size_t kmerSize, const std::unordered_set<uint64_t>& approxHashes, const std::unordered_set<HashType>& exactHashes, F callback)
// {
// 	FastHasher hasher { kmerSize };
// 	for (size_t i = 0; i < kmerSize; i++)
// 	{
// 		hasher.addChar(seq[i]);
// 	}
// 	if (approxHashes.count(hasher.hash()) == 1)
// 	{
// 		HashType exactHash = hash(seq.substr(0, kmerSize));
// 		if (exactHashes.count(exactHash) == 1)
// 		{
// 			callback(0, exactHash);
// 		}
// 	}
// 	for (size_t i = kmerSize; i < seq.size(); i++)
// 	{
// 		hasher.addChar(seq[i]);
// 		hasher.removeChar(seq[i-kmerSize]);
// 		if (approxHashes.count(hasher.hash()) == 0) continue;
// 		HashType exactHash = hash(seq.substr(i-kmerSize+1, kmerSize));
// 		if (exactHashes.count(exactHash) == 0) continue;
// 		callback(i - kmerSize + 1, exactHash);
// 	}
// }

// void writePaths(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::string>& inputReads, const size_t kmerSize, const ReadpartIterator& partIterator, const std::string& outputSequencePaths, const size_t numThreads)
// {
// 	std::unordered_set<uint64_t> approxHashes = collectApproxHashes(hashlist, kmerSize, unitigs);
// 	std::unordered_set<HashType> exactHashes = collectExactHashes(hashlist, kmerSize, unitigs);
// 	std::vector<std::tuple<size_t, bool, size_t>> kmerUnitigPosition;
// 	kmerUnitigPosition.resize(hashlist.size(), std::make_tuple(std::numeric_limits<size_t>::max(), true, std::numeric_limits<size_t>::max()));
// 	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
// 	{
// 		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
// 		{
// 			assert(std::get<0>(kmerUnitigPosition[unitigs.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
// 			std::get<0>(kmerUnitigPosition[unitigs.unitigs[i][j].first]) = i;
// 			std::get<1>(kmerUnitigPosition[unitigs.unitigs[i][j].first]) = unitigs.unitigs[i][j].second;
// 			if (unitigs.unitigs[i][j].second)
// 			{
// 				std::get<2>(kmerUnitigPosition[unitigs.unitigs[i][j].first]) = j;
// 			}
// 			else
// 			{
// 				std::get<2>(kmerUnitigPosition[unitigs.unitigs[i][j].first]) = unitigs.unitigs[i].size() - 1 - j;
// 			}
// 		}
// 	}
// 	std::vector<size_t> unitigLength;
// 	std::vector<size_t> kmerUnitigStart;
// 	std::vector<size_t> kmerUnitigEnd;
// 	unitigLength.resize(unitigs.unitigs.size(), std::numeric_limits<size_t>::max());
// 	kmerUnitigStart.resize(hashlist.size(), std::numeric_limits<size_t>::max());
// 	kmerUnitigEnd.resize(hashlist.size(), std::numeric_limits<size_t>::max());
// 	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
// 	{
// 		std::string seq;
// 		std::vector<uint16_t> lens;
// 		std::tie(seq, lens) = unitigs.getSequenceAndLength(i, hashlist);
// 		std::vector<size_t> RLEposes = getRLEExpandedPositions(seq, lens);
// 		unitigLength[i] = RLEposes.back();
// 		size_t seqPos = 0;
// 		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
// 		{
// 			if (j > 0)
// 			{
// 				size_t overlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
// 				assert(overlap > 0);
// 				assert(overlap < kmerSize);
// 				seqPos += kmerSize - overlap;
// 			}
// 			std::pair<size_t, bool> kmerHere = unitigs.unitigs[i][j];
// 			assert(kmerUnitigStart[kmerHere.first] == std::numeric_limits<size_t>::max());
// 			assert(kmerUnitigEnd[kmerHere.first] == std::numeric_limits<size_t>::max());
// 			assert(seqPos + kmerSize < RLEposes.size());
// 			if (kmerHere.second)
// 			{
// 				kmerUnitigStart[kmerHere.first] = RLEposes[seqPos];
// 				kmerUnitigEnd[kmerHere.first] = RLEposes[seqPos + kmerSize] - 1;
// 				assert(kmerUnitigStart[kmerHere.first] >= 0);
// 				assert(kmerUnitigStart[kmerHere.first] < kmerUnitigEnd[kmerHere.first]);
// 				assert(kmerUnitigEnd[kmerHere.first] < unitigLength[i]);
// 			}
// 			else
// 			{
// 				kmerUnitigEnd[kmerHere.first] = unitigLength[i] - RLEposes[seqPos] - 1;
// 				kmerUnitigStart[kmerHere.first] = unitigLength[i] - RLEposes[seqPos + kmerSize];
// 				assert(kmerUnitigStart[kmerHere.first] >= 0);
// 				assert(kmerUnitigStart[kmerHere.first] < kmerUnitigEnd[kmerHere.first]);
// 				assert(kmerUnitigEnd[kmerHere.first] < unitigLength[i]);
// 			}
// 		}
// 	}
// 	std::ofstream outPaths { outputSequencePaths };
// 	std::mutex pathWriteMutex;
// 	iterateReadsMultithreaded(inputReads, numThreads, [&outPaths, kmerSize, &approxHashes, &exactHashes, &kmerUnitigPosition, &unitigLength, &kmerUnitigStart, &kmerUnitigEnd, &hashlist, &unitigs, &pathWriteMutex, &partIterator](size_t thread, FastQ& read)
// 	{
// 		partIterator.iterateParts(read, [&read, &outPaths, kmerSize, &approxHashes, &exactHashes, &kmerUnitigPosition, &unitigLength, &kmerUnitigStart, &kmerUnitigEnd, &hashlist, &unitigs, &pathWriteMutex](const std::string& seq, const std::vector<uint16_t>& lens) {
// 			auto readExpandedPositions = getRLEExpandedPositions(seq, lens);
// 			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
// 			std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
// 			std::vector<std::tuple<size_t, bool, size_t>> kmerPath;
// 			size_t currentSeqStart = std::numeric_limits<size_t>::max();
// 			size_t currentPathStart = std::numeric_limits<size_t>::max();
// 			size_t currentSeqEnd = std::numeric_limits<size_t>::max();
// 			size_t currentPathEnd = std::numeric_limits<size_t>::max();
// 			std::vector<std::pair<size_t, std::pair<size_t, bool>>> matches;
// 			findCollectedKmers(seq, kmerSize, approxHashes, exactHashes, [&matches, &hashlist](size_t pos, HashType hash)
// 			{
// 				matches.emplace_back(pos, hashlist.hashToNode.at(hash));
// 			});
// 			for (size_t j = matches.size()-1; j < matches.size(); j--)
// 			{
// 				for (size_t k = j+1; k < matches.size(); k++)
// 				{
// 					size_t overlap = kmerSize - (matches[k].first - matches[j].first);
// 					auto canonEdge = canon(matches[j].second, matches[k].second);
// 					if (hashlist.sequenceOverlap[canonEdge.first].count(canonEdge.second) == 0) continue;
// 					if (overlap != hashlist.getOverlap(canonEdge.first, canonEdge.second)) continue;
// 					if (k == j+1) break;
// 					matches.erase(matches.begin()+j+1, matches.begin()+k);
// 					break;
// 				}
// 			}
// 			for (auto pair : matches)
// 			{
// 				size_t pos = pair.first;
// 				std::pair<size_t, bool> current = pair.second;
// 				size_t overlap = lastMinimizerPosition + kmerSize - pos;
// 				lastMinimizerPosition = pos;
// 				assert(current.first != std::numeric_limits<size_t>::max());
// 				assert(current.first < kmerUnitigPosition.size());
// 				std::tuple<size_t, bool, size_t> unitigPosition = kmerUnitigPosition[current.first];
// 				if (std::get<0>(unitigPosition) == std::numeric_limits<size_t>::max())
// 				{
// 					{
// 						std::lock_guard<std::mutex> guard { pathWriteMutex };
// 						outputSequencePathLine(outPaths, read.seq_id, readExpandedPositions, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
// 					}
// 					kmerPath.clear();
// 					currentSeqStart = std::numeric_limits<size_t>::max();
// 					currentSeqEnd = std::numeric_limits<size_t>::max();
// 					currentPathStart = std::numeric_limits<size_t>::max();
// 					currentPathEnd = std::numeric_limits<size_t>::max();
// 					lastKmer.first = std::numeric_limits<size_t>::max();
// 					continue;
// 				}
// 				if (!current.second)
// 				{
// 					std::get<1>(unitigPosition) = !std::get<1>(unitigPosition);
// 					if (std::get<0>(unitigPosition) != std::numeric_limits<size_t>::max())
// 					{
// 						std::get<2>(unitigPosition) = unitigs.unitigs[std::get<0>(unitigPosition)].size() - 1 - std::get<2>(unitigPosition);
// 					}
// 				}
// 				if (kmerPath.size() == 0)
// 				{
// 					currentSeqStart = pos;
// 					currentSeqEnd = pos + kmerSize;
// 					currentPathStart = kmerUnitigStart[current.first];
// 					if (!current.second)
// 					{
// 						currentPathStart = unitigLength[std::get<0>(unitigPosition)] - currentPathStart;
// 					}
// 					currentPathEnd = currentPathStart + kmerSize;
// 					kmerPath.emplace_back(unitigPosition);
// 					lastKmer = current;
// 					continue;
// 				}
// 				assert(lastKmer.first != std::numeric_limits<size_t>::max());
// 				bool hasMaybeValidEdge = true;
// 				auto canonEdge = canon(lastKmer, current);
// 				if (hashlist.sequenceOverlap[canonEdge.first].count(canonEdge.second) == 0)
// 				{
// 					hasMaybeValidEdge = false;
// 				}
// 				if (hasMaybeValidEdge && overlap != hashlist.getOverlap(lastKmer, current))
// 				{
// 					hasMaybeValidEdge = false;
// 				}
// 				if (hasMaybeValidEdge)
// 				{
// 					if (std::get<0>(unitigPosition) == std::get<0>(kmerPath.back()) && std::get<1>(unitigPosition) == std::get<1>(kmerPath.back()) && std::get<2>(unitigPosition) == std::get<2>(kmerPath.back()) + 1)
// 					{
// 						assert(overlap < kmerSize);
// 						currentSeqEnd += kmerSize - overlap;
// 						currentPathEnd += kmerSize - overlap;
// 						kmerPath.emplace_back(unitigPosition);
// 						lastKmer = current;
// 						continue;
// 					}
// 					if (std::get<2>(unitigPosition) == 0 && std::get<2>(kmerPath.back()) == unitigs.unitigs[std::get<0>(kmerPath.back())].size()-1)
// 					{
// 						std::pair<size_t, bool> fromUnitig = std::make_pair(std::get<0>(kmerPath.back()), std::get<1>(kmerPath.back()));
// 						std::pair<size_t, bool> toUnitig = std::make_pair(std::get<0>(unitigPosition), std::get<1>(unitigPosition));
// 						auto canonEdge = canon(fromUnitig, toUnitig);
// 						if (unitigs.edges[canonEdge.first].count(canonEdge.second) == 1)
// 						{
// 							// edge is valid
// 							assert(overlap < kmerSize);
// 							currentSeqEnd += kmerSize - overlap;
// 							currentPathEnd += kmerSize - overlap;
// 							kmerPath.emplace_back(unitigPosition);
// 							lastKmer = current;
// 							continue;
// 						}
// 					}
// 				}
// 				// no valid edge
// 				{
// 					std::lock_guard<std::mutex> guard { pathWriteMutex };
// 					outputSequencePathLine(outPaths, read.seq_id, readExpandedPositions, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
// 				}
// 				kmerPath.clear();
// 				currentSeqStart = pos;
// 				currentSeqEnd = pos + kmerSize;
// 				currentPathStart = kmerUnitigStart[current.first];
// 				if (!current.second)
// 				{
// 					currentPathStart = unitigLength[std::get<0>(unitigPosition)] - currentPathStart;
// 				}
// 				currentPathEnd = currentPathStart + kmerSize;
// 				kmerPath.emplace_back(unitigPosition);
// 				lastKmer = current;
// 			}
// 			{
// 				std::lock_guard<std::mutex> guard { pathWriteMutex };
// 				outputSequencePathLine(outPaths, read.seq_id, readExpandedPositions, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
// 			}
// 		});
// 	});
// }

void startUnitig(UnitigGraph& result, const UnitigGraph& old, std::pair<size_t, bool> start, const VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>>& edges, std::vector<std::pair<size_t, bool>>& belongsToUnitig, VectorWithDirection<bool>& unitigStart, VectorWithDirection<bool>& unitigEnd)
{
	size_t currentUnitig = result.unitigs.size();
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	std::pair<size_t, bool> pos = start;
	unitigStart[start] = true;
	unitigEnd[reverse(start)] = true;
	assert(belongsToUnitig.at(pos.first).first == std::numeric_limits<size_t>::max());
	belongsToUnitig[pos.first] = std::make_pair(currentUnitig, pos.second);
	if (pos.second)
	{
		result.unitigs.back().insert(result.unitigs.back().end(), old.unitigs[pos.first].begin(), old.unitigs[pos.first].end());
		result.unitigCoverage.back().insert(result.unitigCoverage.back().end(), old.unitigCoverage[pos.first].begin(), old.unitigCoverage[pos.first].end());
	}
	else
	{
		for (size_t i = old.unitigs[pos.first].size()-1; i < old.unitigs[pos.first].size(); i--)
		{
			result.unitigs.back().emplace_back(old.unitigs[pos.first][i].first, !old.unitigs[pos.first][i].second);
			result.unitigCoverage.back().emplace_back(old.unitigCoverage[pos.first][i]);
		}
	}
	while (true)
	{
		if (edges.at(pos).size() != 1) break;
		auto newPos = *edges.at(pos).begin();
		auto revPos = std::make_pair(newPos.first, !newPos.second);
		if (edges.at(revPos).size() != 1) break;
		if (newPos == start)
		{
			result.edges[std::make_pair(currentUnitig, true)].emplace(currentUnitig, true);
			result.edgeCoverage(currentUnitig, true, currentUnitig, true) = old.edgeCoverage(pos.first, pos.second, newPos.first, newPos.second);
			break;
		}
		if (belongsToUnitig.at(newPos.first).first != std::numeric_limits<size_t>::max())
		{
			assert(newPos.first == pos.first);
			assert(newPos.second != pos.second);
			assert(belongsToUnitig.at(newPos.first).first == currentUnitig);
			assert(belongsToUnitig.at(newPos.first).second != newPos.second);
			result.edges[std::make_pair(currentUnitig, belongsToUnitig.at(pos.first).second)].emplace(currentUnitig, !belongsToUnitig.at(pos.first).second);
			result.edgeCoverage(currentUnitig, belongsToUnitig.at(pos.first).second, currentUnitig, !belongsToUnitig.at(pos.first).second) = old.edgeCoverage(pos.first, pos.second, newPos.first, newPos.second);
			break;
		}
		pos = newPos;
		assert(belongsToUnitig.at(pos.first).first == std::numeric_limits<size_t>::max());
		belongsToUnitig[pos.first] = std::make_pair(currentUnitig, pos.second);
		if (pos.second)
		{
			result.unitigs.back().insert(result.unitigs.back().end(), old.unitigs[pos.first].begin(), old.unitigs[pos.first].end());
			result.unitigCoverage.back().insert(result.unitigCoverage.back().end(), old.unitigCoverage[pos.first].begin(), old.unitigCoverage[pos.first].end());
		}
		else
		{
			for (size_t i = old.unitigs[pos.first].size()-1; i < old.unitigs[pos.first].size(); i--)
			{
				result.unitigs.back().emplace_back(old.unitigs[pos.first][i].first, !old.unitigs[pos.first][i].second);
				result.unitigCoverage.back().emplace_back(old.unitigCoverage[pos.first][i]);
			}
		}
	}
	unitigEnd[pos] = true;
	unitigStart[reverse(pos)] = true;
}

void startUnitig(UnitigGraph& result, std::pair<size_t, bool> start, const SparseEdgeContainer& edges, std::vector<bool>& belongsToUnitig, const HashList& hashlist, size_t minCoverage)
{
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	std::pair<size_t, bool> pos = start;
	assert(!belongsToUnitig[pos.first]);
	belongsToUnitig[pos.first] = true;
	result.unitigs.back().emplace_back(pos);
	result.unitigCoverage.back().emplace_back(hashlist.coverage[pos.first]);
	while (true)
	{
		if (edges[pos].size() != 1) break;
		auto newPos = edges[pos][0];
		auto revPos = std::make_pair(newPos.first, !newPos.second);
		if (edges[revPos].size() != 1) break;
		if (newPos == start)
		{
			break;
		}
		if (belongsToUnitig[newPos.first])
		{
			assert(newPos.first == pos.first);
			assert(newPos.second != pos.second);
			break;
		}
		pos = newPos;
		assert(!belongsToUnitig[pos.first]);
		belongsToUnitig[pos.first] = true;
		result.unitigs.back().emplace_back(pos);
		result.unitigCoverage.back().emplace_back(hashlist.coverage[pos.first]);
	}
}

SparseEdgeContainer getCoveredEdges(const HashList& hashlist, size_t minCoverage)
{
	SparseEdgeContainer result { hashlist.coverage.size() };
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : hashlist.edgeCoverage.at(fw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(fw, edge.first);
			result.addEdge(reverse(edge.first), reverse(fw));
		}
		for (auto edge : hashlist.edgeCoverage.at(bw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(bw, edge.first);
			result.addEdge(reverse(edge.first), reverse(bw));
		}
	}
	return result;
}

UnitigGraph getUnitigGraph(const HashList& hashlist, size_t minCoverage)
{
	UnitigGraph result;
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(hashlist.coverage.size(), false);
	std::unordered_map<std::pair<size_t, bool>, std::pair<size_t, bool>> unitigTip;
	auto edges = getCoveredEdges(hashlist, minCoverage);
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage[i] < minCoverage) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		auto fwEdges = edges[fw];
		auto bwEdges = edges[bw];
		if (bwEdges.size() != 1)
		{
			if (!belongsToUnitig[i])
			{
				startUnitig(result, fw, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
			for (auto edge : bwEdges)
			{
				if (belongsToUnitig[edge.first]) continue;
				assert(hashlist.coverage[edge.first] >= minCoverage);
				startUnitig(result, edge, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
		}
		if (fwEdges.size() != 1)
		{
			if (!belongsToUnitig[i])
			{
				startUnitig(result, bw, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
			for (auto edge : fwEdges)
			{
				if (belongsToUnitig[edge.first]) continue;
				assert(hashlist.coverage[edge.first] >= minCoverage);
				startUnitig(result, edge, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
		}
	}
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (belongsToUnitig[i]) continue;
		if (hashlist.coverage[i] < minCoverage) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		auto fwEdges = edges[fw];
		auto bwEdges = edges[bw];
		assert(fwEdges.size() == 1);
		assert(bwEdges.size() == 1);
		startUnitig(result, fw, edges, belongsToUnitig, hashlist, minCoverage);
		assert(result.unitigs.size() > 0);
		assert(result.unitigs.back().back() == reverse(bwEdges[0]));
		unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
		unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
	}
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage[i] < minCoverage) continue;
		assert(belongsToUnitig[i]);
	}
	for (auto tip : unitigTip)
	{
		auto fromNode = tip.first;
		auto fromUnitig = tip.second;
		for (auto edge : edges[fromNode])
		{
			auto toNodeFw = edge;
			auto toNodeRev = reverse(edge);
			assert(unitigTip.count(toNodeRev) == 1);
			auto toUnitig = reverse(unitigTip.at(toNodeRev));
			assert(hashlist.coverage[fromNode.first] >= minCoverage);
			assert(hashlist.coverage[toNodeFw.first] >= minCoverage);
			result.edges[fromUnitig].emplace(toUnitig);
			result.edges[reverse(toUnitig)].emplace(reverse(fromUnitig));
			result.edgeCoverage(fromUnitig, toUnitig) = hashlist.getEdgeCoverage(fromNode, toNodeFw);
		}
	}
	return result;
}

UnitigGraph getUnitigs(const UnitigGraph& oldgraph)
{
	UnitigGraph result;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	edges.resize(oldgraph.unitigs.size());
	for (size_t i = 0; i < oldgraph.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : oldgraph.edges[fw])
		{
			edges[fw].emplace(to);
			edges[reverse(to)].emplace(reverse(fw));
		}
		for (auto to : oldgraph.edges[bw])
		{
			edges[bw].emplace(to);
			edges[reverse(to)].emplace(reverse(bw));
		}
	}
	std::vector<std::pair<size_t, bool>> belongsToUnitig;
	belongsToUnitig.resize(oldgraph.unitigs.size(), std::make_pair(std::numeric_limits<size_t>::max(), true));
	VectorWithDirection<bool> unitigStart;
	VectorWithDirection<bool> unitigEnd;
	unitigStart.resize(oldgraph.unitigs.size(), false);
	unitigEnd.resize(oldgraph.unitigs.size(), false);
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		auto bw = std::make_pair(node, false);
		if (edges.at(fw).size() != 1)
		{
			for (auto start : edges.at(fw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig, unitigStart, unitigEnd);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, bw, edges, belongsToUnitig, unitigStart, unitigEnd);
		}
		if (edges.at(bw).size() != 1)
		{
			for (auto start : edges.at(bw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig, unitigStart, unitigEnd);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig, unitigStart, unitigEnd);
		}
	}
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig, unitigStart, unitigEnd);
	}
	for (size_t i = 0; i < oldgraph.edges.size(); i++)
	{
		auto fw = std::make_pair(i, true);
		if (unitigEnd[fw])
		{
			for (auto curr : oldgraph.edges.at(fw))
			{
				if (!unitigStart[curr]) continue;
				auto from = belongsToUnitig.at(fw.first);
				bool prevFw = fw.second;
				auto to = belongsToUnitig.at(curr.first);
				bool currFw = curr.second;
				from.second = !(from.second ^ prevFw);
				to.second = !(to.second ^ currFw);
				result.edges[from].emplace(to);
				result.edgeCoverage(from, to) = oldgraph.edgeCoverage(fw, curr);
			}
		}
		auto bw = std::make_pair(i, false);
		if (unitigEnd[bw])
		{
			for (auto curr : oldgraph.edges.at(bw))
			{
				if (!unitigStart[curr]) continue;
				auto from = belongsToUnitig.at(bw.first);
				bool prevFw = bw.second;
				auto to = belongsToUnitig.at(curr.first);
				bool currFw = curr.second;
				from.second = !(from.second ^ prevFw);
				to.second = !(to.second ^ currFw);
				result.edges[from].emplace(to);
				result.edgeCoverage(from, to) = oldgraph.edgeCoverage(bw, curr);
			}
		}
	}
	return result;
}

std::string getSequence(const std::string& rle)
{
	std::string result;
	for (size_t i = 0; i < rle.size(); i++)
	{
		result += "-ACGT"[(int)rle[i]];
	}
	return result;
}

std::string getSequence(const std::string& rle, const std::vector<uint16_t>& characterLength)
{
	std::string result;
	assert(rle.size() == characterLength.size());
	for (size_t i = 0; i < rle.size(); i++)
	{
		for (size_t j = 0; j < characterLength[i]; j++)
		{
			result += "-ACGT"[(int)rle[i]];
		}
	}
	return result;
}

size_t getOverlapFromRLE(const std::vector<std::pair<std::string, std::vector<uint16_t>>>& unitigSequences, std::pair<size_t, bool> fromUnitig, size_t rleOverlap)
{
	assert(fromUnitig.first < unitigSequences.size());
	assert(rleOverlap < unitigSequences[fromUnitig.first].first.size());
	size_t result = 0;
	for (size_t i = 0; i < rleOverlap; i++)
	{
		size_t index = i;
		if (fromUnitig.second) index = unitigSequences[fromUnitig.first].first.size() - i - 1;
		result += unitigSequences[fromUnitig.first].second[index];
	}
	return result;
}

size_t getN50(std::vector<size_t>& nodeSizes, size_t totalSize)
{
	std::sort(nodeSizes.begin(), nodeSizes.end());
	size_t sizeSum = 0;
	sizeSum += nodeSizes.back();
	while (sizeSum * 2 < totalSize)
	{
		nodeSizes.pop_back();
		sizeSum += nodeSizes.back();
	}
	return nodeSizes.back();
}

AssemblyStats writeGraph(const BluntGraph& graph, const std::string& filename)
{
	AssemblyStats stats;
	stats.size = 0;
	stats.nodes = graph.nodes.size();
	stats.edges = graph.edges.size();
	stats.N50 = 0;
	stats.approxKmers = 0;
	std::vector<size_t> nodeSizes;
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		file << "S\t" << (i+1) << "\t" << graph.nodes[i] << "\tll:f:" << graph.nodeAvgCoverage[i] << "\tFC:f:" << (graph.nodes[i].size() * graph.nodeAvgCoverage[i]) << std::endl;
		stats.size += graph.nodes[i].size();
		nodeSizes.push_back(graph.nodes[i].size());
	}
	for (auto edge : graph.edges)
	{
		file << "L\t" << (std::get<0>(edge)+1) << "\t" << (std::get<1>(edge) ? "+" : "-") << "\t" << (std::get<2>(edge)+1) << "\t" << (std::get<3>(edge) ? "+" : "-") << "\t0M\tec:i:" << std::get<4>(edge) << std::endl;
	}
	stats.N50 = getN50(nodeSizes, stats.size);
	return stats;
}

AssemblyStats writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist, const std::vector<std::pair<std::string, std::vector<uint16_t>>>& unitigSequences, const size_t kmerSize)
{
	AssemblyStats stats;
	stats.size = 0;
	stats.nodes = unitigs.unitigs.size();
	stats.edges = 0;
	stats.N50 = 0;
	stats.approxKmers = 0;
	std::vector<size_t> nodeSizes;
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		file << "S\t" << (i+1) << "\t";
		std::string realSequence = getHPCExpanded(unitigSequences[i].first, unitigSequences[i].second);
		file << realSequence;
		file << "\tll:f:" << unitigs.averageCoverage(i);
		file << "\tFC:f:" << (unitigs.averageCoverage(i) * realSequence.size());
		file << std::endl;
		stats.size += realSequence.size();
		nodeSizes.push_back(realSequence.size());
	}
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : unitigs.edges[fw])
		{
			if (canon(fw, to).first == fw) stats.edges += 1;
			std::pair<size_t, bool> last = unitigs.unitigs[i].back();
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t rleOverlap = hashlist.getOverlap(last, first);
			size_t overlap = getOverlapFromRLE(unitigSequences, fw, rleOverlap);
			file << "L\t" << (fw.first+1) << "\t" << (fw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		for (auto to : unitigs.edges[bw])
		{
			if (canon(bw, to).first == bw) stats.edges += 1;
			std::pair<size_t, bool> last = reverse(unitigs.unitigs[i][0]);
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t rleOverlap = hashlist.getOverlap(last, first);
			size_t overlap = getOverlapFromRLE(unitigSequences, bw, rleOverlap);
			file << "L\t" << (bw.first+1) << "\t" << (bw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(bw, to) << std::endl;
		}
	}
	stats.N50 = getN50(nodeSizes, stats.size);
	stats.approxKmers = stats.size - stats.nodes * kmerSize;
	return stats;
}

auto getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

AssemblyStats getSizeAndN50(const BluntGraph& graph)
{
	AssemblyStats result;
	result.nodes = graph.nodes.size();
	result.edges = graph.edges.size();
	std::vector<size_t> sizes;
	size_t totalSize = 0;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		sizes.push_back(graph.nodes[i].size());
		totalSize += graph.nodes[i].size();
	}
	result.size = totalSize;
	result.approxKmers = 0;
	std::sort(sizes.begin(), sizes.end());
	result.N50 = 0;
	size_t partialSum = 0;
	for (size_t i = sizes.size()-1; i < sizes.size(); i--)
	{
		partialSum += sizes[i];
		if (partialSum >= totalSize * 0.5)
		{
			result.N50 = sizes[i];
			break;
		}
	}
	return result;
}

std::pair<NodeType, size_t> find(std::unordered_map<NodeType, std::vector<std::pair<NodeType, size_t>>>& parent, std::pair<NodeType, size_t> key)
{
	while (true)
	{
		assert(parent.count(key.first) == 1);
		auto foundParent = parent[key.first][key.second];
		assert(parent.count(foundParent.first) == 1);
		if (parent[key.first][key.second] == parent[foundParent.first][foundParent.second]) return parent[key.first][key.second];
		parent[key.first][key.second] = parent[foundParent.first][foundParent.second];
	}
}

void merge(std::unordered_map<NodeType, std::vector<std::pair<NodeType, size_t>>>& parent, std::unordered_map<NodeType, std::vector<size_t>>& rank, std::pair<NodeType, size_t> left, std::pair<NodeType, size_t> right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent[left.first][left.second] == left);
	assert(parent[right.first][right.second] == right);
	assert(rank.count(left.first) == 1);
	assert(rank.count(right.first) == 1);
	if (rank[right.first][right.second] > rank[left.first][left.second]) std::swap(left, right);
	parent[right.first][right.second] = left;
	assert(rank[right.first][right.second] <= rank[left.first][left.second]);
	if (rank[right.first][right.second] == rank[left.first][left.second]) rank[left.first][left.second] += 1;
}

void merge(std::unordered_map<NodeType, std::vector<std::pair<NodeType, size_t>>>& parent, std::unordered_map<NodeType, std::vector<size_t>>& rank, NodeType leftNode, size_t leftOffset, NodeType rightNode, size_t rightOffset)
{
	merge(parent, rank, std::make_pair(leftNode, leftOffset), std::make_pair(rightNode, rightOffset));
}

// void forceEdgeConsistency(const UnitigGraph& unitigs, HashList& hashlist, const size_t kmerSize)
// {
// 	std::unordered_map<NodeType, std::vector<std::pair<NodeType, size_t>>> parent;
// 	std::unordered_map<NodeType, std::vector<size_t>> rank;
// 	for (const auto& unitig : unitigs.unitigs)
// 	{
// 		assert(unitig.size() > 0);
// 		assert(parent.count(unitig[0].first) == 0);
// 		assert(parent.count(unitig.back().first) == 0);
// 		assert(rank.count(unitig[0].first) == 0);
// 		assert(rank.count(unitig.back().first) == 0);
// 		for (size_t i = 0; i < kmerSize; i++)
// 		{
// 			parent[unitig[0].first].emplace_back(unitig[0].first, i);
// 			rank[unitig[0].first].emplace_back(0);
// 			if (unitig.size() > 1)
// 			{
// 				parent[unitig.back().first].emplace_back(unitig.back().first, i);
// 				rank[unitig.back().first].emplace_back(0);
// 			}
// 		}
// 		if (unitig.size() == 1) continue;
// 		size_t firstToLastOverlap = kmerSize;
// 		for (size_t i = 1; i < unitig.size(); i++)
// 		{
// 			size_t overlap = hashlist.getOverlap(unitig[i-1], unitig[i]);
// 			assert(overlap < kmerSize);
// 			size_t removeOverlap = kmerSize - overlap;
// 			if (removeOverlap >= firstToLastOverlap)
// 			{
// 				firstToLastOverlap = 0;
// 				break;
// 			}
// 			firstToLastOverlap -= removeOverlap;
// 			if (i == unitig.size()-1) break;
// 		}
// 		if (firstToLastOverlap > 0)
// 		{
// 			assert(firstToLastOverlap < kmerSize);
// 			for (size_t i = 0; i < firstToLastOverlap; i++)
// 			{
// 				size_t firstOffset = kmerSize - firstToLastOverlap + i;
// 				assert(firstOffset < kmerSize);
// 				if (!unitig[0].second) firstOffset = kmerSize - 1 - firstOffset;
// 				assert(firstOffset < kmerSize);
// 				size_t secondOffset = i;
// 				if (!unitig.back().second) secondOffset = kmerSize - 1 - secondOffset;
// 				assert(secondOffset < kmerSize);
// 				merge(parent, rank, unitig[0].first, firstOffset, unitig.back().first, secondOffset);
// 			}
// 		}
// 	}
// 	for (size_t unitig = 0; unitig < unitigs.edges.size(); unitig++)
// 	{
// 		std::pair<size_t, bool> fw { unitig, true };
// 		std::pair<size_t, bool> bw { unitig, false };
// 		for (auto target : unitigs.edges[fw])
// 		{
// 			std::pair<NodeType, bool> from = unitigs.unitigs[unitig].back();
// 			std::pair<NodeType, bool> to = unitigs.unitigs[target.first][0];
// 			if (!target.second) to = reverse(unitigs.unitigs[target.first].back());
// 			size_t overlap = hashlist.getOverlap(from, to);
// 			for (size_t i = 0; i < overlap; i++)
// 			{
// 				size_t firstOffset = kmerSize - overlap + i;
// 				assert(firstOffset < kmerSize);
// 				if (!from.second) firstOffset = kmerSize - 1 - firstOffset;
// 				assert(firstOffset < kmerSize);
// 				size_t secondOffset = i;
// 				if (!to.second) secondOffset = kmerSize - 1 - secondOffset;
// 				assert(secondOffset < kmerSize);
// 				merge(parent, rank, from.first, firstOffset, to.first, secondOffset);
// 			}
// 		}
// 		for (auto target : unitigs.edges[bw])
// 		{
// 			std::pair<NodeType, bool> from = reverse(unitigs.unitigs[unitig][0]);
// 			std::pair<NodeType, bool> to = unitigs.unitigs[target.first][0];
// 			if (!target.second) to = reverse(unitigs.unitigs[target.first].back());
// 			size_t overlap = hashlist.getOverlap(from, to);
// 			for (size_t i = 0; i < overlap; i++)
// 			{
// 				size_t firstOffset = kmerSize - overlap + i;
// 				assert(firstOffset < kmerSize);
// 				if (!from.second) firstOffset = kmerSize - 1 - firstOffset;
// 				assert(firstOffset < kmerSize);
// 				size_t secondOffset = i;
// 				if (!to.second) secondOffset = kmerSize - 1 - secondOffset;
// 				assert(secondOffset < kmerSize);
// 				merge(parent, rank, from.first, firstOffset, to.first, secondOffset);
// 			}
// 		}
// 	}
// 	std::unordered_set<NodeType> keys;
// 	for (auto pair : parent)
// 	{
// 		keys.emplace(pair.first);
// 	}
// 	for (auto node : keys)
// 	{
// 		for (size_t i = 0; i < kmerSize; i++)
// 		{
// 			auto found = find(parent, std::make_pair(node, i));
// 			size_t length = hashlist.getRunLength(found.first, found.second);
// 			assert(length > 0);
// 			hashlist.setRunLength(node, i, length);
// 		}
// 	}
// }

// void collectApproxNeighbors(std::unordered_map<uint64_t, unsigned char>& approxNeighbors, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
// {
// 	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
// 	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
// 	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
// 	if (!toKmer.second) toSequence = revCompRLE(toSequence);
// 	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
// 	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
// 	assert(edgeSequence.size() > kmerSize);
// 	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
// 	assert(edgeSequence.size() < 2 * kmerSize);
// 	if (edgeSequence.size() == kmerSize + 1) return;
// 	FastHasher fwkmerHasher { kmerSize };
// 	for (size_t i = 0; i < kmerSize; i++)
// 	{
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 	}
// 	for (size_t i = kmerSize; i < edgeSequence.size(); i++)
// 	{
// 		auto oldFwHash = fwkmerHasher.getFwHash();
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
// 		auto newBwHash = fwkmerHasher.getBwHash();
// 		if (approxNeighbors.count(oldFwHash) == 1 && approxNeighbors.at(oldFwHash) != edgeSequence[i]) approxNeighbors[oldFwHash] = std::numeric_limits<unsigned char>::max();
// 		if (approxNeighbors.count(oldFwHash) == 0) approxNeighbors[oldFwHash] = edgeSequence[i];
// 		if (approxNeighbors.count(newBwHash) == 1 && approxNeighbors.at(newBwHash) != revCompRLE(edgeSequence[i-kmerSize])) approxNeighbors[newBwHash] = std::numeric_limits<unsigned char>::max();
// 		if (approxNeighbors.count(newBwHash) == 0) approxNeighbors[newBwHash] = revCompRLE(edgeSequence[i-kmerSize]);
// 	}
// }

// void collectExactNeighbors(const std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::unordered_map<HashType, unsigned char>& exactNeighbors, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
// {
// 	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
// 	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
// 	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
// 	if (!toKmer.second) toSequence = revCompRLE(toSequence);
// 	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
// 	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
// 	assert(edgeSequence.size() > kmerSize);
// 	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
// 	assert(edgeSequence.size() < 2 * kmerSize);
// 	if (edgeSequence.size() == kmerSize + 1) return;
// 	FastHasher fwkmerHasher { kmerSize };
// 	for (size_t i = 0; i < kmerSize; i++)
// 	{
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 	}
// 	for (size_t i = kmerSize; i < edgeSequence.size(); i++)
// 	{
// 		auto oldFwHash = fwkmerHasher.getFwHash();
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
// 		auto newBwHash = fwkmerHasher.getBwHash();
// 		assert(approxNeighbors.count(oldFwHash) == 1);
// 		assert(approxNeighbors.count(newBwHash) == 1);
// 		// neither is a junction node if both not in approx
// 		if (approxNeighbors.at(oldFwHash) != std::numeric_limits<unsigned char>::max() && approxNeighbors.at(newBwHash) != std::numeric_limits<unsigned char>::max()) continue;
// 		HashType oldFwExactHash = hash(edgeSequence.substr(i-kmerSize, kmerSize));
// 		HashType newBwExactHash = hash(revCompRLE(edgeSequence.substr(i-kmerSize+1, kmerSize)));
// 		if (exactNeighbors.count(oldFwExactHash) == 1 && exactNeighbors.at(oldFwExactHash) != edgeSequence[i]) exactNeighbors[oldFwExactHash] = std::numeric_limits<unsigned char>::max();
// 		if (exactNeighbors.count(oldFwExactHash) == 0) exactNeighbors[oldFwExactHash] = edgeSequence[i];
// 		if (exactNeighbors.count(newBwExactHash) == 1 && exactNeighbors.at(newBwExactHash) != revCompRLE(edgeSequence[i-kmerSize])) exactNeighbors[newBwExactHash] = std::numeric_limits<unsigned char>::max();
// 		if (exactNeighbors.count(newBwExactHash) == 0) exactNeighbors[newBwExactHash] = revCompRLE(edgeSequence[i-kmerSize]);
// 	}
// }

// void collectJunctionHashes(const std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::unordered_set<uint64_t>& approxHashes, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
// {
// 	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
// 	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
// 	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
// 	if (!toKmer.second) toSequence = revCompRLE(toSequence);
// 	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
// 	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
// 	assert(edgeSequence.size() > kmerSize);
// 	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
// 	assert(edgeSequence.size() < 2 * kmerSize);
// 	if (edgeSequence.size() == kmerSize + 1) return;
// 	FastHasher fwkmerHasher { kmerSize };
// 	for (size_t i = 0; i < kmerSize; i++)
// 	{
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 	}
// 	for (size_t i = kmerSize; i < edgeSequence.size(); i++)
// 	{
// 		auto oldFwHash = fwkmerHasher.getFwHash();
// 		auto oldHash = fwkmerHasher.hash();
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
// 		auto newBwHash = fwkmerHasher.getBwHash();
// 		assert(approxNeighbors.count(oldFwHash) == 1);
// 		assert(approxNeighbors.count(newBwHash) == 1);
// 		// neither is a junction node if both not in approx
// 		if (approxNeighbors.at(oldFwHash) != std::numeric_limits<unsigned char>::max() && approxNeighbors.at(newBwHash) != std::numeric_limits<unsigned char>::max()) continue;
// 		auto newHash = fwkmerHasher.hash();
// 		approxHashes.insert(oldHash);
// 		approxHashes.insert(newHash);
// 	}
// }

// bool addJunctionsToHashes(const std::unordered_set<uint64_t>& approxHashes, std::unordered_map<NodeType, std::pair<size_t, bool>>& belongsToUnitig, std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>& addEdges, HashList& hashlist, UnitigGraph& unitigs, const size_t kmerSize, std::pair<size_t, bool> unitigPosition, std::pair<size_t, bool> toUnitigPosition, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
// {
// 	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
// 	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
// 	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
// 	if (!toKmer.second) toSequence = revCompRLE(toSequence);
// 	size_t bigOverlap = hashlist.getOverlap(fromKmer, toKmer); // actually the smallest overlap in the entire function
// 	assert(bigOverlap < kmerSize);
// 	std::string edgeSequence = fromSequence + toSequence.substr(bigOverlap);
// 	assert(edgeSequence.size() > kmerSize);
// 	assert(edgeSequence.size() < 2 * kmerSize);
// 	if (edgeSequence.size() == kmerSize + 1) return false;
// 	std::string revSeq = revCompRLE(edgeSequence);
// 	std::vector<uint16_t> lens;
// 	lens = hashlist.getHashCharacterLength(fromKmer.first);
// 	if (!fromKmer.second) std::reverse(lens.begin(), lens.end());
// 	std::vector<uint16_t> tolens;
// 	tolens = hashlist.getHashCharacterLength(toKmer.first);
// 	if (!toKmer.second) std::reverse(tolens.begin(), tolens.end());
// 	lens.insert(lens.end(), tolens.begin() + bigOverlap, tolens.end());
// 	assert(lens.size() == edgeSequence.size());
// 	FastHasher fwkmerHasher { kmerSize };
// 	for (size_t i = 0; i < kmerSize; i++)
// 	{
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 	}
// 	std::pair<NodeType, bool> lastKmer = fromKmer;
// 	std::pair<NodeType, bool> lastUnitig = unitigPosition;
// 	size_t lastPosition = 0;
// 	bool addedAny = false;
// 	size_t edgeCoverage = unitigs.edgeCoverage(unitigPosition, toUnitigPosition);
// 	assert(edgeCoverage > 0);
// 	for (size_t i = kmerSize; i < edgeSequence.size()-1; i++)
// 	{
// 		fwkmerHasher.addChar(edgeSequence[i]);
// 		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
// 		if (approxHashes.count(fwkmerHasher.hash()) == 0) continue;
// 		size_t kmerStartPos = i - kmerSize + 1;

// 		size_t overlap = kmerSize - kmerStartPos + lastPosition;
// 		assert(overlap < kmerSize);
// 		HashType lastHash = 0;
// 		uint64_t minHash = 0;

// 		std::string_view minimizerSequence { edgeSequence.data() + kmerStartPos, kmerSize };
// 		size_t revPos = edgeSequence.size() - (kmerStartPos + kmerSize);
// 		std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
// 		std::pair<size_t, bool> current;
// 		std::tie(current, lastHash) = hashlist.addNode(minimizerSequence, revMinimizerSequence, lens, kmerStartPos, kmerStartPos+kmerSize, lastHash, overlap, minHash);
// 		hashlist.addSequenceOverlap(lastKmer, current, overlap);
// 		std::pair<size_t, bool> newUnitigPosition { std::numeric_limits<size_t>::max(), true };
// 		if (belongsToUnitig.count(current.first) == 0)
// 		{
// 			unitigs.unitigs.emplace_back();
// 			unitigs.unitigs.back().emplace_back(current.first, true);
// 			unitigs.unitigCoverage.emplace_back();
// 			unitigs.unitigCoverage.back().emplace_back(0);
// 			unitigs.unitigCoverage.back()[0] += edgeCoverage;
// 			unitigs.edges.emplace_back();
// 			unitigs.edgeCov.emplace_back();
// 			belongsToUnitig[current.first] = std::make_pair(unitigs.unitigs.size()-1, true);
// 			newUnitigPosition = std::make_pair(unitigs.unitigs.size()-1, true);
// 			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
// 		}
// 		else
// 		{
// 			newUnitigPosition = belongsToUnitig.at(current.first);
// 			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
// 			assert(unitigs.unitigCoverage[newUnitigPosition.first].size() == 1);
// 			unitigs.unitigCoverage[newUnitigPosition.first][0] += edgeCoverage;
// 		}
// 		assert(newUnitigPosition.first != std::numeric_limits<size_t>::max());
// 		unitigs.edgeCoverage(lastUnitig, newUnitigPosition) += edgeCoverage;
// 		auto key = canon(lastUnitig, newUnitigPosition);
// 		addEdges.insert(key);
// 		lastUnitig = newUnitigPosition;
// 		lastKmer = current;
// 		lastPosition = kmerStartPos;
// 		addedAny = true;
// 	}
// 	if (lastPosition != 0)
// 	{
// 		assert(lastKmer != fromKmer);
// 		assert(lastUnitig != unitigPosition);
// 		assert(lastPosition != 0);
// 		assert(addedAny);
// 		size_t overlap = kmerSize - (edgeSequence.size() - kmerSize) + lastPosition;
// 		assert(overlap < kmerSize);
// 		hashlist.addSequenceOverlap(lastKmer, toKmer, overlap);
// 		unitigs.edgeCoverage(lastUnitig, toUnitigPosition) += edgeCoverage;
// 		auto key = canon(lastUnitig, toUnitigPosition);
// 		addEdges.insert(key);
// 	}
// 	return addedAny;
// }

// void forceEdgeDeterminism(HashList& reads, UnitigGraph& unitigs, const size_t kmerSize, const double minUnitigCoverage)
// {
// 	std::unordered_map<uint64_t, unsigned char> approxNeighbors;
// 	for (size_t i = 0; i < unitigs.edges.size(); i++)
// 	{
// 		std::pair<size_t, bool> fw { i, true };
// 		auto fromKmer = unitigs.unitigs[i].back();
// 		for (auto edge : unitigs.edges[fw])
// 		{
// 			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			collectApproxNeighbors(approxNeighbors, reads, kmerSize, fromKmer, toKmer);
// 		}
// 		std::pair<size_t, bool> bw { i, false };
// 		fromKmer = reverse(unitigs.unitigs[i][0]);
// 		for (auto edge : unitigs.edges[bw])
// 		{
// 			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			collectApproxNeighbors(approxNeighbors, reads, kmerSize, fromKmer, toKmer);
// 		}
// 	}
// 	std::unordered_set<uint64_t> approxJunctionHashes;
// 	for (size_t i = 0; i < unitigs.edges.size(); i++)
// 	{
// 		std::pair<size_t, bool> fw { i, true };
// 		auto fromKmer = unitigs.unitigs[i].back();
// 		for (auto edge : unitigs.edges[fw])
// 		{
// 			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			collectJunctionHashes(approxNeighbors, approxJunctionHashes, reads, kmerSize, fromKmer, toKmer);
// 		}
// 		std::pair<size_t, bool> bw { i, false };
// 		fromKmer = reverse(unitigs.unitigs[i][0]);
// 		for (auto edge : unitigs.edges[bw])
// 		{
// 			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			collectJunctionHashes(approxNeighbors, approxJunctionHashes, reads, kmerSize, fromKmer, toKmer);
// 		}
// 	}
// 	std::unordered_map<NodeType, std::pair<size_t, bool>> belongsToUnitig;
// 	std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> removeEdges;
// 	std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> addEdges;
// 	for (size_t i = 0; i < unitigs.edges.size(); i++)
// 	{
// 		std::pair<size_t, bool> fw { i, true };
// 		auto fromKmer = unitigs.unitigs[i].back();
// 		for (auto edge : unitigs.edges[fw])
// 		{
// 			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, addEdges, reads, unitigs, kmerSize, fw, edge, fromKmer, toKmer))
// 			{
// 				removeEdges.emplace(fw, edge);
// 			}
// 		}
// 		std::pair<size_t, bool> bw { i, false };
// 		fromKmer = reverse(unitigs.unitigs[i][0]);
// 		for (auto edge : unitigs.edges[bw])
// 		{
// 			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
// 			auto toKmer = unitigs.unitigs[edge.first][0];
// 			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
// 			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, addEdges, reads, unitigs, kmerSize, bw, edge, fromKmer, toKmer))
// 			{
// 				removeEdges.emplace(bw, edge);
// 			}
// 		}
// 	}
// 	assert((addEdges.size() > 0) == (removeEdges.size() > 0));
// 	bool addedAny = removeEdges.size() > 0;
// 	for (auto pair : addEdges)
// 	{
// 		unitigs.edges[pair.first].insert(pair.second);
// 		unitigs.edges[reverse(pair.second)].insert(reverse(pair.first));
// 	}
// 	for (auto pair : removeEdges)
// 	{
// 		if (unitigs.edges[pair.first].count(pair.second) == 1)
// 		{
// 			unitigs.edges[pair.first].erase(pair.second);
// 		}
// 		if (unitigs.edges[reverse(pair.second)].count(reverse(pair.first)) == 1)
// 		{
// 			unitigs.edges[reverse(pair.second)].erase(reverse(pair.first));
// 		}
// 	}
// 	if (addedAny)
// 	{
// 		unitigs = getUnitigs(unitigs);
// 		if (minUnitigCoverage > 1)
// 		{
// 			size_t oldSize = unitigs.unitigs.size();
// 			unitigs = unitigs.filterUnitigsByCoverage(minUnitigCoverage);
// 			if (unitigs.unitigs.size() != oldSize) unitigs = getUnitigs(unitigs);
// 		}
// 	}
// }

void printUnitigKmerCount(const UnitigGraph& unitigs)
{
	size_t unitigKmers = 0;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		unitigKmers += unitigs.unitigs[i].size();
	}
	std::cerr << unitigKmers << " distinct selected k-mers in unitigs after filtering" << std::endl;
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool blunt, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths)
{
	auto beforeReading = getTime();
	// check that all files actually exist
	for (const std::string& name : inputReads)
	{
		std::ifstream file { name };
		if (!file.good())
		{
			std::cerr << "Input file " << name << " can't be read!" << std::endl;
			std::exit(1);
		}
	}
	ReadpartIterator partIterator { kmerSize, windowSize, hpc };
	if (includeEndKmers)
	{
		// 2^32 arbitrarily, should be big enough for a human genome
		// a 30x coverage error-free hifi dataset will have an edge s-mer about every 250-300bp
		// so a human genome will have about 10000000-12000000 set bits
		// so about ~0.3% filled -> ~0.3% of kmers have a hash collision and are picked even though they are not edge k-mers
		// note that this also limits the effective window size to 250-300 regardless of what parameter is given
		std::vector<bool> endSmer;
		endSmer.resize(4294967296, false);
		collectEndSmers(endSmer, inputReads, kmerSize, windowSize, partIterator);
		std::swap(endSmer, partIterator.endSmers);
	}
	HashList reads { kmerSize };
	loadReadsAsHashesMultithread(reads, inputReads, kmerSize, partIterator, numThreads);
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage) unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	printUnitigKmerCount(unitigs);
	auto beforeSequences = getTime();
	std::cerr << "Getting unitig sequences" << std::endl;
	std::vector<std::pair<std::string, std::vector<uint16_t>>> unitigSequences;
	unitigSequences = getHPCUnitigSequences(reads, unitigs, inputReads, kmerSize, partIterator);
	auto beforeDeterminism = getTime();
	auto beforeConsistency = getTime();
	auto beforeWrite = getTime();
	AssemblyStats stats;
	if (blunt)
	{
		BluntGraph blunt { reads, unitigs, unitigSequences };
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(blunt, outputGraph);
	}
	else
	{
		//todo fix
		// if (windowSize > 1) forceEdgeDeterminism(reads, unitigs, unitigSequences, kmerSize, minUnitigCoverage);
		beforeConsistency = getTime();
		//todo fix
		// if (hpc) forceEdgeConsistency(unitigs, reads, unitigSequences, kmerSize);
		beforeWrite = getTime();
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(unitigs, outputGraph, reads, unitigSequences, kmerSize);
	}
	auto afterWrite = getTime();
	if (outputSequencePaths != "")
	{
		//todo fix
		assert(false);
		// assert(!blunt);
		// std::cerr << "Writing paths to " << outputSequencePaths << std::endl;
		// writePaths(reads, unitigs, inputReads, kmerSize, partIterator, outputSequencePaths, numThreads);
	}
	auto afterPaths = getTime();
	std::cerr << "reading and hashing graph topology took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeSequences) << std::endl;
	std::cerr << "building unitig sequences took " << formatTime(beforeSequences, beforeDeterminism) << std::endl;
	if (!blunt && windowSize > 1) std::cerr << "forcing edge determinism took " << formatTime(beforeDeterminism, beforeConsistency) << std::endl;
	if (!blunt && hpc) std::cerr << "forcing edge consistency took " << formatTime(beforeConsistency, beforeWrite) << std::endl;
	std::cerr << "writing the graph and calculating stats took " << formatTime(beforeWrite, afterWrite) << std::endl;
	if (outputSequencePaths != "") std::cerr << "writing sequence paths took " << formatTime(afterWrite, afterPaths) << std::endl;
	std::cerr << "nodes: " << stats.nodes << std::endl;
	std::cerr << "edges: " << stats.edges << std::endl;
	std::cerr << "assembly size " << stats.size << " bp, N50 " << stats.N50 << std::endl;
	if (stats.approxKmers > 0) std::cerr << "approximate number of k-mers ~ " << stats.approxKmers << std::endl;
	if (blunt && stats.nodes > 1 && stats.edges >= 2)
	{
		std::cerr << std::endl << "The output graph has redundant nodes and edges." << std::endl;
		std::cerr << "It is recommended to clean the graph using vg (https://github.com/vgteam/vg) with the command:" << std::endl << std::endl;
		std::cerr << "vg view -Fv " << outputGraph << " | vg mod -n -U 100 - | vg view - > cleaned." << outputGraph << std::endl << std::endl;
	}
}
