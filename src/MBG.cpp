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
#include <phmap.h>
#include <thread>
#include "fastqloader.h"
#include "CommonUtils.h"
#include "MBGCommon.h"
#include "MBG.h"
#include "VectorWithDirection.h"
#include "FastHasher.h"
#include "SparseEdgeContainer.h"
#include "HashList.h"
#include "UnitigGraph.h"
#include "BluntGraph.h"
#include "ReadHelper.h"
#include "HPCConsensus.h"
#include "ErrorMaskHelper.h"
#include "VectorView.h"
#include "StringIndex.h"

struct AssemblyStats
{
public:
	size_t nodes;
	size_t edges;
	size_t size;
	size_t N50;
	size_t approxKmers;
};

void collectEndSmers(std::vector<bool>& endSmer, const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const ReadpartIterator& partIterator)
{
	size_t smerSize = kmerSize - windowSize + 1;
	size_t addedEndSmers = 0;
	iterateReadsMultithreaded(files, 1, [&endSmer, smerSize, &addedEndSmers, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iterateParts(read, [&endSmer, smerSize, &addedEndSmers](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq) {
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
	iterateReadsMultithreaded(files, numThreads, [&result, &totalNodes, kmerSize, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &totalNodes, kmerSize](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, uint64_t minHash, const std::vector<size_t>& positions)
		{
			SequenceCharType revSeq = revCompRLE(seq);
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			HashType lastHash = 0;
			for (auto pos : positions)
			{
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				VectorView<CharType> revMinimizerSequence { revSeq, revPos, revPos + kmerSize };
				std::pair<size_t, bool> current;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				std::tie(current, lastHash) = result.addNode(minimizerSequence, revMinimizerSequence, lastHash, overlap, minHash);
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

void outputSequencePathLine(std::ofstream& outPaths, const std::string& readName, const std::vector<size_t>& poses, const std::string& rawSeq, const size_t currentSeqStart, const size_t currentSeqEnd, const std::vector<std::tuple<size_t, bool, size_t>>& kmerPath, const size_t currentPathStart, const size_t currentPathEnd, const std::vector<size_t>& unitigLength, const std::vector<size_t>& kmerUnitigStart, const std::vector<size_t>& kmerUnitigEnd, const UnitigGraph& unitigs, const HashList& hashlist)
{
	if (kmerPath.size() == 0) return;
	assert(poses.size() >= currentSeqEnd);
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
	size_t RLExpandedSeqSize = rawSeq.size();
	size_t RLExpandedSeqStart = poses[currentSeqStart];
	size_t RLExpandedSeqEnd = poses[currentSeqEnd];
	size_t RLExpandedPathSize = endAddition + unitigLength[lastUnitig.first];
	assert(RLExpandedSeqStart < RLExpandedSeqEnd);
	assert(RLExpandedSeqEnd <= RLExpandedSeqSize);
	assert(RLExpandedPathStart < RLExpandedPathEnd);
	assert(RLExpandedPathEnd <= RLExpandedPathSize);
	outPaths << readName << "\t" << RLExpandedSeqSize << "\t" << RLExpandedSeqStart << "\t" << RLExpandedSeqEnd << "\t" << "+" << "\t" << pathStr << "\t" << RLExpandedPathSize << "\t" << RLExpandedPathStart << "\t" << RLExpandedPathEnd << "\t" << (RLExpandedPathEnd - RLExpandedPathStart) << "\t" << (RLExpandedPathEnd - RLExpandedPathStart) << "\t" << "60" << std::endl;
}

std::unordered_set<uint64_t> collectApproxHashes(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const std::vector<std::tuple<size_t, bool, size_t>>& kmerUnitigOffset)
{
	std::unordered_set<uint64_t> approxHashes;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			FastHasher hasher { kmerSize };
			assert(std::get<0>(kmerUnitigOffset[unitigs.unitigs[i][j].first]) == i);
			std::tuple<size_t, bool, size_t> unitigpos = kmerUnitigOffset.at(unitigs.unitigs[i][j].first);
			size_t unitig = std::get<0>(unitigpos);
			bool fw = std::get<1>(unitigpos);
			size_t offset = std::get<2>(unitigpos);
			for (size_t k = 0; k < kmerSize; k++)
			{
				size_t pos = offset + k;
				if (!fw) pos = offset + kmerSize - 1 - k;
				int chr = unitigSequences[unitig].getCompressed(pos);
				if (!fw) chr = complement(chr);
				hasher.addChar(chr);
			}
			approxHashes.insert(hasher.hash());
		}
	}
	return approxHashes;
}

std::unordered_set<HashType> collectExactHashes(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const std::vector<std::tuple<size_t, bool, size_t>>& kmerUnitigOffset)
{
	std::unordered_set<HashType> exactHashes;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			assert(std::get<0>(kmerUnitigOffset[unitigs.unitigs[i][j].first]) == i);
			std::tuple<size_t, bool, size_t> unitigpos = kmerUnitigOffset.at(unitigs.unitigs[i][j].first);
			size_t unitig = std::get<0>(unitigpos);
			size_t offset = std::get<2>(unitigpos);
			SequenceCharType seq = unitigSequences[unitig].compressedSubstr(offset, kmerSize);
			assert(seq.size() == kmerSize);
			exactHashes.insert(hash(seq));
			seq = revCompRLE(seq);
			exactHashes.insert(hash(seq));
		}
	}
	return exactHashes;
}

template <typename F>
void findCollectedKmers(const SequenceCharType& seq, const size_t kmerSize, const std::unordered_set<uint64_t>& approxHashes, const std::unordered_set<HashType>& exactHashes, F callback)
{
	FastHasher hasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		hasher.addChar(seq[i]);
	}
	if (approxHashes.count(hasher.hash()) == 1)
	{
		VectorView view { seq, 0, kmerSize };
		HashType exactHash = hash(view);
		if (exactHashes.count(exactHash) == 1)
		{
			callback(0, exactHash);
		}
	}
	for (size_t i = kmerSize; i < seq.size(); i++)
	{
		hasher.addChar(seq[i]);
		hasher.removeChar(seq[i-kmerSize]);
		if (approxHashes.count(hasher.hash()) == 0) continue;
		VectorView view { seq, i - kmerSize + 1, i + 1 };
		HashType exactHash = hash(view);
		if (exactHashes.count(exactHash) == 0) continue;
		callback(i - kmerSize + 1, exactHash);
	}
}

void writePaths(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const std::vector<std::string>& inputReads, const size_t kmerSize, const ReadpartIterator& partIterator, const std::string& outputSequencePaths, const size_t numThreads)
{
	std::vector<std::tuple<size_t, bool, size_t>> kmerUnitigIndex;
	std::vector<std::tuple<size_t, bool, size_t>> kmerUnitigOffset;
	kmerUnitigIndex.resize(hashlist.size(), std::make_tuple(std::numeric_limits<size_t>::max(), true, std::numeric_limits<size_t>::max()));
	kmerUnitigOffset.resize(hashlist.size(), std::make_tuple(std::numeric_limits<size_t>::max(), true, std::numeric_limits<size_t>::max()));
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t offset = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0) offset += kmerSize - hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			assert(std::get<0>(kmerUnitigIndex[unitigs.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			std::get<0>(kmerUnitigIndex[unitigs.unitigs[i][j].first]) = i;
			std::get<1>(kmerUnitigIndex[unitigs.unitigs[i][j].first]) = unitigs.unitigs[i][j].second;
			std::get<0>(kmerUnitigOffset[unitigs.unitigs[i][j].first]) = i;
			std::get<1>(kmerUnitigOffset[unitigs.unitigs[i][j].first]) = unitigs.unitigs[i][j].second;
			std::get<2>(kmerUnitigOffset[unitigs.unitigs[i][j].first]) = offset;
			if (unitigs.unitigs[i][j].second)
			{
				std::get<2>(kmerUnitigIndex[unitigs.unitigs[i][j].first]) = j;
			}
			else
			{
				std::get<2>(kmerUnitigIndex[unitigs.unitigs[i][j].first]) = unitigs.unitigs[i].size() - 1 - j;
			}
		}
	}
	std::unordered_set<uint64_t> approxHashes = collectApproxHashes(hashlist, kmerSize, unitigs, unitigSequences, kmerUnitigOffset);
	std::unordered_set<HashType> exactHashes = collectExactHashes(hashlist, kmerSize, unitigs, unitigSequences, kmerUnitigOffset);
	std::vector<size_t> unitigLength;
	std::vector<size_t> kmerUnitigStart;
	std::vector<size_t> kmerUnitigEnd;
	unitigLength.resize(unitigs.unitigs.size(), std::numeric_limits<size_t>::max());
	kmerUnitigStart.resize(hashlist.size(), std::numeric_limits<size_t>::max());
	kmerUnitigEnd.resize(hashlist.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::vector<size_t> RLEposes = unitigSequences[i].getExpandedPositions(stringIndex);

		unitigLength[i] = RLEposes.back();
		size_t seqPos = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0)
			{
				size_t overlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
				assert(overlap > 0);
				assert(overlap < kmerSize);
				seqPos += kmerSize - overlap;
			}
			std::pair<size_t, bool> kmerHere = unitigs.unitigs[i][j];
			assert(kmerUnitigStart[kmerHere.first] == std::numeric_limits<size_t>::max());
			assert(kmerUnitigEnd[kmerHere.first] == std::numeric_limits<size_t>::max());
			assert(seqPos + kmerSize < RLEposes.size());
			if (kmerHere.second)
			{
				kmerUnitigStart[kmerHere.first] = RLEposes[seqPos];
				kmerUnitigEnd[kmerHere.first] = RLEposes[seqPos + kmerSize] - 1;
				assert(kmerUnitigStart[kmerHere.first] < kmerUnitigEnd[kmerHere.first]);
				assert(kmerUnitigEnd[kmerHere.first] < unitigLength[i]);
			}
			else
			{
				kmerUnitigEnd[kmerHere.first] = unitigLength[i] - RLEposes[seqPos] - 1;
				kmerUnitigStart[kmerHere.first] = unitigLength[i] - RLEposes[seqPos + kmerSize];
				assert(kmerUnitigStart[kmerHere.first] < kmerUnitigEnd[kmerHere.first]);
				assert(kmerUnitigEnd[kmerHere.first] < unitigLength[i]);
			}
		}
	}
	std::ofstream outPaths { outputSequencePaths };
	std::mutex pathWriteMutex;
	iterateReadsMultithreaded(inputReads, numThreads, [&outPaths, kmerSize, &approxHashes, &exactHashes, &kmerUnitigIndex, &unitigLength, &kmerUnitigStart, &kmerUnitigEnd, &hashlist, &unitigs, &pathWriteMutex, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iterateParts(read, [&read, &outPaths, kmerSize, &approxHashes, &exactHashes, &kmerUnitigIndex, &unitigLength, &kmerUnitigStart, &kmerUnitigEnd, &hashlist, &unitigs, &pathWriteMutex](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq) {
			if (seq.size() < kmerSize) return;
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
			std::vector<std::tuple<size_t, bool, size_t>> kmerPath;
			size_t currentSeqStart = std::numeric_limits<size_t>::max();
			size_t currentPathStart = std::numeric_limits<size_t>::max();
			size_t currentSeqEnd = std::numeric_limits<size_t>::max();
			size_t currentPathEnd = std::numeric_limits<size_t>::max();
			std::vector<std::pair<size_t, std::pair<size_t, bool>>> matches;
			findCollectedKmers(seq, kmerSize, approxHashes, exactHashes, [&matches, &hashlist](size_t pos, HashType hash)
			{
				matches.emplace_back(pos, hashlist.hashToNode.at(hash));
			});
			for (size_t j = matches.size()-1; j < matches.size(); j--)
			{
				for (size_t k = j+1; k < matches.size(); k++)
				{
					size_t overlap = kmerSize - (matches[k].first - matches[j].first);
					auto canonEdge = canon(matches[j].second, matches[k].second);
					if (hashlist.sequenceOverlap[canonEdge.first].count(canonEdge.second) == 0) continue;
					if (overlap != hashlist.getOverlap(canonEdge.first, canonEdge.second)) continue;
					if (k == j+1) break;
					matches.erase(matches.begin()+j+1, matches.begin()+k);
					break;
				}
			}
			for (auto pair : matches)
			{
				size_t pos = pair.first;
				std::pair<size_t, bool> current = pair.second;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				lastMinimizerPosition = pos;
				assert(current.first != std::numeric_limits<size_t>::max());
				assert(current.first < kmerUnitigIndex.size());
				std::tuple<size_t, bool, size_t> unitigPosition = kmerUnitigIndex[current.first];
				if (std::get<0>(unitigPosition) == std::numeric_limits<size_t>::max())
				{
					{
						std::lock_guard<std::mutex> guard { pathWriteMutex };
						outputSequencePathLine(outPaths, read.seq_id, poses, rawSeq, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
					}
					kmerPath.clear();
					currentSeqStart = std::numeric_limits<size_t>::max();
					currentSeqEnd = std::numeric_limits<size_t>::max();
					currentPathStart = std::numeric_limits<size_t>::max();
					currentPathEnd = std::numeric_limits<size_t>::max();
					lastKmer.first = std::numeric_limits<size_t>::max();
					continue;
				}
				if (!current.second)
				{
					std::get<1>(unitigPosition) = !std::get<1>(unitigPosition);
					if (std::get<0>(unitigPosition) != std::numeric_limits<size_t>::max())
					{
						std::get<2>(unitigPosition) = unitigs.unitigs[std::get<0>(unitigPosition)].size() - 1 - std::get<2>(unitigPosition);
					}
				}
				if (kmerPath.size() == 0)
				{
					currentSeqStart = pos;
					currentSeqEnd = pos + kmerSize;
					currentPathStart = kmerUnitigStart[current.first];
					if (!current.second)
					{
						currentPathStart = unitigLength[std::get<0>(unitigPosition)] - currentPathStart;
					}
					currentPathEnd = currentPathStart + kmerSize;
					kmerPath.emplace_back(unitigPosition);
					lastKmer = current;
					continue;
				}
				assert(lastKmer.first != std::numeric_limits<size_t>::max());
				bool hasMaybeValidEdge = true;
				auto canonEdge = canon(lastKmer, current);
				if (hashlist.sequenceOverlap[canonEdge.first].count(canonEdge.second) == 0)
				{
					hasMaybeValidEdge = false;
				}
				if (hasMaybeValidEdge && overlap != hashlist.getOverlap(lastKmer, current))
				{
					hasMaybeValidEdge = false;
				}
				if (hasMaybeValidEdge)
				{
					if (std::get<0>(unitigPosition) == std::get<0>(kmerPath.back()) && std::get<1>(unitigPosition) == std::get<1>(kmerPath.back()) && std::get<2>(unitigPosition) == std::get<2>(kmerPath.back()) + 1)
					{
						assert(overlap < kmerSize);
						currentSeqEnd += kmerSize - overlap;
						currentPathEnd += kmerSize - overlap;
						kmerPath.emplace_back(unitigPosition);
						lastKmer = current;
						continue;
					}
					if (std::get<2>(unitigPosition) == 0 && std::get<2>(kmerPath.back()) == unitigs.unitigs[std::get<0>(kmerPath.back())].size()-1)
					{
						std::pair<size_t, bool> fromUnitig = std::make_pair(std::get<0>(kmerPath.back()), std::get<1>(kmerPath.back()));
						std::pair<size_t, bool> toUnitig = std::make_pair(std::get<0>(unitigPosition), std::get<1>(unitigPosition));
						auto canonEdge = canon(fromUnitig, toUnitig);
						if (unitigs.edges[canonEdge.first].count(canonEdge.second) == 1)
						{
							// edge is valid
							assert(overlap < kmerSize);
							currentSeqEnd += kmerSize - overlap;
							currentPathEnd += kmerSize - overlap;
							kmerPath.emplace_back(unitigPosition);
							lastKmer = current;
							continue;
						}
					}
				}
				// no valid edge
				{
					std::lock_guard<std::mutex> guard { pathWriteMutex };
					outputSequencePathLine(outPaths, read.seq_id, poses, rawSeq, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
				}
				kmerPath.clear();
				currentSeqStart = pos;
				currentSeqEnd = pos + kmerSize;
				currentPathStart = kmerUnitigStart[current.first];
				if (!current.second)
				{
					currentPathStart = unitigLength[std::get<0>(unitigPosition)] - currentPathStart;
				}
				currentPathEnd = currentPathStart + kmerSize;
				kmerPath.emplace_back(unitigPosition);
				lastKmer = current;
			}
			{
				std::lock_guard<std::mutex> guard { pathWriteMutex };
				outputSequencePathLine(outPaths, read.seq_id, poses, rawSeq, currentSeqStart, currentSeqEnd, kmerPath, currentPathStart, currentPathEnd, unitigLength, kmerUnitigStart, kmerUnitigEnd, unitigs, hashlist);
			}
		});
	});
}

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
			result.edges[std::make_pair(currentUnitig, false)].emplace(currentUnitig, false);
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
				assert(unitigStart[curr]);
				auto from = belongsToUnitig.at(fw.first);
				bool prevFw = fw.second;
				auto to = belongsToUnitig.at(curr.first);
				bool currFw = curr.second;
				from.second = !(from.second ^ prevFw);
				to.second = !(to.second ^ currFw);
				result.edges[from].emplace(to);
				result.edges[reverse(to)].emplace(reverse(from));
				result.edgeCoverage(from, to) = oldgraph.edgeCoverage(fw, curr);
			}
		}
		auto bw = std::make_pair(i, false);
		if (unitigEnd[bw])
		{
			for (auto curr : oldgraph.edges.at(bw))
			{
				assert(unitigStart[curr]);
				auto from = belongsToUnitig.at(bw.first);
				bool prevFw = bw.second;
				auto to = belongsToUnitig.at(curr.first);
				bool currFw = curr.second;
				from.second = !(from.second ^ prevFw);
				to.second = !(to.second ^ currFw);
				result.edges[from].emplace(to);
				result.edges[reverse(to)].emplace(reverse(from));
				result.edgeCoverage(from, to) = oldgraph.edgeCoverage(bw, curr);
			}
		}
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

size_t getOverlapFromRLE(const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, std::pair<size_t, bool> fromUnitig, size_t rleOverlap)
{
	size_t overlap = 0;
	if (fromUnitig.second)
	{
		for (size_t i = 0; i < rleOverlap; i++)
		{
			overlap += unitigSequences[fromUnitig.first].getExpandedStr(unitigSequences[fromUnitig.first].compressedSize() - 1 - i, stringIndex).size();
		}
	}
	else
	{
		for (size_t i = 0; i < rleOverlap; i++)
		{
			overlap += unitigSequences[fromUnitig.first].getExpandedStr(i, stringIndex).size();
		}
	}
	return overlap;
}

AssemblyStats writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const size_t kmerSize)
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
		std::string realSequence = unitigSequences[i].getExpandedSequence(stringIndex);
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
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, fw, rleOverlap);
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
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, bw, rleOverlap);
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

std::tuple<NodeType, size_t, bool> find(std::vector<std::vector<std::tuple<NodeType, size_t, bool>>>& parent, std::pair<NodeType, size_t> key)
{
	while (true)
	{
		assert(key.first < parent.size());
		auto foundParent = parent[key.first][key.second];
		assert(std::get<0>(foundParent) < parent.size());
		if (std::get<0>(parent[key.first][key.second]) == std::get<0>(parent[std::get<0>(foundParent)][std::get<1>(foundParent)]) && std::get<1>(parent[key.first][key.second]) == std::get<1>(parent[std::get<0>(foundParent)][std::get<1>(foundParent)]))
		{
			assert(std::get<2>(parent[std::get<0>(foundParent)][std::get<1>(foundParent)]) == true);
			return parent[key.first][key.second];
		}
		bool keepDirection = std::get<2>(parent[key.first][key.second]);
		parent[key.first][key.second] = parent[std::get<0>(foundParent)][std::get<1>(foundParent)];
		if (!keepDirection) std::get<2>(parent[key.first][key.second]) = !std::get<2>(parent[key.first][key.second]);
	}
}

void merge(std::vector<std::vector<std::tuple<NodeType, size_t, bool>>>& parent, std::vector<std::vector<size_t>>& rank, std::pair<NodeType, size_t> left, std::pair<NodeType, size_t> right, bool keepDirection)
{
	auto foundLeft = find(parent, left);
	auto foundRight = find(parent, right);
	left.first = std::get<0>(foundLeft);
	left.second = std::get<1>(foundLeft);
	right.first = std::get<0>(foundRight);
	right.second = std::get<1>(foundRight);
	bool keepOnce = std::get<2>(foundLeft) == std::get<2>(foundRight);
	std::get<2>(foundLeft) = true;
	std::get<2>(foundRight) = true;
	assert(std::get<0>(parent[left.first][left.second]) == left.first);
	assert(std::get<1>(parent[left.first][left.second]) == left.second);
	assert(std::get<2>(parent[left.first][left.second]) == true);
	assert(std::get<0>(parent[right.first][right.second]) == right.first);
	assert(std::get<1>(parent[right.first][right.second]) == right.second);
	assert(std::get<2>(parent[right.first][right.second]) == true);
	if (left == right) return;
	if (rank[right.first][right.second] > rank[left.first][left.second])
	{
		std::swap(left, right);
		std::swap(foundLeft, foundRight);
	}
	parent[right.first][right.second] = foundLeft;
	if (!keepOnce) std::get<2>(parent[right.first][right.second]) = !std::get<2>(parent[right.first][right.second]);
	if (!keepDirection) std::get<2>(parent[right.first][right.second]) = !std::get<2>(parent[right.first][right.second]);
	assert(rank[right.first][right.second] <= rank[left.first][left.second]);
	if (rank[right.first][right.second] == rank[left.first][left.second]) rank[left.first][left.second] += 1;
}

void merge(std::vector<std::vector<std::tuple<NodeType, size_t, bool>>>& parent, std::vector<std::vector<size_t>>& rank, NodeType leftNode, size_t leftOffset, NodeType rightNode, size_t rightOffset, bool keepDirection)
{
	merge(parent, rank, std::make_pair(leftNode, leftOffset), std::make_pair(rightNode, rightOffset), keepDirection);
}

void forceEdgeConsistency(const UnitigGraph& unitigs, HashList& hashlist, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize)
{
	std::vector<std::vector<std::tuple<NodeType, size_t, bool>>> parent;
	std::vector<std::vector<size_t>> rank;
	parent.resize(unitigSequences.size());
	rank.resize(unitigSequences.size());
	for (size_t i = 0; i < unitigSequences.size(); i++)
	{
		if (unitigSequences[i].compressedSize() >= 2 * kmerSize)
		{
			for (size_t j = 0; j < 2 * kmerSize; j++)
			{
				parent[i].emplace_back(i, j, true);
				rank[i].emplace_back(0);
			}
		}
		else
		{
			for (size_t j = 0; j < unitigSequences[i].compressedSize(); j++)
			{
				parent[i].emplace_back(i, j, true);
				rank[i].emplace_back(0);
			}
		}
	}
	for (size_t unitig = 0; unitig < unitigs.edges.size(); unitig++)
	{
		std::pair<size_t, bool> fw { unitig, true };
		std::pair<size_t, bool> bw { unitig, false };
		for (auto target : unitigs.edges[fw])
		{
			std::pair<NodeType, bool> from = unitigs.unitigs[unitig].back();
			std::pair<NodeType, bool> to = unitigs.unitigs[target.first][0];
			if (!target.second) to = reverse(unitigs.unitigs[target.first].back());
			size_t overlap = hashlist.getOverlap(from, to);
			for (size_t i = 0; i < overlap; i++)
			{
				size_t firstOffset = parent[unitig].size() - overlap + i;
				size_t secondOffset = i;
				if (!target.second) secondOffset = parent[target.first].size() - 1 - i;
				merge(parent, rank, unitig, firstOffset, target.first, secondOffset, fw.second == target.second);
			}
		}
		for (auto target : unitigs.edges[bw])
		{
			std::pair<NodeType, bool> from = reverse(unitigs.unitigs[unitig][0]);
			std::pair<NodeType, bool> to = unitigs.unitigs[target.first][0];
			if (!target.second) to = reverse(unitigs.unitigs[target.first].back());
			size_t overlap = hashlist.getOverlap(from, to);
			for (size_t i = 0; i < overlap; i++)
			{
				size_t firstOffset = overlap - 1 - i;
				size_t secondOffset = i;
				if (!target.second) secondOffset = parent[target.first].size() - 1 - i;
				merge(parent, rank, unitig, firstOffset, target.first, secondOffset, bw.second == target.second);
			}
		}
	}
	for (size_t i = 0; i < parent.size(); i++)
	{
		for (size_t j = 0; j < parent[i].size(); j++)
		{
			auto found = find(parent, std::make_pair(i, j));
			size_t secondPos = std::get<1>(found);
			size_t secondUnitig = std::get<0>(found);
			assert(secondPos < 2 * kmerSize);
			if (secondPos >= kmerSize)
			{
				size_t secondFromEnd = parent[secondUnitig].size() - 1 - secondPos;
				assert(secondFromEnd < kmerSize);
				secondPos = unitigSequences[secondUnitig].compressedSize() - 1 - secondFromEnd;
			}
			size_t firstPos = j;
			assert(firstPos < 2 * kmerSize);
			if (firstPos >= kmerSize)
			{
				size_t firstFromEnd = parent[i].size() - 1 - firstPos;
				assert(firstFromEnd < kmerSize);
				firstPos = unitigSequences[i].compressedSize() - 1 - firstFromEnd;
			}
			uint32_t seq = unitigSequences[secondUnitig].getExpanded(secondPos);
			if (std::get<2>(found))
			{
				assert(unitigSequences[i].getCompressed(firstPos) == unitigSequences[secondUnitig].getCompressed(secondPos));
			}
			else
			{
				assert(unitigSequences[i].getCompressed(firstPos) == complement(unitigSequences[secondUnitig].getCompressed(secondPos)));
			}
			assert(secondUnitig < unitigSequences.size());
			assert(i < unitigSequences.size());
			unitigSequences[i].setExpanded(firstPos, seq);
		}
	}
}

SequenceCharType getEdgeSequence(const std::vector<CompressedSequenceType>& unitigSequences, const std::pair<size_t, bool> fromUnitig, const std::pair<size_t, bool> toUnitig, const size_t kmerSize, const size_t overlap)
{
	SequenceCharType edgeSequence;
	if (fromUnitig.second)
	{
		edgeSequence = unitigSequences[fromUnitig.first].compressedSubstr(unitigSequences[fromUnitig.first].compressedSize() - kmerSize, kmerSize);
	}
	else
	{
		edgeSequence = unitigSequences[fromUnitig.first].compressedSubstr(0, kmerSize);
		edgeSequence = revCompRLE(edgeSequence);
	}
	if (toUnitig.second)
	{
		SequenceCharType add = unitigSequences[toUnitig.first].compressedSubstr(overlap, kmerSize - overlap);
		edgeSequence.insert(edgeSequence.end(), add.begin(), add.end());
	}
	else
	{
		SequenceCharType add = unitigSequences[toUnitig.first].compressedSubstr(unitigSequences[toUnitig.first].compressedSize() - kmerSize, kmerSize - overlap);
		add = revCompRLE(add);
		edgeSequence.insert(edgeSequence.end(), add.begin(), add.end());
	}
	assert(edgeSequence.size() > kmerSize);
	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
	assert(edgeSequence.size() < 2 * kmerSize);
	return edgeSequence;
}

void collectApproxNeighbors(std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize, std::pair<NodeType, bool> fromUnitig, std::pair<NodeType, bool> toUnitig, const size_t overlap)
{
	SequenceCharType edgeSequence = getEdgeSequence(unitigSequences, fromUnitig, toUnitig, kmerSize, overlap);
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(edgeSequence[i]);
	}
	for (size_t i = kmerSize; i < edgeSequence.size(); i++)
	{
		auto oldFwHash = fwkmerHasher.getFwHash();
		fwkmerHasher.addChar(edgeSequence[i]);
		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
		auto newBwHash = fwkmerHasher.getBwHash();
		if (approxNeighbors.count(oldFwHash) == 1 && approxNeighbors.at(oldFwHash) != edgeSequence[i]) approxNeighbors[oldFwHash] = std::numeric_limits<unsigned char>::max();
		if (approxNeighbors.count(oldFwHash) == 0) approxNeighbors[oldFwHash] = edgeSequence[i];
		if (approxNeighbors.count(newBwHash) == 1 && approxNeighbors.at(newBwHash) != complement(edgeSequence[i-kmerSize])) approxNeighbors[newBwHash] = std::numeric_limits<unsigned char>::max();
		if (approxNeighbors.count(newBwHash) == 0) approxNeighbors[newBwHash] = complement(edgeSequence[i-kmerSize]);
	}
}

void collectJunctionHashes(const std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::unordered_set<uint64_t>& approxHashes, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize, std::pair<NodeType, bool> fromUnitig, std::pair<NodeType, bool> toUnitig, const size_t overlap)
{
	SequenceCharType edgeSequence = getEdgeSequence(unitigSequences, fromUnitig, toUnitig, kmerSize, overlap);
	if (edgeSequence.size() == kmerSize + 1) return;
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(edgeSequence[i]);
	}
	for (size_t i = kmerSize; i < edgeSequence.size(); i++)
	{
		auto oldFwHash = fwkmerHasher.getFwHash();
		auto oldHash = fwkmerHasher.hash();
		fwkmerHasher.addChar(edgeSequence[i]);
		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
		auto newBwHash = fwkmerHasher.getBwHash();
		assert(approxNeighbors.count(oldFwHash) == 1);
		assert(approxNeighbors.count(newBwHash) == 1);
		// neither is a junction node if both not in approx
		if (approxNeighbors.at(oldFwHash) != std::numeric_limits<unsigned char>::max() && approxNeighbors.at(newBwHash) != std::numeric_limits<unsigned char>::max()) continue;
		auto newHash = fwkmerHasher.hash();
		approxHashes.insert(oldHash);
		approxHashes.insert(newHash);
	}
}

bool addJunctionsToHashes(const std::unordered_set<uint64_t>& approxHashes, std::unordered_map<NodeType, std::pair<size_t, bool>>& belongsToUnitig, HashList& hashlist, UnitigGraph& newUnitigs, std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>>& addedEdges, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize, const std::pair<size_t, bool> fromUnitig, const std::pair<size_t, bool> toUnitig, const std::pair<size_t, bool> fromKmer, const std::pair<size_t, bool> toKmer, const size_t edgeCoverage, const size_t unitigOverlap, const size_t zeroKmer, std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>>& addedKmerSequences)
{
	SequenceCharType edgeSequence = getEdgeSequence(unitigSequences, fromUnitig, toUnitig, kmerSize, unitigOverlap);
	if (edgeSequence.size() == kmerSize + 1) return false;
	SequenceCharType revSeq = revCompRLE(edgeSequence);
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(edgeSequence[i]);
	}
	std::pair<NodeType, bool> lastKmer = fromKmer;
	std::pair<NodeType, bool> lastUnitig = fromUnitig;
	size_t lastPosition = 0;
	bool addedAny = false;
	assert(edgeCoverage > 0);
	for (size_t i = kmerSize; i < edgeSequence.size()-1; i++)
	{
		fwkmerHasher.addChar(edgeSequence[i]);
		fwkmerHasher.removeChar(edgeSequence[i-kmerSize]);
		if (approxHashes.count(fwkmerHasher.hash()) == 0) continue;
		size_t kmerStartPos = i - kmerSize + 1;

		size_t overlap = kmerSize - kmerStartPos + lastPosition;
		assert(overlap < kmerSize);
		HashType lastHash = 0;
		uint64_t minHash = 0;

		VectorView<CharType> minimizerSequence { edgeSequence, kmerStartPos, kmerStartPos + kmerSize };
		size_t revPos = edgeSequence.size() - (kmerStartPos + kmerSize);
		VectorView<CharType> revMinimizerSequence { revSeq, revPos, revPos + kmerSize };
		std::pair<size_t, bool> current;
		std::tie(current, lastHash) = hashlist.addNode(minimizerSequence, revMinimizerSequence, lastHash, overlap, minHash);
		assert(current.first >= zeroKmer);
		assert(current.first <= zeroKmer + addedKmerSequences.size());
		if (current.first == zeroKmer + addedKmerSequences.size())
		{
			assert(current.second);
			addedKmerSequences.emplace_back();
			std::get<0>(addedKmerSequences.back()) = fromKmer;
			std::get<1>(addedKmerSequences.back()) = toKmer;
			std::get<2>(addedKmerSequences.back()) = kmerStartPos;
		}
		hashlist.addSequenceOverlap(lastKmer, current, overlap);
		std::pair<size_t, bool> newUnitigPosition { std::numeric_limits<size_t>::max(), true };
		if (belongsToUnitig.count(current.first) == 0)
		{
			newUnitigs.unitigs.emplace_back();
			newUnitigs.unitigs.back().emplace_back(current.first, true);
			newUnitigs.unitigCoverage.emplace_back();
			newUnitigs.unitigCoverage.back().emplace_back(0);
			newUnitigs.unitigCoverage.back()[0] += edgeCoverage;
			newUnitigs.edges.emplace_back();
			newUnitigs.edgeCov.emplace_back();
			belongsToUnitig[current.first] = std::make_pair(newUnitigs.unitigs.size()-1, true);
			newUnitigPosition = std::make_pair(newUnitigs.unitigs.size()-1, true);
			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
		}
		else
		{
			newUnitigPosition = belongsToUnitig.at(current.first);
			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
			assert(newUnitigs.unitigCoverage[newUnitigPosition.first].size() == 1);
			newUnitigs.unitigCoverage[newUnitigPosition.first][0] += edgeCoverage;
		}
		assert(newUnitigPosition.first != std::numeric_limits<size_t>::max());
		if (lastUnitig.first != std::numeric_limits<size_t>::max())
		{
			auto key = canon(lastUnitig, newUnitigPosition);
			addedEdges.emplace_back(key.first, key.second, edgeCoverage);
		}
		lastUnitig = newUnitigPosition;
		lastKmer = current;
		lastPosition = kmerStartPos;
		addedAny = true;
	}
	if (lastPosition != 0)
	{
		addedEdges.emplace_back(lastUnitig, toUnitig, edgeCoverage);
		size_t kmerStartPos = edgeSequence.size() - 1 - kmerSize + 1;
		size_t overlap = kmerSize - kmerStartPos + lastPosition;
		assert(overlap < kmerSize);
		hashlist.addSequenceOverlap(lastKmer, toKmer, overlap);
	}
	return addedAny;
}

CompressedSequenceType getKmerSequence(const std::vector<CompressedSequenceType>& unitigSequences, const std::vector<std::tuple<size_t, bool, size_t>>& kmerSequencePosition, const std::pair<size_t, bool> kmer, const size_t kmerSize)
{
	size_t unitig = std::get<0>(kmerSequencePosition[kmer.first]);
	bool fw = std::get<1>(kmerSequencePosition[kmer.first]) == kmer.second;
	size_t offset = std::get<2>(kmerSequencePosition[kmer.first]);
	assert(unitig != std::numeric_limits<size_t>::max());
	CompressedSequenceType result = unitigSequences[unitig].substr(offset, kmerSize);
	if (!fw)
	{
		result = result.revComp();
	}
	return result;
}

void forceEdgeDeterminism(HashList& reads, UnitigGraph& unitigs, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize, const double minUnitigCoverage)
{
	std::unordered_map<uint64_t, unsigned char> approxNeighbors;
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		auto fromKmer = unitigs.unitigs[i].back();
		for (auto edge : unitigs.edges[fw])
		{
			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			collectApproxNeighbors(approxNeighbors, unitigSequences, kmerSize, fw, edge, overlap);
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			collectApproxNeighbors(approxNeighbors, unitigSequences, kmerSize, bw, edge, overlap);
		}
	}
	std::unordered_set<uint64_t> approxJunctionHashes;
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		auto fromKmer = unitigs.unitigs[i].back();
		for (auto edge : unitigs.edges[fw])
		{
			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			collectJunctionHashes(approxNeighbors, approxJunctionHashes, unitigSequences, kmerSize, fw, edge, overlap);
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			collectJunctionHashes(approxNeighbors, approxJunctionHashes, unitigSequences, kmerSize, bw, edge, overlap);
		}
	}
	std::vector<std::tuple<size_t, bool, size_t>> kmerSequencePosition;
	kmerSequencePosition.resize(reads.coverage.size(), std::make_tuple(std::numeric_limits<size_t>::max(), true, 0));
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t offset = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0)
			{
				size_t overlap = reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
				offset += kmerSize - overlap;
			}
			assert(std::get<0>(kmerSequencePosition[unitigs.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			kmerSequencePosition[unitigs.unitigs[i][j].first] = std::make_tuple(i, unitigs.unitigs[i][j].second, offset);
		}
	}
	std::unordered_map<NodeType, std::pair<size_t, bool>> belongsToUnitig;
	std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> removeEdges;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> addedEdges;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> addedKmerSequences;
	size_t zeroKmer = reads.size();
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		auto fromKmer = unitigs.unitigs[i].back();
		for (auto edge : unitigs.edges[fw])
		{
			if (canon(fw, edge) != std::make_pair(fw, edge))
			{
				assert(unitigs.edges[reverse(edge)].count(reverse(fw)) == 1);
				continue;
			}
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, reads, unitigs, addedEdges, unitigSequences, kmerSize, fw, edge, fromKmer, toKmer, unitigs.edgeCoverage(fw, edge), overlap, zeroKmer, addedKmerSequences))
			{
				removeEdges.emplace(fw, edge);
			}
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge))
			{
				assert(unitigs.edges[reverse(edge)].count(reverse(bw)) == 1);
				continue;
			}
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			size_t overlap = reads.getOverlap(fromKmer, toKmer);
			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, reads, unitigs, addedEdges, unitigSequences, kmerSize, bw, edge, fromKmer, toKmer, unitigs.edgeCoverage(bw, edge), overlap, zeroKmer, addedKmerSequences))
			{
				removeEdges.emplace(bw, edge);
			}
		}
	}
	assert((addedEdges.size() == 0) == (removeEdges.size() == 0));
	if (removeEdges.size() == 0) return;
	for (auto e : removeEdges)
	{
		unitigs.edges[e.first].erase(e.second);
		unitigs.edges[reverse(e.second)].erase(reverse(e.first));
	}
	for (auto t : addedEdges)
	{
		std::pair<size_t, bool> from = std::get<0>(t);
		std::pair<size_t, bool> to = std::get<1>(t);
		size_t coverage = std::get<2>(t);
		unitigs.edges[from].insert(to);
		unitigs.edges[reverse(to)].insert(reverse(from));
		unitigs.edgeCoverage(from, to) += coverage;
	}
	unitigs = getUnitigs(unitigs);
	if (minUnitigCoverage > 1)
	{
		size_t oldSize = unitigs.unitigs.size();
		unitigs = unitigs.filterUnitigsByCoverage(minUnitigCoverage);
		if (unitigs.unitigs.size() != oldSize) unitigs = getUnitigs(unitigs);
	}
	std::vector<CompressedSequenceType> newUnitigSequences;
	newUnitigSequences.resize(unitigs.unitigs.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t length = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			size_t overlap = 0;
			if (j > 0) overlap = reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			length += kmerSize - overlap;
		}
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			CompressedSequenceType seqHere;
			if (unitigs.unitigs[i][j].first < zeroKmer)
			{
				seqHere = getKmerSequence(unitigSequences, kmerSequencePosition, unitigs.unitigs[i][j], kmerSize);
			}
			else
			{
				size_t k = unitigs.unitigs[i][j].first - zeroKmer;
				bool fw = unitigs.unitigs[i][j].second;
				std::pair<size_t, bool> fromKmer = std::get<0>(addedKmerSequences[k]);
				std::pair<size_t, bool> toKmer = std::get<1>(addedKmerSequences[k]);
				seqHere = getKmerSequence(unitigSequences, kmerSequencePosition, fromKmer, kmerSize);
				CompressedSequenceType otherSeq = getKmerSequence(unitigSequences, kmerSequencePosition, toKmer, kmerSize);
				size_t overlap = reads.getOverlap(fromKmer, toKmer);
				seqHere.insertEnd(otherSeq.substr(overlap, otherSeq.compressedSize()-overlap));
				size_t seqStart = std::get<2>(addedKmerSequences[k]);
				seqHere = seqHere.substr(seqStart, kmerSize);
				if (!fw)
				{
					seqHere = seqHere.revComp();
				}
			}
			size_t overlap = 0;
			if (j > 0) overlap = reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			newUnitigSequences[i].insertEnd(seqHere.substr(overlap, seqHere.compressedSize() - overlap));
		}
		assert(newUnitigSequences[i].compressedSize() == length);
	}
	std::swap(unitigSequences, newUnitigSequences);
}

void printUnitigKmerCount(const UnitigGraph& unitigs)
{
	size_t unitigKmers = 0;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		unitigKmers += unitigs.unitigs[i].size();
	}
	std::cerr << unitigKmers << " distinct selected k-mers in unitigs after filtering" << std::endl;
}

HashList filterHashes(const HashList& reads, const size_t kmerSize, const std::vector<size_t>& newIndex, size_t maxCount)
{
	HashList resultHashlist { kmerSize };
	resultHashlist.coverage.resize(maxCount, 0);
	resultHashlist.sequenceOverlap.resize(maxCount);
	resultHashlist.edgeCoverage.resize(maxCount);
	for (auto pair : reads.hashToNode)
	{
		if (newIndex[pair.second.first] == std::numeric_limits<size_t>::max()) continue;
		resultHashlist.hashToNode[pair.first] = std::make_pair(newIndex[pair.second.first], pair.second.second);
	}
	assert(reads.coverage.size() == reads.sequenceOverlap.size());
	assert(reads.coverage.size() == reads.edgeCoverage.size());
	for (size_t i = 0; i < reads.coverage.size(); i++)
	{
		if (newIndex[i] == std::numeric_limits<size_t>::max()) continue;
		size_t newI = newIndex[i];
		resultHashlist.coverage[newI] = reads.coverage[i];
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> newFw { newI, true };
		for (auto pair : reads.sequenceOverlap[fw])
		{
			std::pair<size_t, bool> to = pair.first;
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, bool> newTo { newIndex[to.first], to.second };
			resultHashlist.addSequenceOverlap(newFw, newTo, pair.second);
		}
		for (auto pair : reads.edgeCoverage[fw])
		{
			std::pair<size_t, bool> to = pair.first;
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, bool> newTo { newIndex[to.first], to.second };
			auto key = canon(newFw, newTo);
			resultHashlist.edgeCoverage[key.first][key.second] = pair.second;
		}
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newBw { newI, false };
		for (auto pair : reads.sequenceOverlap[bw])
		{
			std::pair<size_t, bool> to = pair.first;
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, bool> newTo { newIndex[to.first], to.second };
			resultHashlist.addSequenceOverlap(newBw, newTo, pair.second);
		}
		for (auto pair : reads.edgeCoverage[bw])
		{
			std::pair<size_t, bool> to = pair.first;
			if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, bool> newTo { newIndex[to.first], to.second };
			auto key = canon(newBw, newTo);
			resultHashlist.edgeCoverage[key.first][key.second] = pair.second;
		}
	}
	return resultHashlist;
}

std::pair<UnitigGraph, HashList> filterKmersToUnitigKmers(const UnitigGraph& unitigs, const HashList& reads, const size_t kmerSize)
{
	UnitigGraph resultUnitigs;
	resultUnitigs.edges = unitigs.edges;
	resultUnitigs.edgeCov = unitigs.edgeCov;
	resultUnitigs.unitigs.resize(unitigs.unitigs.size());
	resultUnitigs.unitigCoverage.resize(unitigs.unitigCoverage.size());
	std::vector<size_t> newIndex;
	newIndex.resize(reads.coverage.size(), std::numeric_limits<size_t>::max());
	size_t count = 0;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> newOverlaps;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t lastStart = 0;
		size_t currentPos = 0;
		std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0) currentPos += kmerSize - reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			bool skip = true;
			if (j == 0 || j == unitigs.unitigs[i].size()-1)
			{
				skip = false;
			}
			else
			{
				assert(j+1 < unitigs.unitigs[i].size());
				if (currentPos + (kmerSize - reads.getOverlap(unitigs.unitigs[i][j], unitigs.unitigs[i][j+1])) >= lastStart + kmerSize)
				{
					skip = false;
				}
			}
			assert(newIndex[unitigs.unitigs[i][j].first] == std::numeric_limits<size_t>::max());
			if (!skip)
			{
				newIndex[unitigs.unitigs[i][j].first] = count;
				count += 1;
				std::pair<size_t, bool> thisKmer { newIndex[unitigs.unitigs[i][j].first], unitigs.unitigs[i][j].second };
				if (lastKmer.first != std::numeric_limits<size_t>::max())
				{
					assert(currentPos > lastStart);
					assert(currentPos - lastStart < kmerSize);
					newOverlaps.emplace_back(lastKmer, thisKmer, kmerSize - (currentPos - lastStart));
				}
				resultUnitigs.unitigs[i].push_back(thisKmer);
				resultUnitigs.unitigCoverage[i].push_back(unitigs.unitigCoverage[i][j]);
				lastStart = currentPos;
				lastKmer = thisKmer;
			}
		}
	}
	HashList resultHashlist = filterHashes(reads, kmerSize, newIndex, count);
	for (auto t : newOverlaps)
	{
		resultHashlist.addSequenceOverlap(std::get<0>(t), std::get<1>(t), std::get<2>(t));
	}
	return std::make_pair(resultUnitigs, resultHashlist);
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const ErrorMasking errorMasking, const bool blunt, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths)
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
	ReadpartIterator partIterator { kmerSize, windowSize, errorMasking };
	if (includeEndKmers)
	{
		std::cerr << "Collecting end k-mers" << std::endl;
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
	std::cerr << "Collecting selected k-mers" << std::endl;
	loadReadsAsHashesMultithread(reads, inputReads, kmerSize, partIterator, numThreads);
	auto beforeUnitigs = getTime();
	std::cerr << "Unitigifying" << std::endl;
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage)
	{
		std::cerr << "Filtering by unitig coverage" << std::endl;
		unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	}
	std::tie(unitigs, reads) = filterKmersToUnitigKmers(unitigs, reads, kmerSize);
	printUnitigKmerCount(unitigs);
	auto beforeSequences = getTime();
	std::cerr << "Getting unitig sequences" << std::endl;
	std::vector<CompressedSequenceType> unitigSequences;
	StringIndex stringIndex;
	std::tie(unitigSequences, stringIndex) = getHPCUnitigSequences(reads, unitigs, inputReads, kmerSize, partIterator, numThreads);
	assert(unitigSequences.size() == unitigs.unitigs.size());
	auto beforeDeterminism = getTime();
	auto beforeConsistency = getTime();
	auto beforeWrite = getTime();
	AssemblyStats stats;
	if (blunt)
	{
		BluntGraph blunt { reads, unitigs, unitigSequences, stringIndex };
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(blunt, outputGraph);
	}
	else
	{
		if (windowSize > 1)
		{
			std::cerr << "Determinizing edges" << std::endl;
			forceEdgeDeterminism(reads, unitigs, unitigSequences, kmerSize, minUnitigCoverage);
		}
		beforeConsistency = getTime();
		if (errorMasking != ErrorMasking::No && errorMasking != ErrorMasking::Collapse)
		{
			std::cerr << "Forcing edge consistency" << std::endl;
			forceEdgeConsistency(unitigs, reads, unitigSequences, kmerSize);
		}
		beforeWrite = getTime();
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(unitigs, outputGraph, reads, unitigSequences, stringIndex, kmerSize);
	}
	auto afterWrite = getTime();
	if (outputSequencePaths != "")
	{
		assert(!blunt);
		std::cerr << "Writing paths to " << outputSequencePaths << std::endl;
		writePaths(reads, unitigs, unitigSequences, stringIndex, inputReads, kmerSize, partIterator, outputSequencePaths, numThreads);
	}
	auto afterPaths = getTime();
	std::cerr << "selecting k-mers and building graph topology took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeSequences) << std::endl;
	std::cerr << "building unitig sequences took " << formatTime(beforeSequences, beforeDeterminism) << std::endl;
	if (!blunt && windowSize > 1) std::cerr << "forcing edge determinism took " << formatTime(beforeDeterminism, beforeConsistency) << std::endl;
	if (!blunt && errorMasking != ErrorMasking::No && errorMasking != ErrorMasking::Collapse) std::cerr << "forcing edge consistency took " << formatTime(beforeConsistency, beforeWrite) << std::endl;
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
