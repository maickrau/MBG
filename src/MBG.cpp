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
#include "RankBitvector.h"
#include "UnitigResolver.h"

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

std::vector<HashPath> loadReadsAsHashesMultithread(HashList& result, const std::vector<std::string>& files, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	std::atomic<size_t> totalNodes = 0;
	std::vector<HashPath> paths;
	// todo multithreading with paths
	iterateReadsMultithreaded(files, 1, [&result, &totalNodes, kmerSize, &partIterator, &paths](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &totalNodes, kmerSize, &paths, read](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, uint64_t minHash, const std::vector<size_t>& positions)
		{
			SequenceCharType revSeq = revCompRLE(seq);
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			paths.emplace_back();
			paths.back().readName = read.seq_id;
			for (auto pos : positions)
			{
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				VectorView<CharType> revMinimizerSequence { revSeq, revPos, revPos + kmerSize };
				std::pair<size_t, bool> current;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				HashType hash;
				std::tie(current, hash) = result.addNode(minimizerSequence, revMinimizerSequence);
				paths.back().hashes.emplace_back(hash);
				paths.back().hashPoses.emplace_back(pos);
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
			if (paths.back().hashes.size() == 0) paths.pop_back();
		});
	});
	std::cerr << totalNodes << " total selected k-mers in reads" << std::endl;
	std::cerr << result.size() << " distinct selected k-mers in reads" << std::endl;
	std::cerr << paths.size() << " input paths" << std::endl;
	return paths;
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
				matches.emplace_back(pos, hashlist.getHashNode(hash));
			});
			for (size_t j = matches.size()-1; j < matches.size(); j--)
			{
				for (size_t k = j+1; k < matches.size(); k++)
				{
					size_t overlap = kmerSize - (matches[k].first - matches[j].first);
					auto canonEdge = canon(matches[j].second, matches[k].second);
					if (!hashlist.hasSequenceOverlap(canonEdge.first,canonEdge.second)) continue;
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
				if (!hashlist.hasSequenceOverlap(canonEdge.first, canonEdge.second))
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
	result.edgeOvlp.emplace_back();
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
			result.edgeOverlap(currentUnitig, true, currentUnitig, true) = old.edgeOverlap(pos.first, pos.second, newPos.first, newPos.second);
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
			result.edgeOverlap(currentUnitig, belongsToUnitig.at(pos.first).second, currentUnitig, !belongsToUnitig.at(pos.first).second) = old.edgeOverlap(pos.first, pos.second, newPos.first, newPos.second);
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

std::vector<std::pair<size_t, bool>> getUnitigHashes(const std::pair<size_t, bool> start, const SparseEdgeContainer& edges)
{
	std::vector<std::pair<size_t, bool>> result;
	result.emplace_back(start);
	while (true)
	{
		if (edges[result.back()].size() != 1) break;
		auto newPos = edges[result.back()][0];
		auto revPos = std::make_pair(newPos.first, !newPos.second);
		if (edges[revPos].size() != 1) break;
		if (newPos == start) break;
		if (newPos.first == result.back().first) break;
		result.emplace_back(newPos);
	}
	return result;
}

void checkUnitigHashes(const std::pair<size_t, bool> start, const SparseEdgeContainer& edges, std::vector<bool>& checked, const HashList& hashlist, const double minUnitigCoverage, RankBitvector& kept)
{
	std::vector<std::pair<size_t, bool>> hashes = getUnitigHashes(start, edges);
	for (auto node : hashes)
	{
		assert(!checked[node.first]);
		checked[node.first] = true;
	}
	assert(hashes.size() > 0);
	double totalCoverage = 0;
	for (auto pos : hashes)
	{
		totalCoverage += hashlist.coverage.get(pos.first);
	}
	totalCoverage /= (double)hashes.size();
	if (totalCoverage >= minUnitigCoverage)
	{
		for (auto pos : hashes)
		{
			kept.set(pos.first, true);
		}
	}
}

void startUnitig(UnitigGraph& result, std::pair<size_t, bool> start, const SparseEdgeContainer& edges, std::vector<bool>& belongsToUnitig, const HashList& hashlist, size_t minCoverage)
{
	std::vector<std::pair<size_t, bool>> hashes = getUnitigHashes(start, edges);
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	result.edgeOvlp.emplace_back();
	result.unitigs.back().reserve(hashes.size());
	result.unitigCoverage.back().reserve(hashes.size());
	for (auto pos : hashes)
	{
		result.unitigs.back().emplace_back(pos);
		result.unitigCoverage.back().emplace_back(hashlist.coverage.get(pos.first));
		assert(!belongsToUnitig[pos.first]);
		belongsToUnitig[pos.first] = true;
	}
}

SparseEdgeContainer getCoveredEdges(const HashList& hashlist, size_t minCoverage)
{
	SparseEdgeContainer result { hashlist.coverage.size() };
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : hashlist.getEdgeCoverages(fw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(fw, edge.first);
			result.addEdge(reverse(edge.first), reverse(fw));
		}
		for (auto edge : hashlist.getEdgeCoverages(bw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(bw, edge.first);
			result.addEdge(reverse(edge.first), reverse(bw));
		}
	}
	return result;
}

UnitigGraph getUnitigGraph(HashList& hashlist, const size_t minCoverage, const double minUnitigCoverage)
{
	{
		auto edges = getCoveredEdges(hashlist, minCoverage);
		RankBitvector kept { hashlist.size() };
		std::vector<bool> checked;
		checked.resize(hashlist.size(), false);
		for (size_t i = 0; i < hashlist.size(); i++)
		{
			std::pair<size_t, bool> fw { i, true };
			std::pair<size_t, bool> bw { i, false };
			auto fwEdges = edges[fw];
			auto bwEdges = edges[bw];
			if (bwEdges.size() != 1 || bwEdges[0] == fw)
			{
				if (!checked[i])
				{
					checkUnitigHashes(fw, edges, checked, hashlist, minUnitigCoverage, kept);
				}
				for (auto edge : bwEdges)
				{
					if (checked[edge.first]) continue;
					assert(hashlist.coverage.get(edge.first) >= minCoverage);
					checkUnitigHashes(edge, edges, checked, hashlist, minUnitigCoverage, kept);
				}
			}
			if (fwEdges.size() != 1 || fwEdges[0] == bw)
			{
				if (!checked[i])
				{
					checkUnitigHashes(bw, edges, checked, hashlist, minUnitigCoverage, kept);
				}
				for (auto edge : fwEdges)
				{
					if (checked[edge.first]) continue;
					assert(hashlist.coverage.get(edge.first) >= minCoverage);
					checkUnitigHashes(edge, edges, checked, hashlist, minUnitigCoverage, kept);
				}
			}
		}
		for (size_t i = 0; i < hashlist.size(); i++)
		{
			if (checked[i]) continue;
			std::pair<size_t, bool> fw { i, true };
			std::pair<size_t, bool> bw { i, false };
			auto fwEdges = edges[fw];
			auto bwEdges = edges[bw];
			assert(fwEdges.size() == 1);
			assert(bwEdges.size() == 1);
			checkUnitigHashes(fw, edges, checked, hashlist, minUnitigCoverage, kept);
		}
		kept.buildRanks();
		hashlist.filter(kept);
	}
	UnitigGraph result;
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(hashlist.coverage.size(), false);
	std::unordered_map<std::pair<size_t, bool>, std::pair<size_t, bool>> unitigTip;
	auto edges = getCoveredEdges(hashlist, minCoverage);
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage.get(i) < minCoverage) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		auto fwEdges = edges[fw];
		auto bwEdges = edges[bw];
		if (bwEdges.size() != 1 || bwEdges[0] == fw)
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
				assert(hashlist.coverage.get(edge.first) >= minCoverage);
				startUnitig(result, edge, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
		}
		if (fwEdges.size() != 1 || fwEdges[0] == bw)
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
				assert(hashlist.coverage.get(edge.first) >= minCoverage);
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
		if (hashlist.coverage.get(i) < minCoverage) continue;
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
		if (hashlist.coverage.get(i) < minCoverage) continue;
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
			assert(hashlist.coverage.get(fromNode.first) >= minCoverage);
			assert(hashlist.coverage.get(toNodeFw.first) >= minCoverage);
			result.edges[fromUnitig].emplace(toUnitig);
			result.edges[reverse(toUnitig)].emplace(reverse(fromUnitig));
			result.edgeCoverage(fromUnitig, toUnitig) = hashlist.getEdgeCoverage(fromNode, toNodeFw);
			result.edgeOverlap(fromUnitig, toUnitig) = 0;
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
				result.edgeOverlap(from, to) = oldgraph.edgeOverlap(fw, curr);
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
				result.edgeOverlap(from, to) = oldgraph.edgeOverlap(bw, curr);
			}
		}
	}
	return result;
}

size_t getN50(std::vector<size_t>& nodeSizes, size_t totalSize)
{
	if (nodeSizes.size() == 0)
	{
		assert(totalSize == 0);
		return 0;
	}
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

size_t getUnitigOverlap(const HashList& hashlist, const size_t kmerSize, const UnitigGraph& unitigs, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	auto kmerOverlap = unitigs.edgeOverlap(from, to);
	if (kmerOverlap == 0)
	{
		std::pair<size_t, bool> fromKmer = unitigs.unitigs[from.first].back();
		if (!from.second) fromKmer = reverse(unitigs.unitigs[from.first][0]);
		std::pair<size_t, bool> toKmer = unitigs.unitigs[to.first][0];
		if (!to.second) toKmer = reverse(unitigs.unitigs[to.first].back());
		return hashlist.getOverlap(fromKmer, toKmer);
	}
	assert(unitigs.unitigs[from.first].size() > kmerOverlap);
	size_t result = kmerSize;
	for (size_t i = 1; i < kmerOverlap; i++)
	{
		std::pair<size_t, bool> prev = unitigs.unitigs[from.first][unitigs.unitigs[from.first].size() - kmerOverlap + i - 1];
		std::pair<size_t, bool> curr = unitigs.unitigs[from.first][unitigs.unitigs[from.first].size() - kmerOverlap + i];
		if (!from.second)
		{
			prev = unitigs.unitigs[from.first][i-1];
			curr = unitigs.unitigs[from.first][i];
		}
		result += kmerSize - hashlist.getOverlap(prev, curr);
	}
	return result;
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
			size_t rleOverlap = getUnitigOverlap(hashlist, kmerSize, unitigs, fw, to);
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, fw, rleOverlap);
			file << "L\t" << (fw.first+1) << "\t" << (fw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		for (auto to : unitigs.edges[bw])
		{
			if (canon(bw, to).first == bw) stats.edges += 1;
			size_t rleOverlap = getUnitigOverlap(hashlist, kmerSize, unitigs, bw, to);
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

void forceEdgeConsistency(const UnitigGraph& unitigs, HashList& hashlist, const StringIndex& stringIndex, std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize)
{
	std::vector<std::vector<std::tuple<NodeType, size_t, bool>>> parent;
	std::vector<std::vector<size_t>> rank;
	parent.resize(unitigSequences.size());
	rank.resize(unitigSequences.size());
	for (size_t i = 0; i < unitigSequences.size(); i++)
	{
		// if (unitigSequences[i].compressedSize() >= 2 * kmerSize)
		// {
		// 	for (size_t j = 0; j < 2 * kmerSize; j++)
		// 	{
		// 		parent[i].emplace_back(i, j, true);
		// 		rank[i].emplace_back(0);
		// 	}
		// }
		// else
		// {
			for (size_t j = 0; j < unitigSequences[i].compressedSize(); j++)
			{
				parent[i].emplace_back(i, j, true);
				rank[i].emplace_back(0);
			}
		// }
	}
	for (size_t unitig = 0; unitig < unitigs.edges.size(); unitig++)
	{
		std::pair<size_t, bool> fw { unitig, true };
		std::pair<size_t, bool> bw { unitig, false };
		for (auto target : unitigs.edges[fw])
		{
			size_t overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, fw, target);
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
			size_t overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, bw, target);
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
			// assert(secondPos < 2 * kmerSize);
			// if (secondPos >= kmerSize)
			// {
			// 	size_t secondFromEnd = parent[secondUnitig].size() - 1 - secondPos;
			// 	assert(secondFromEnd < kmerSize);
			// 	secondPos = unitigSequences[secondUnitig].compressedSize() - 1 - secondFromEnd;
			// }
			size_t firstPos = j;
			// assert(firstPos < 2 * kmerSize);
			// if (firstPos >= kmerSize)
			// {
			// 	size_t firstFromEnd = parent[i].size() - 1 - firstPos;
			// 	assert(firstFromEnd < kmerSize);
			// 	firstPos = unitigSequences[i].compressedSize() - 1 - firstFromEnd;
			// }
			uint32_t seq = unitigSequences[secondUnitig].getExpanded(secondPos);
			if (std::get<2>(found))
			{
				assert(unitigSequences[i].getCompressed(firstPos) == unitigSequences[secondUnitig].getCompressed(secondPos));
			}
			else
			{
				assert(unitigSequences[i].getCompressed(firstPos) == complement(unitigSequences[secondUnitig].getCompressed(secondPos)));
				seq = stringIndex.getReverseIndex(unitigSequences[i].getCompressed(firstPos), seq);
			}
			assert(secondUnitig < unitigSequences.size());
			assert(i < unitigSequences.size());
			unitigSequences[i].setExpanded(firstPos, seq);
		}
	}
}

CompressedSequenceType getKmerSequence(const std::vector<CompressedSequenceType>& unitigSequences, const std::vector<std::tuple<size_t, bool, size_t>>& kmerSequencePosition, const std::pair<size_t, bool> kmer, const size_t kmerSize, const StringIndex& stringIndex)
{
	size_t unitig = std::get<0>(kmerSequencePosition[kmer.first]);
	bool fw = std::get<1>(kmerSequencePosition[kmer.first]) == kmer.second;
	size_t offset = std::get<2>(kmerSequencePosition[kmer.first]);
	assert(unitig != std::numeric_limits<size_t>::max());
	CompressedSequenceType result = unitigSequences[unitig].substr(offset, kmerSize);
	if (!fw)
	{
		result = result.revComp(stringIndex);
	}
	return result;
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

void filterKmersToUnitigKmers(UnitigGraph& unitigs, HashList& reads, const size_t kmerSize)
{
	RankBitvector kept { reads.coverage.size() };
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> newOverlaps;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::vector<std::pair<size_t, bool>> newUnitig;
		std::vector<size_t> newCoverage;
		newUnitig.reserve(unitigs.unitigs[i].size());
		newCoverage.reserve(unitigs.unitigs[i].size());
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
			assert(!kept.get(unitigs.unitigs[i][j].first));
			if (!skip)
			{
				kept.set(unitigs.unitigs[i][j].first, true);
				assert(kept.get(unitigs.unitigs[i][j].first));
				std::pair<size_t, bool> thisKmer = unitigs.unitigs[i][j];
				if (lastKmer.first != std::numeric_limits<size_t>::max())
				{
					assert(currentPos > lastStart);
					assert(currentPos - lastStart < kmerSize);
					newOverlaps.emplace_back(lastKmer, thisKmer, kmerSize - (currentPos - lastStart));
				}
				newUnitig.push_back(thisKmer);
				newCoverage.push_back(unitigs.unitigCoverage[i][j]);
				lastStart = currentPos;
				lastKmer = thisKmer;
			}
		}
		newUnitig.shrink_to_fit();
		newCoverage.shrink_to_fit();
		std::swap(unitigs.unitigs[i], newUnitig);
		std::swap(unitigs.unitigCoverage[i], newCoverage);
	}
	kept.buildRanks();
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			unitigs.unitigs[i][j].first = kept.getRank(unitigs.unitigs[i][j].first);
		}
	}
	reads.filter(kept);
	for (auto t : newOverlaps)
	{
		std::pair<size_t, bool> from { kept.getRank(std::get<0>(t).first), std::get<0>(t).second };
		std::pair<size_t, bool> to { kept.getRank(std::get<1>(t).first), std::get<1>(t).second };
		reads.addSequenceOverlap(from, to, std::get<2>(t));
	}
}

void verifyEdgeConsistency(const UnitigGraph& unitigs, const HashList& hashlist, const StringIndex& stringIndex, const std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	assert(unitigs.edges[from].count(to) == 1);
	assert(unitigs.edges[reverse(to)].count(reverse(from)) == 1);
	size_t overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, from, to);
	assert(overlap > 0);
	for (size_t i = 0; i < overlap; i++)
	{
		size_t fromIndex = unitigSequences[from.first].compressedSize() - overlap + i;
		if (!from.second) fromIndex = overlap - 1 - i;
		size_t toIndex = i;
		if (!to.second) toIndex = unitigSequences[to.first].compressedSize() - 1 - i;
		bool fw = from.second == to.second;
		if (fw)
		{
			assert(unitigSequences[from.first].getCompressed(fromIndex) == unitigSequences[to.first].getCompressed(toIndex));
			assert(unitigSequences[from.first].getExpanded(fromIndex) == unitigSequences[to.first].getExpanded(toIndex));
		}
		else
		{
			assert(unitigSequences[from.first].getCompressed(fromIndex) == complement(unitigSequences[to.first].getCompressed(toIndex)));
			assert(unitigSequences[from.first].getExpandedStr(fromIndex, stringIndex) == revCompRaw(unitigSequences[to.first].getExpandedStr(toIndex, stringIndex)));
		}
	}
}

void verifyEdgeConsistency(const UnitigGraph& unitigs, const HashList& hashlist, const StringIndex& stringIndex, const std::vector<CompressedSequenceType>& unitigSequences, const size_t kmerSize)
{
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		for (auto edge : unitigs.edges[std::make_pair(i, true)])
		{
			verifyEdgeConsistency(unitigs, hashlist, stringIndex, unitigSequences, kmerSize, std::make_pair(i, true), edge);
		}
		for (auto edge : unitigs.edges[std::make_pair(i, false)])
		{
			verifyEdgeConsistency(unitigs, hashlist, stringIndex, unitigSequences, kmerSize, std::make_pair(i, false), edge);
		}
	}
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const ErrorMasking errorMasking, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths)
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
	auto paths = loadReadsAsHashesMultithread(reads, inputReads, kmerSize, partIterator, numThreads);
	auto beforeUnitigs = getTime();
	std::cerr << "Unitigifying" << std::endl;
	auto unitigs = getUnitigGraph(reads, minCoverage, minUnitigCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage)
	{
		std::cerr << "Filtering by unitig coverage" << std::endl;
		unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	}
	filterKmersToUnitigKmers(unitigs, reads, kmerSize);
	printUnitigKmerCount(unitigs);
	auto beforeResolve = getTime();
	std::cerr << "Resolving unitigs" << std::endl;
	std::cerr << unitigs.unitigs.size() << " unitigs before resolving" << std::endl;
	std::vector<ReadPath> readPaths;
	std::tie(unitigs, readPaths) = resolveUnitigs(unitigs, reads, paths, minUnitigCoverage, kmerSize);
	std::cerr << unitigs.unitigs.size() << " unitigs after resolving" << std::endl;
	auto beforeSequences = getTime();
	std::cerr << "Getting unitig sequences" << std::endl;
	std::vector<CompressedSequenceType> unitigSequences;
	StringIndex stringIndex;
	std::tie(unitigSequences, stringIndex) = getHPCUnitigSequences(reads, unitigs, inputReads, readPaths, kmerSize, partIterator, numThreads);
	assert(unitigSequences.size() == unitigs.unitigs.size());
	auto beforeConsistency = getTime();
	AssemblyStats stats;
	beforeConsistency = getTime();
	if (errorMasking != ErrorMasking::No && errorMasking != ErrorMasking::Collapse)
	{
		std::cerr << "Forcing edge consistency" << std::endl;
		forceEdgeConsistency(unitigs, reads, stringIndex, unitigSequences, kmerSize);
	}
	verifyEdgeConsistency(unitigs, reads, stringIndex, unitigSequences, kmerSize);
	auto beforeWrite = getTime();
	std::cerr << "Writing graph to " << outputGraph << std::endl;
	stats = writeGraph(unitigs, outputGraph, reads, unitigSequences, stringIndex, kmerSize);
	auto afterWrite = getTime();
	if (outputSequencePaths != "")
	{
		std::cerr << "Writing paths to " << outputSequencePaths << std::endl;
		writePaths(reads, unitigs, unitigSequences, stringIndex, inputReads, kmerSize, partIterator, outputSequencePaths, numThreads);
	}
	auto afterPaths = getTime();
	std::cerr << "selecting k-mers and building graph topology took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeResolve) << std::endl;
	std::cerr << "resolving unitigs took " << formatTime(beforeResolve, beforeSequences) << std::endl;
	std::cerr << "building unitig sequences took " << formatTime(beforeSequences, beforeConsistency) << std::endl;
	if (errorMasking != ErrorMasking::No && errorMasking != ErrorMasking::Collapse) std::cerr << "forcing edge consistency took " << formatTime(beforeConsistency, beforeWrite) << std::endl;
	std::cerr << "writing the graph and calculating stats took " << formatTime(beforeWrite, afterWrite) << std::endl;
	if (outputSequencePaths != "") std::cerr << "writing sequence paths took " << formatTime(afterWrite, afterPaths) << std::endl;
	std::cerr << "nodes: " << stats.nodes << std::endl;
	std::cerr << "edges: " << stats.edges << std::endl;
	std::cerr << "assembly size " << stats.size << " bp, N50 " << stats.N50 << std::endl;
	if (stats.approxKmers > 0) std::cerr << "approximate number of k-mers ~ " << stats.approxKmers << std::endl;
}
