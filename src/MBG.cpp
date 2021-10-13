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
#include "UnitigHelper.h"

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
			if (seq.size() < smerSize) return;
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
				assert(pos+kmerSize < poses.size());
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

void updatePathRemaining(size_t& rleRemaining, size_t& expanded, bool fw, const std::vector<size_t>& expandedPoses, size_t overlap)
{
	if (rleRemaining == std::numeric_limits<size_t>::max()) return;
	size_t zeroIndex = overlap;
	if (!fw) zeroIndex = expandedPoses.size() - 1 - overlap;
	if (rleRemaining < expandedPoses.size() - overlap)
	{
		size_t checkIndex = overlap + rleRemaining;
		if (fw)
		{
			assert(checkIndex >= zeroIndex);
			expanded += expandedPoses[checkIndex] - expandedPoses[zeroIndex];
		}
		else
		{
			checkIndex = expandedPoses.size() - 1 - checkIndex;
			assert(checkIndex <= zeroIndex);
			expanded += expandedPoses[zeroIndex] - expandedPoses[checkIndex];
		}
		rleRemaining = std::numeric_limits<size_t>::max();
	}
	else
	{
		assert(expandedPoses.size() >= overlap + 2);
		rleRemaining -= expandedPoses.size() - 1 - overlap;
		if (fw)
		{
			expanded += expandedPoses.back() - expandedPoses[overlap];
		}
		else
		{
			expanded += expandedPoses[expandedPoses.size()-1-overlap];
		}
		assert(rleRemaining > 0);
	}
}

void writePaths(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const std::vector<ReadPath>& readPaths, const size_t kmerSize, const std::string& outputSequencePaths)
{
	std::ofstream outPaths { outputSequencePaths };
	std::vector<std::vector<size_t>> unitigExpandedPoses;
	unitigExpandedPoses.resize(unitigSequences.size());
	for (size_t i = 0; i < unitigExpandedPoses.size(); ++i)
	{
		unitigExpandedPoses[i].push_back(0);
		for (size_t j = 0; j < unitigSequences[i].compressedSize(); j++)
		{
			unitigExpandedPoses[i].push_back(unitigExpandedPoses[i].back() + unitigSequences[i].getExpandedStr(j, stringIndex).size());
		}
	}
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		size_t readLength = path.readLength;
		size_t readStart = path.expandedReadPosStart;
		size_t readEnd = path.expandedReadPosEnd;
		assert(readEnd > readStart);
		assert(readEnd <= readLength);
		std::string pathStr;
		std::vector<std::pair<size_t, bool>> kmerPath;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			size_t skip = 0;
			if (i > 0) skip = unitigs.edgeOverlap(path.path[i-1], path.path[i]);
			for (size_t j = skip; j < unitigs.unitigs[path.path[i].first].size(); j++)
			{
				if (!path.path[i].second)
				{
					kmerPath.push_back(reverse(unitigs.unitigs[path.path[i].first][unitigs.unitigs[path.path[i].first].size() - 1 - j]));
				}
				else
				{
					kmerPath.push_back(unitigs.unitigs[path.path[i].first][j]);
				}
			}
		}
		size_t pathLengthRLE = 0;
		for (size_t i = 0; i < kmerPath.size(); i++)
		{
			pathLengthRLE += kmerSize;
			if (i > 0) pathLengthRLE -= hashlist.getOverlap(kmerPath[i-1], kmerPath[i]);
		}
		size_t pathLeftClipRLE = 0;
		for (size_t i = 1; i < path.leftClip; i++)
		{
			pathLeftClipRLE += kmerSize - hashlist.getOverlap(kmerPath[i-1], kmerPath[i]);
		}
		size_t pathRightClipRLE = 0;
		for (size_t i = 1; i < path.rightClip; i++)
		{
			pathRightClipRLE += kmerSize - hashlist.getOverlap(kmerPath[kmerPath.size()-1-i], kmerPath[kmerPath.size()-i]);
		}
		size_t pathUnitigLeftClipRLE = path.path[0].second ? unitigs.leftClip[path.path[0].first] : unitigs.rightClip[path.path[0].first];
		size_t pathUnitigRightClipRLE = path.path.back().second ? unitigs.rightClip[path.path.back().first] : unitigs.leftClip[path.path.back().first];
		assert(pathLengthRLE > pathUnitigRightClipRLE + pathUnitigLeftClipRLE);
		pathLengthRLE -= pathUnitigLeftClipRLE + pathUnitigRightClipRLE;
		size_t readLeftClip = 0;
		size_t readRightClip = 0;
		if (pathLeftClipRLE > pathUnitigLeftClipRLE)
		{
			pathLeftClipRLE -= pathUnitigLeftClipRLE;
		}
		else
		{
			readLeftClip = pathUnitigLeftClipRLE - pathLeftClipRLE;
			pathLeftClipRLE = 0;
		}
		if (pathRightClipRLE > pathUnitigRightClipRLE)
		{
			pathRightClipRLE -= pathUnitigRightClipRLE;
		}
		else
		{
			readRightClip = pathUnitigRightClipRLE - pathRightClipRLE;
			pathRightClipRLE = 0;
		}
		for (size_t i = 0; i < path.path.size(); i++)
		{
			pathStr += (path.path[i].second ? ">" : "<");
			pathStr += std::to_string(path.path[i].first+1);
		}
		assert(readEnd - readStart > readRightClip + readLeftClip);
		// todo fix: readLeftClip and readRightClip are rle, readStart and readEnd are not
		readStart += readLeftClip;
		readEnd -= readRightClip;
		assert(pathLengthRLE >= kmerSize);
		assert(pathLeftClipRLE + pathRightClipRLE < pathLengthRLE);
		size_t pathStartExpanded = 0;
		size_t pathEndExpanded = 0;
		size_t pathLengthExpanded = 0;
		size_t startRemaining = pathLeftClipRLE;
		size_t endRemaining = pathLengthRLE - pathRightClipRLE;
		size_t lenRemaining = pathLengthRLE;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			size_t overlap = 0;
			if (i > 0)
			{
				overlap = getUnitigOverlap(hashlist, kmerSize, unitigs, path.path[i-1], path.path[i]);
				size_t fromClip = path.path[i-1].second ? unitigs.rightClip[path.path[i-1].first] : unitigs.leftClip[path.path[i-1].first];
				size_t toClip = path.path[i].second ? unitigs.leftClip[path.path[i].first] : unitigs.rightClip[path.path[i].first];
				assert(overlap > fromClip + toClip);
				overlap -= fromClip + toClip;
			}
			updatePathRemaining(startRemaining, pathStartExpanded, path.path[i].second, unitigExpandedPoses[path.path[i].first], overlap);
			updatePathRemaining(endRemaining, pathEndExpanded, path.path[i].second, unitigExpandedPoses[path.path[i].first], overlap);
			updatePathRemaining(lenRemaining, pathLengthExpanded, path.path[i].second, unitigExpandedPoses[path.path[i].first], overlap);
		}
		assert(startRemaining == std::numeric_limits<size_t>::max());
		assert(endRemaining == std::numeric_limits<size_t>::max());
		assert(lenRemaining == std::numeric_limits<size_t>::max());
		assert(pathLeftClipRLE == 0 || pathStartExpanded != 0);
		assert(pathEndExpanded != 0);
		assert(pathLengthExpanded != 0);
		assert(pathEndExpanded > pathStartExpanded);
		assert(pathEndExpanded <= pathLengthExpanded);
		assert(pathLengthExpanded >= pathLengthRLE);
		size_t mapq = 60;
		outPaths << path.readName << "\t" << readLength << "\t" << readStart << "\t" << readEnd << "\t+\t" << pathStr << "\t" << pathLengthExpanded << "\t" << pathStartExpanded << "\t" << pathEndExpanded << "\t" << (pathEndExpanded - pathStartExpanded) << "\t" << (pathEndExpanded - pathStartExpanded) << "\t" << mapq << std::endl;
	}
}

void startUnitig(UnitigGraph& result, const UnitigGraph& old, std::pair<size_t, bool> start, const VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>>& edges, std::vector<std::pair<size_t, bool>>& belongsToUnitig, VectorWithDirection<bool>& unitigStart, VectorWithDirection<bool>& unitigEnd)
{
	size_t currentUnitig = result.unitigs.size();
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	result.edgeOvlp.emplace_back();
	result.leftClip.emplace_back(0);
	result.rightClip.emplace_back(0);
	assert(old.leftClip[start.first] == 0);
	assert(old.rightClip[start.first] == 0);
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
		assert(old.leftClip[pos.first] == 0);
		assert(old.rightClip[pos.first] == 0);
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
	result.leftClip.emplace_back(0);
	result.rightClip.emplace_back(0);
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
			size_t fromClip = fw.second ? unitigs.rightClip[fw.first] : unitigs.leftClip[fw.first];
			size_t toClip = to.second ? unitigs.leftClip[to.first] : unitigs.rightClip[to.first];
			assert(rleOverlap > fromClip + toClip);
			rleOverlap -= fromClip + toClip;
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, fw, rleOverlap);
			file << "L\t" << (fw.first+1) << "\t" << (fw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		for (auto to : unitigs.edges[bw])
		{
			if (canon(bw, to).first == bw) stats.edges += 1;
			size_t rleOverlap = getUnitigOverlap(hashlist, kmerSize, unitigs, bw, to);
			size_t fromClip = bw.second ? unitigs.rightClip[bw.first] : unitigs.leftClip[bw.first];
			size_t toClip = to.second ? unitigs.leftClip[to.first] : unitigs.rightClip[to.first];
			assert(rleOverlap > fromClip + toClip);
			rleOverlap -= fromClip + toClip;
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
	size_t fromClip = from.second ? unitigs.rightClip[from.first] : unitigs.leftClip[from.first];
	size_t toClip = to.second ? unitigs.leftClip[to.first] : unitigs.rightClip[to.first];
	assert(overlap > fromClip + toClip);
	overlap -= fromClip + toClip;
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

std::vector<ReadPath> getReadPaths(const UnitigGraph& graph, const HashList& hashlist, const std::vector<std::string>& readFiles, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize)
{
	size_t maxKmer = 0;
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			maxKmer = std::max(maxKmer, graph.unitigs[i][j].first);
		}
	}
	std::vector<std::tuple<size_t, size_t, bool>> kmerLocator;
	kmerLocator.resize(maxKmer+1, std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true));
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			assert(std::get<0>(kmerLocator[graph.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			if (graph.unitigs[i][j].second)
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, j, true);
			}
			else
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, graph.unitigs[i].size()-1-j, false);
			}
		}
	}
	std::vector<ReadPath> result;
	std::mutex resultMutex;
	iterateReadsMultithreaded(readFiles, numThreads, [&result, &resultMutex, &kmerLocator, kmerSize, &graph, &hashlist, &partIterator](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &resultMutex, &kmerLocator, kmerSize, &graph, &hashlist, &read](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, uint64_t minHash, const std::vector<size_t>& positions)
		{
			ReadPath current;
			current.readName = read.seq_id;
			current.readLength = rawSeq.size();
			current.readLengthHPC = seq.size();
			size_t lastReadPos = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
			std::tuple<size_t, size_t, bool> lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			for (const size_t readPos : positions)
			{
				VectorView<CharType> minimizerSequence { seq, readPos, readPos + kmerSize };
				assert(readPos + kmerSize <= seq.size());
				std::pair<size_t, bool> kmer = hashlist.getNodeOrNull(minimizerSequence);
				if (kmer.first == std::numeric_limits<size_t>::max()) continue;
				size_t readPosExpandedStart = poses[readPos];
				size_t readPosExpandedEnd = poses[readPos+kmerSize];
				if (kmer.first >= kmerLocator.size() || std::get<0>(kmerLocator[kmer.first]) == std::numeric_limits<size_t>::max())
				{
					if (current.path.size() > 0)
					{
						assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
						std::lock_guard<std::mutex> lock { resultMutex };
						result.emplace_back();
						std::swap(result.back(), current);
					}
					current.path.clear();
					current.readName = read.seq_id;
					current.readLength = rawSeq.size();
					current.readLengthHPC = seq.size();
					current.leftClip = 0;
					current.rightClip = 0;
					lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
					lastReadPos = readPos;
					lastKmer = kmer;
					continue;
				}
				auto pos = kmerLocator[kmer.first];
				assert(std::get<0>(pos) < graph.unitigs.size());
				if (!kmer.second)
				{
					std::get<2>(pos) = !std::get<2>(pos);
					std::get<1>(pos) = graph.unitigs[std::get<0>(pos)].size()-1-std::get<1>(pos);
				}
				if (std::get<0>(lastPos) == std::numeric_limits<size_t>::max())
				{
					assert(current.path.size() == 0);
					current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
					current.readPoses.push_back(readPos);
					current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
					current.leftClip = std::get<1>(pos);
					assert(current.leftClip + current.rightClip + 1 == graph.unitigs[std::get<0>(pos)].size());
					lastPos = pos;
					lastReadPos = readPos;
					lastKmer = kmer;
					continue;
				}
				if (lastKmer.first != std::numeric_limits<size_t>::max())
				{
					assert(lastReadPos != std::numeric_limits<size_t>::max());
					if (!hashlist.hasSequenceOverlap(lastKmer, kmer) || readPos - lastReadPos != kmerSize - hashlist.getOverlap(lastKmer, kmer))
					{
						assert(current.path.size() > 0);
						{
							assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
							std::lock_guard<std::mutex> lock { resultMutex };
							result.emplace_back();
							std::swap(result.back(), current);
						}
						current.path.clear();
						current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
						current.readPoses.clear();
						current.readPoses.push_back(readPos);
						current.leftClip = std::get<1>(pos);
						current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
						current.readName = read.seq_id;
						current.readLength = rawSeq.size();
						current.readLengthHPC = seq.size();
						assert(current.leftClip + current.rightClip + 1 == graph.unitigs[std::get<0>(pos)].size());
						lastPos = pos;
						lastReadPos = readPos;
						lastKmer = kmer;
						continue;
					}
				}
				assert(current.path.size() > 0);
				if (std::get<0>(pos) == std::get<0>(lastPos) && std::get<2>(pos) == std::get<2>(lastPos) && std::get<1>(pos) == std::get<1>(lastPos)+1)
				{
					lastPos = pos;
					lastReadPos = readPos;
					lastKmer = kmer;
					current.rightClip -= 1;
					current.readPoses.push_back(readPos);
					assert(current.rightClip == graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos));
					continue;
				}
				std::pair<size_t, bool> fromEdge { std::get<0>(lastPos), std::get<2>(lastPos) };
				std::pair<size_t, bool> toEdge { std::get<0>(pos), std::get<2>(pos) };
				assert(graph.edges[fromEdge].count(toEdge) == graph.edges[reverse(toEdge)].count(reverse(fromEdge)));
				if (graph.edges[fromEdge].count(toEdge) == 1 && std::get<1>(pos) == 0 && std::get<1>(lastPos) == graph.unitigs[std::get<0>(lastPos)].size()-1)
				{
					assert(current.rightClip == 0);
					current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
					current.readPoses.push_back(readPos);
					current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1;
					lastPos = pos;
					lastReadPos = readPos;
					lastKmer = kmer;
					continue;
				}
				assert(current.path.size() > 0);
				{
					assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
					std::lock_guard<std::mutex> lock { resultMutex };
					result.emplace_back();
					std::swap(result.back(), current);
				}
				current.path.clear();
				current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
				current.readPoses.clear();
				current.readPoses.push_back(readPos);
				current.leftClip = std::get<1>(pos);
				current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
				current.readName = read.seq_id;
				current.readLength = rawSeq.size();
				current.readLengthHPC = seq.size();
				assert(current.leftClip + current.rightClip + 1 == graph.unitigs[std::get<0>(pos)].size());
				lastPos = pos;
				lastReadPos = readPos;
				lastKmer = kmer;
				continue;
			}
			if (current.path.size() > 0)
			{
				assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
				std::lock_guard<std::mutex> lock { resultMutex };
				result.emplace_back();
				std::swap(result.back(), current);
			}
		});
	});

	return result;
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const ErrorMasking errorMasking, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths, const size_t maxResolveLength, const bool blunt)
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
	auto unitigs = getUnitigGraph(reads, minCoverage, minUnitigCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage)
	{
		std::cerr << "Filtering by unitig coverage" << std::endl;
		unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	}
	filterKmersToUnitigKmers(unitigs, reads, kmerSize);
	printUnitigKmerCount(unitigs);
	auto beforePaths = getTime();
	std::vector<ReadPath> readPaths;
	std::cerr << "Getting read paths" << std::endl;
	readPaths = getReadPaths(unitigs, reads, inputReads, numThreads, partIterator, kmerSize);
	auto beforeResolve = getTime();
	if (maxResolveLength > 0)
	{
		std::cerr << "Resolving unitigs" << std::endl;
		std::cerr << unitigs.unitigs.size() << " unitigs before resolving" << std::endl;
		std::tie(unitigs, readPaths) = resolveUnitigs(unitigs, reads, readPaths, partIterator, minUnitigCoverage, kmerSize, maxResolveLength);
		std::cerr << unitigs.unitigs.size() << " unitigs after resolving" << std::endl;
	}
	auto beforeSequences = getTime();
	std::cerr << "Building unitig sequences" << std::endl;
	std::vector<CompressedSequenceType> unitigSequences;
	StringIndex stringIndex;
	std::tie(unitigSequences, stringIndex) = getHPCUnitigSequences(reads, unitigs, inputReads, readPaths, kmerSize, partIterator, numThreads);
	assert(unitigSequences.size() == unitigs.unitigs.size());
	auto beforeConsistency = getTime();
	AssemblyStats stats;
	beforeConsistency = getTime();
	verifyEdgeConsistency(unitigs, reads, stringIndex, unitigSequences, kmerSize);
	auto beforeWrite = getTime();
	if (blunt)
	{
		BluntGraph bluntGraph { reads, unitigs, unitigSequences, stringIndex, kmerSize };
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(bluntGraph, outputGraph);
	}
	else
	{
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(unitigs, outputGraph, reads, unitigSequences, stringIndex, kmerSize);
	}
	auto afterWrite = getTime();
	if (outputSequencePaths != "")
	{
		std::cerr << "Writing paths to " << outputSequencePaths << std::endl;
		writePaths(reads, unitigs, unitigSequences, stringIndex, readPaths, kmerSize, outputSequencePaths);
	}
	auto afterPaths = getTime();
	std::cerr << "selecting k-mers and building graph topology took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforePaths) << std::endl;
	std::cerr << "getting read paths took " << formatTime(beforePaths, beforeResolve) << std::endl;
	if (maxResolveLength > 0) std::cerr << "resolving unitigs took " << formatTime(beforeResolve, beforeSequences) << std::endl;
	std::cerr << "building unitig sequences took " << formatTime(beforeSequences, beforeConsistency) << std::endl;
	if (errorMasking != ErrorMasking::No && errorMasking != ErrorMasking::Collapse) std::cerr << "forcing edge consistency took " << formatTime(beforeConsistency, beforeWrite) << std::endl;
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
