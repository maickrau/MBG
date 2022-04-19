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
#include "DumbSelect.h"

struct AssemblyStats
{
public:
	size_t nodes;
	size_t edges;
	size_t size;
	size_t N50;
	size_t approxKmers;
};

void loadReadsAsHashesMultithread(HashList& result, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	std::atomic<size_t> totalNodes = 0;
	partIterator.iterateHashes([&result, &totalNodes, kmerSize](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
		std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
		assert(positions.size() == hashes.size());
		for (size_t i = 0; i < positions.size(); i++)
		{
			const auto pos = positions[i];
			const HashType fwHash = hashes[i];
			assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
			std::pair<size_t, bool> current = result.addNode(fwHash);
			size_t overlap = lastMinimizerPosition + kmerSize - pos;
			assert(pos+kmerSize <= poses.size());
			assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition < kmerSize);
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
	std::cerr << totalNodes << " total selected k-mers in reads" << std::endl;
	std::cerr << result.size() << " distinct selected k-mers in reads" << std::endl;
}

void updatePathRemaining(size_t& rleRemaining, size_t& expanded, bool fw, const DumbSelect& expandedPoses, size_t overlap)
{
	if (rleRemaining == std::numeric_limits<size_t>::max()) return;
	size_t zeroIndex = overlap;
	if (!fw) zeroIndex = expandedPoses.countOnes() - 1 - overlap;
	if (rleRemaining < expandedPoses.countOnes() - overlap)
	{
		size_t checkIndex = overlap + rleRemaining;
		if (fw)
		{
			assert(checkIndex >= zeroIndex);
			expanded += expandedPoses.selectOne(checkIndex) - expandedPoses.selectOne(zeroIndex);
		}
		else
		{
			checkIndex = expandedPoses.countOnes() - 1 - checkIndex;
			assert(checkIndex <= zeroIndex);
			expanded += expandedPoses.selectOne(zeroIndex) - expandedPoses.selectOne(checkIndex);
		}
		rleRemaining = std::numeric_limits<size_t>::max();
	}
	else
	{
		assert(expandedPoses.countOnes() >= overlap + 2);
		rleRemaining -= expandedPoses.countOnes() - 1 - overlap;
		if (fw)
		{
			expanded += expandedPoses.selectOne(expandedPoses.countOnes()-1) - expandedPoses.selectOne(overlap);
		}
		else
		{
			expanded += expandedPoses.selectOne(expandedPoses.countOnes()-1-overlap);
		}
		assert(rleRemaining > 0);
	}
}

void writePaths(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const std::vector<ReadPath>& readPaths, const size_t kmerSize, const std::string& outputSequencePaths, const std::string& nodeNamePrefix)
{
	std::ofstream outPaths { outputSequencePaths };
	std::vector<DumbSelect> unitigExpandedPoses;
	unitigExpandedPoses.reserve(unitigSequences.size());
	for (size_t i = 0; i < unitigSequences.size(); ++i)
	{
		size_t sizeHere = 0;
		for (size_t j = 0; j < unitigSequences[i].compressedSize(); j++)
		{
			sizeHere += unitigSequences[i].getExpandedStr(j, stringIndex).size();
		}
		unitigExpandedPoses.emplace_back(sizeHere+1);
		size_t posNow = 0;
		for (size_t j = 0; j < unitigSequences[i].compressedSize(); j++)
		{
			posNow += unitigSequences[i].getExpandedStr(j, stringIndex).size();
			unitigExpandedPoses.back().bitvector.set(posNow, true);
		}
		assert(posNow == sizeHere);
		unitigExpandedPoses.back().build();
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
		for (size_t i = 0; i < path.leftClip; i++)
		{
			pathLeftClipRLE += kmerSize - hashlist.getOverlap(kmerPath[i], kmerPath[i+1]);
		}
		size_t pathRightClipRLE = 0;
		for (size_t i = 0; i < path.rightClip; i++)
		{
			pathRightClipRLE += kmerSize - hashlist.getOverlap(kmerPath[kmerPath.size()-2-i], kmerPath[kmerPath.size()-1-i]);
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
			pathStr += nodeNamePrefix;
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
			result.edges[std::make_pair(currentUnitig, belongsToUnitig.at(pos.first).second == pos.second)].emplace(currentUnitig, !(belongsToUnitig.at(pos.first).second == pos.second));
			result.edgeCoverage(currentUnitig, belongsToUnitig.at(pos.first).second == pos.second, currentUnitig, !(belongsToUnitig.at(pos.first).second == pos.second)) = old.edgeCoverage(pos.first, pos.second, newPos.first, newPos.second);
			result.edgeOverlap(currentUnitig, belongsToUnitig.at(pos.first).second == pos.second, currentUnitig, !(belongsToUnitig.at(pos.first).second == pos.second)) = old.edgeOverlap(pos.first, pos.second, newPos.first, newPos.second);
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


std::unordered_set<std::pair<size_t, bool>> findReachableNewTips(const RankBitvector& kept, const SparseEdgeContainer& edges, const HashList& hashlist, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start)
{
	assert(kept.get(start.first));
	std::vector<std::pair<size_t, bool>> stack;
	stack.push_back(start);
	std::unordered_set<std::pair<size_t, bool>> visited;
	std::unordered_set<std::pair<size_t, bool>> result;
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		if (kept.get(top.first) && top != start)
		{
			if (newlyTip[reverse(top)])
			{
				result.insert(top);
			}
			continue;
		}
		visited.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	return result;
}

void keepReachableNewTips(const RankBitvector& kept, const SparseEdgeContainer& edges, const HashList& hashlist, const VectorWithDirection<bool>& newlyTip, const std::pair<size_t, bool> start, const std::unordered_set<std::pair<size_t, bool>>& reachableTips, std::unordered_set<size_t>& newlyKept)
{
	std::unordered_set<std::pair<size_t, bool>> reachableFw;
	std::unordered_set<std::pair<size_t, bool>> reachableBw;
	std::vector<std::pair<size_t, bool>> stack;
	stack.push_back(start);
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (reachableFw.count(top) == 1) continue;
		if (kept.get(top.first) && top != start) continue;
		reachableFw.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	for (auto node : reachableTips)
	{
		stack.push_back(reverse(node));
	}
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (reachableBw.count(top) == 1) continue;
		if (kept.get(top.first) && reachableTips.count(reverse(top)) == 0) continue;
		reachableBw.insert(top);
		for (auto edge : edges[top])
		{
			stack.push_back(edge);
		}
	}
	for (auto node : reachableFw)
	{
		if (reachableBw.count(reverse(node)) == 1)
		{
			newlyKept.insert(node.first);
		}
	}
}

bool isTipGap(const RankBitvector& kept, const SparseEdgeContainer& edges, const HashList& hashlist, const std::pair<size_t, bool> start)
{
	if (!kept.get(start.first)) return false;
	if (edges[start].size() == 0) return false;
	for (auto edge : edges[start])
	{
		if (kept.get(edge.first)) return false;
	}
	return true;
}

void keepTipGaps(RankBitvector& kept, const SparseEdgeContainer& edges, const HashList& hashlist)
{
	VectorWithDirection<bool> newlyTip;
	newlyTip.resize(hashlist.size(), false);
	size_t newlyTipped = 0;
	for (size_t i = 0; i < hashlist.size(); i++)
	{
		if (!kept.get(i)) continue;
		std::pair<size_t, bool> fw { i, true };
		if (isTipGap(kept, edges, hashlist, fw))
		{
			newlyTip[fw] = true;
			newlyTipped += 1;
		}
		std::pair<size_t, bool> bw { i, false };
		if (isTipGap(kept, edges, hashlist, bw))
		{
			newlyTip[bw] = true;
			newlyTipped += 1;
		}
	}
	std::unordered_set<size_t> newlyKept;
	size_t keptCheck = 0;
	for (size_t i = 0; i < hashlist.size(); i++)
	{
		if (!kept.get(i)) continue;
		std::pair<size_t, bool> fw { i, true };
		if (edges[fw].size() > 0 && newlyTip[fw])
		{
			auto reachableTips = findReachableNewTips(kept, edges, hashlist, newlyTip, fw);
			if (reachableTips.size() > 0)
			{
				keptCheck += 1;
				keepReachableNewTips(kept, edges, hashlist, newlyTip, fw, reachableTips, newlyKept);
			}
		}
		std::pair<size_t, bool> bw { i, false };
		if (edges[bw].size() > 0 && newlyTip[bw])
		{
			auto reachableTips = findReachableNewTips(kept, edges, hashlist, newlyTip, bw);
			if (reachableTips.size() > 0)
			{
				keptCheck += 1;
				keepReachableNewTips(kept, edges, hashlist, newlyTip, bw, reachableTips, newlyKept);
			}
		}
	}
	for (auto node : newlyKept)
	{
		assert(!kept.get(node));
		kept.set(node, true);
	}

}

UnitigGraph getUnitigGraph(HashList& hashlist, const size_t minCoverage, const double minUnitigCoverage, const bool keepGaps)
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
		if (keepGaps) keepTipGaps(kept, edges, hashlist);
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

AssemblyStats writeGraph(const BluntGraph& graph, const std::string& filename, const std::string& nodeNamePrefix)
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
		file << "S\t" << nodeNamePrefix << (i+1) << "\t" << graph.nodes[i] << "\tll:f:" << graph.nodeAvgCoverage[i] << "\tFC:f:" << (graph.nodes[i].size() * graph.nodeAvgCoverage[i]) << std::endl;
		stats.size += graph.nodes[i].size();
		nodeSizes.push_back(graph.nodes[i].size());
	}
	for (auto edge : graph.edges)
	{
		file << "L\t" << nodeNamePrefix << (std::get<0>(edge)+1) << "\t" << (std::get<1>(edge) ? "+" : "-") << "\t" << nodeNamePrefix << (std::get<2>(edge)+1) << "\t" << (std::get<3>(edge) ? "+" : "-") << "\t0M\tec:i:" << std::get<4>(edge) << std::endl;
	}
	stats.N50 = getN50(nodeSizes, stats.size);
	return stats;
}

AssemblyStats writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist, const std::vector<CompressedSequenceType>& unitigSequences, const StringIndex& stringIndex, const size_t kmerSize, const std::vector<double>& unitigRawKmerCoverages, const std::string& nodeNamePrefix)
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
		file << "S\t" << nodeNamePrefix << (i+1) << "\t";
		std::string realSequence = unitigSequences[i].getExpandedSequence(stringIndex);
		file << realSequence;
		file << "\tll:f:" << unitigs.averageCoverage(i);
		file << "\tFC:f:" << (unitigs.averageCoverage(i) * realSequence.size());
		if (unitigRawKmerCoverages[i] != -1)
		{
			file << "\tkl:f:" << unitigRawKmerCoverages[i];
		}
		file << std::endl;
		stats.size += realSequence.size();
		nodeSizes.push_back(realSequence.size());
	}
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		std::vector<std::pair<size_t, bool>> fwEdges { unitigs.edges[fw].begin(), unitigs.edges[fw].end() };
		std::sort(fwEdges.begin(), fwEdges.end());
		for (auto to : fwEdges)
		{
			if (canon(fw, to).first == fw) stats.edges += 1;
			size_t rleOverlap = getUnitigOverlap(hashlist, kmerSize, unitigs, fw, to);
			size_t fromClip = fw.second ? unitigs.rightClip[fw.first] : unitigs.leftClip[fw.first];
			size_t toClip = to.second ? unitigs.leftClip[to.first] : unitigs.rightClip[to.first];
			assert(rleOverlap > fromClip + toClip);
			rleOverlap -= fromClip + toClip;
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, fw, rleOverlap);
			file << "L\t" << nodeNamePrefix << (fw.first+1) << "\t" << (fw.second ? "+" : "-") << "\t" << nodeNamePrefix << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		std::vector<std::pair<size_t, bool>> bwEdges { unitigs.edges[bw].begin(), unitigs.edges[bw].end() };
		std::sort(bwEdges.begin(), bwEdges.end());
		for (auto to : bwEdges)
		{
			if (canon(bw, to).first == bw) stats.edges += 1;
			size_t rleOverlap = getUnitigOverlap(hashlist, kmerSize, unitigs, bw, to);
			size_t fromClip = bw.second ? unitigs.rightClip[bw.first] : unitigs.leftClip[bw.first];
			size_t toClip = to.second ? unitigs.leftClip[to.first] : unitigs.rightClip[to.first];
			assert(rleOverlap > fromClip + toClip);
			rleOverlap -= fromClip + toClip;
			size_t overlap = getOverlapFromRLE(unitigSequences, stringIndex, bw, rleOverlap);
			file << "L\t" << nodeNamePrefix << (bw.first+1) << "\t" << (bw.second ? "+" : "-") << "\t" << nodeNamePrefix << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(bw, to) << std::endl;
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

void sortKmersByHashes(UnitigGraph& unitigs, HashList& reads)
{
	std::vector<size_t> kmerMapping = reads.sortByHash();
	unitigs.sort(kmerMapping);
}

void filterKmersToUnitigKmers(UnitigGraph& unitigs, HashList& reads, const size_t kmerSize)
{
	std::vector<bool> reverseContig;
	reverseContig.resize(unitigs.unitigs.size(), false);
	{
		std::unordered_map<size_t, std::pair<size_t, bool>> kmerPosition;
		std::vector<HashType> firstHash;
		std::vector<HashType> lastHash;
		for (size_t i = 0; i < unitigs.unitigs.size(); i++)
		{
			kmerPosition[unitigs.unitigs[i][0].first] = std::make_pair(i, false);
			kmerPosition[unitigs.unitigs[i].back().first] = std::make_pair(i, true);
		}
		firstHash.resize(unitigs.unitigs.size());
		lastHash.resize(unitigs.unitigs.size());
		for (auto pair : reads.hashToNode)
		{
			if (kmerPosition.count(pair.second) == 0) continue;
			auto pos = kmerPosition.at(pair.second);
			if (pos.second)
			{
				lastHash[pos.first] = pair.first;
			}
			else
			{
				firstHash[pos.first] = pair.first;
			}
		}
		for (size_t i = 0; i < unitigs.unitigs.size(); i++)
		{
			reverseContig[i] = lastHash[i] < firstHash[i];
		}
	}
	RankBitvector kept { reads.coverage.size() };
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> newOverlaps;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		if (reverseContig[i])
		{
			std::reverse(unitigs.unitigs[i].begin(), unitigs.unitigs[i].end());
			std::reverse(unitigs.unitigCoverage[i].begin(), unitigs.unitigCoverage[i].end());
			for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
			{
				unitigs.unitigs[i][j] = reverse(unitigs.unitigs[i][j]);
			}
		}
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
		if (reverseContig[i])
		{
			std::reverse(newUnitig.begin(), newUnitig.end());
			std::reverse(newCoverage.begin(), newCoverage.end());
			for (size_t j = 0; j < newUnitig.size(); j++)
			{
				newUnitig[j] = reverse(newUnitig[j]);
			}
		}
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

std::vector<ReadPath> getReadPaths(const UnitigGraph& graph, const HashList& hashlist, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize)
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
	partIterator.iterateHashes([&result, &resultMutex, &kmerLocator, kmerSize, &graph, &hashlist](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		ReadPath current;
		current.readName = read.readName;
		current.readLength = rawSeq.size();
		current.readLengthHPC = seq.size();
		size_t lastReadPos = std::numeric_limits<size_t>::max();
		std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
		std::tuple<size_t, size_t, bool> lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
		assert(positions.size() == hashes.size());
		for (size_t i = 0; i < positions.size(); i++)
		{
			const size_t readPos = positions[i];
			const HashType fwHash = hashes[i];
			assert(readPos + kmerSize <= seq.size());
			std::pair<size_t, bool> kmer = hashlist.getNodeOrNull(fwHash);
			if (kmer.first == std::numeric_limits<size_t>::max()) continue;
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
				current.readName = read.readName;
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
					current.readName = read.readName;
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
			current.readName = read.readName;
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

	return result;
}

std::vector<double> getRawKmerCoverages(const UnitigGraph& unitigs, const std::vector<CompressedSequenceType>& unitigSequences, const HashList& reads, const size_t kmerSize)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> kmerIntervals;
	kmerIntervals.resize(reads.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t unitigLength = unitigSequences[i].compressedSize() + unitigs.leftClip[i] + unitigs.rightClip[i];
		size_t currentPos = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0) currentPos -= reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			size_t kmerStart = 0;
			if (unitigs.leftClip[i] > currentPos)
			{
				kmerStart = unitigs.leftClip[i] - currentPos;
			}
			size_t kmerEnd = kmerSize;
			assert(currentPos < unitigLength - unitigs.rightClip[i]);
			if (currentPos + kmerSize > unitigLength - unitigs.rightClip[i])
			{
				kmerEnd = (unitigLength - unitigs.rightClip[i]) - currentPos;
			}
			currentPos += kmerSize;
			if (kmerEnd <= kmerStart) continue;
			if (!unitigs.unitigs[i][j].second)
			{
				std::swap(kmerStart, kmerEnd);
				kmerStart = kmerSize - kmerStart;
				kmerEnd = kmerSize - kmerEnd;
				assert(kmerEnd > kmerStart);
			}
			kmerIntervals[unitigs.unitigs[i][j].first].emplace_back(kmerStart, kmerEnd);
		}
		assert(currentPos == unitigLength);
	}
	std::vector<std::vector<std::pair<size_t, size_t>>> forbiddenIntervals;
	forbiddenIntervals.resize(reads.size());
	for (size_t i = 0; i < kmerIntervals.size(); i++)
	{
		if (kmerIntervals[i].size() == 0) continue;
		std::sort(kmerIntervals[i].begin(), kmerIntervals[i].end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return left.first < right.first; });
		size_t intervalStart = kmerIntervals[i][0].first;
		size_t forbiddenStart = kmerIntervals[i][0].second;
		size_t forbiddenEnd = kmerIntervals[i][0].second;
		size_t intervalEnd = kmerIntervals[i][0].second;
		for (size_t j = 1; j < kmerIntervals[i].size(); j++)
		{
			if (kmerIntervals[i][j].first > intervalEnd)
			{
				if (forbiddenEnd > forbiddenStart) forbiddenIntervals[i].emplace_back(forbiddenStart, forbiddenEnd);
				intervalStart = kmerIntervals[i][j].first;
				forbiddenStart = kmerIntervals[i][j].second;
				forbiddenEnd = kmerIntervals[i][j].second;
				intervalEnd = kmerIntervals[i][j].second;
			}
			else
			{
				assert(kmerIntervals[i][j].first >= intervalStart);
				assert(kmerIntervals[i][j].second > kmerIntervals[i][j].first);
				forbiddenStart = std::min(forbiddenStart, kmerIntervals[i][j].first);
				assert(forbiddenStart >= intervalStart);
				forbiddenEnd = std::max(forbiddenEnd, std::min(kmerIntervals[i][j].second, intervalEnd));
				intervalEnd = std::max(intervalEnd, kmerIntervals[i][j].second);
				assert(intervalEnd >= forbiddenEnd);
			}
			assert(forbiddenStart >= intervalStart);
			assert(forbiddenEnd >= forbiddenStart);
			assert(intervalEnd >= forbiddenEnd);
		}
		if (forbiddenEnd > forbiddenStart) forbiddenIntervals[i].emplace_back(forbiddenStart, forbiddenEnd);
	}
	for (size_t i = 0; i < forbiddenIntervals.size(); i++)
	{
		for (size_t j = 1; j < forbiddenIntervals[i].size(); j++)
		{
			assert(forbiddenIntervals[i][j].second > forbiddenIntervals[i][j].first);
			assert(forbiddenIntervals[i][j].first > forbiddenIntervals[i][j-1].second);
		}
	}
	std::vector<double> result;
	result.resize(unitigs.unitigs.size(), -1);
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t coverageCount = 0;
		size_t coverageSum = 0;
		size_t unitigLength = unitigSequences[i].compressedSize() + unitigs.leftClip[i] + unitigs.rightClip[i];
		size_t currentPos = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0) currentPos -= reads.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			size_t kmerStart = 0;
			if (unitigs.leftClip[i] > currentPos)
			{
				kmerStart = unitigs.leftClip[i] - currentPos;
			}
			size_t kmerEnd = kmerSize;
			assert(currentPos < unitigLength - unitigs.rightClip[i]);
			if (currentPos + kmerSize > unitigLength - unitigs.rightClip[i])
			{
				kmerEnd = (unitigLength - unitigs.rightClip[i]) - currentPos;
			}
			currentPos += kmerSize;
			if (kmerEnd <= kmerStart) continue;
			if (!unitigs.unitigs[i][j].second)
			{
				std::swap(kmerStart, kmerEnd);
				kmerStart = kmerSize - kmerStart;
				kmerEnd = kmerSize - kmerEnd;
				assert(kmerEnd > kmerStart);
			}
			for (auto forbidden : forbiddenIntervals[unitigs.unitigs[i][j].first])
			{
				if (forbidden.first <= kmerStart && forbidden.second >= kmerEnd)
				{
					kmerStart = 0;
					kmerEnd = 0;
					break;
				}
				assert(forbidden.first == kmerStart || forbidden.second == kmerEnd);
				if (forbidden.first == kmerStart)
				{
					assert(forbidden.second > kmerStart);
					kmerStart = forbidden.second;
					if (kmerEnd <= kmerStart) break;
				}
				else
				{
					assert(forbidden.first < kmerEnd);
					kmerEnd = forbidden.first;
					if (kmerEnd <= kmerStart) break;
				}
			}
			if (kmerEnd <= kmerStart) continue;
			coverageCount += kmerEnd - kmerStart;
			coverageSum += (kmerEnd - kmerStart) * reads.coverage.get(unitigs.unitigs[i][j].first);
		}
		assert(currentPos == unitigLength);
		if (coverageCount != 0)
		{
			result[i] = (double)coverageSum / (double)coverageCount;
		}
	}
	return result;
}

void sortPaths(std::vector<ReadPath>& readPaths)
{
	std::stable_sort(readPaths.begin(), readPaths.end(), [](const ReadPath& left, const ReadPath& right)
	{
		if (&left == &right) return false;
		if (left.readName < right.readName) return true;
		if (left.readName > right.readName) return false;
		assert(left.readName == right.readName);
		assert(left.expandedReadPosStart < left.readLength);
		assert(right.expandedReadPosStart < left.readLength);
		if (left.expandedReadPosStart < right.expandedReadPosStart) return true;
		if (left.expandedReadPosStart > right.expandedReadPosStart) return false;
		assert(false);
	});
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const ErrorMasking errorMasking, const size_t numThreads, const bool includeEndKmers, const std::string& outputSequencePaths, const size_t maxResolveLength, const bool blunt, const size_t maxUnconditionalResolveLength, const std::string& nodeNamePrefix, const std::string& sequenceCacheFile, const bool keepGaps, const double hpcVariantOnecopyCoverage, const bool guesswork)
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
	ReadpartIterator partIterator { kmerSize, windowSize, errorMasking, numThreads, inputReads, includeEndKmers, sequenceCacheFile };
	HashList reads { kmerSize };
	auto beforeVariants = getTime();
	if (hpcVariantOnecopyCoverage != 0)
	{
		std::cerr << "Collecting hpc variant k-mers" << std::endl;
		loadReadsAsHashesMultithread(reads, kmerSize, partIterator, numThreads);
		auto unitigs = getUnitigGraph(reads, minCoverage, minUnitigCoverage, keepGaps);
		if (minUnitigCoverage > minCoverage)
		{
			unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage, keepGaps));
		}
		filterKmersToUnitigKmers(unitigs, reads, kmerSize);
		sortKmersByHashes(unitigs, reads);
		auto readPaths = getReadPaths(unitigs, reads, numThreads, partIterator, kmerSize);
		getHpcVariants(reads, unitigs, readPaths, kmerSize, partIterator, numThreads, hpcVariantOnecopyCoverage * 1.5, hpcVariantOnecopyCoverage * 0.5);
		reads.clear();
		partIterator.clearCache();
	}
	std::cerr << "Collecting selected k-mers" << std::endl;
	loadReadsAsHashesMultithread(reads, kmerSize, partIterator, numThreads);
	auto beforeUnitigs = getTime();
	std::cerr << "Unitigifying" << std::endl;
	auto unitigs = getUnitigGraph(reads, minCoverage, minUnitigCoverage, keepGaps);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage)
	{
		std::cerr << "Filtering by unitig coverage" << std::endl;
		unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage, keepGaps));
	}
	filterKmersToUnitigKmers(unitigs, reads, kmerSize);
	sortKmersByHashes(unitigs, reads);
	printUnitigKmerCount(unitigs);
	auto beforePaths = getTime();
	std::vector<ReadPath> readPaths;
	std::cerr << "Getting read paths" << std::endl;
	readPaths = getReadPaths(unitigs, reads, numThreads, partIterator, kmerSize);
	auto beforeResolve = getTime();
	if (maxResolveLength > 0)
	{
		std::cerr << "Resolving unitigs" << std::endl;
		std::cerr << unitigs.unitigs.size() << " unitigs before resolving" << std::endl;
		std::tie(unitigs, readPaths) = resolveUnitigs(unitigs, reads, readPaths, partIterator, minUnitigCoverage, kmerSize, maxResolveLength, maxUnconditionalResolveLength, keepGaps, guesswork);
		std::cerr << unitigs.unitigs.size() << " unitigs after resolving" << std::endl;
	}
	auto beforeSequences = getTime();
	std::cerr << "Building unitig sequences" << std::endl;
	std::vector<CompressedSequenceType> unitigSequences;
	StringIndex stringIndex;
	std::tie(unitigSequences, stringIndex) = getHPCUnitigSequences(reads, unitigs, readPaths, kmerSize, partIterator, numThreads);
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
		stats = writeGraph(bluntGraph, outputGraph, nodeNamePrefix);
	}
	else
	{
		std::vector<double> unitigRawKmerCoverages = getRawKmerCoverages(unitigs, unitigSequences, reads, kmerSize);
		std::cerr << "Writing graph to " << outputGraph << std::endl;
		stats = writeGraph(unitigs, outputGraph, reads, unitigSequences, stringIndex, kmerSize, unitigRawKmerCoverages, nodeNamePrefix);
	}
	auto afterWrite = getTime();
	if (outputSequencePaths != "")
	{
		std::cerr << "Writing paths to " << outputSequencePaths << std::endl;
		sortPaths(readPaths);
		writePaths(reads, unitigs, unitigSequences, stringIndex, readPaths, kmerSize, outputSequencePaths, nodeNamePrefix);
	}
	auto afterPaths = getTime();
	if (hpcVariantOnecopyCoverage != 0) std::cerr << "selecting hpc variant k-mers took " << formatTime(beforeVariants, beforeReading) << std::endl;
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
