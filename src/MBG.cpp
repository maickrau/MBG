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
#include "fastqloader.h"
#include "CommonUtils.h"
#include "MBGCommon.h"
#include "MBG.h"
#include "VectorWithDirection.h"
#include "TwobitString.h"
#include "FastHasher.h"
#include "SparseEdgeContainer.h"
#include "HashList.h"
#include "TransitiveCleaner.h"
#include "UnitigGraph.h"
#include "BluntGraph.h"

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
				if (resultStr.size() == 0) continue;
				result.emplace_back(std::move(resultStr), std::move(lens));
				resultStr.clear();
				lens.clear();
				break;
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

std::string strFromRevComp(const std::string& revComp)
{
	std::string result;
	result.reserve(revComp.size());
	for (size_t i = 0; i < revComp.size(); i++)
	{
		result.push_back("-ACGT"[revComp[i]]);
	}
	return result;
}

template <typename F>
void findMinimizerPositions(const std::string& sequence, size_t kmerSize, size_t windowSize, F callback)
{
	if (sequence.size() < kmerSize + windowSize) return;
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	std::deque<std::tuple<size_t, uint64_t, uint64_t, uint64_t>> minimizerOrder;
	minimizerOrder.emplace_back(0, fwkmerHasher.hash(), fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	for (size_t i = 0; i < windowSize-1; i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && std::get<1>(minimizerOrder.back()) > hash) minimizerOrder.pop_back();
		minimizerOrder.emplace_back(i+1, hash, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	}
	auto pos = minimizerOrder.begin();
	while (pos != minimizerOrder.end() && std::get<1>(*pos) == std::get<1>(minimizerOrder.front()))
	{
		callback(std::get<0>(*pos), std::get<2>(*pos), std::get<3>(*pos));
		++pos;
	}
	for (size_t i = windowSize-1; kmerSize+i < sequence.size(); i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		auto oldMinimizer = std::get<1>(minimizerOrder.front());
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && std::get<0>(minimizerOrder.front()) <= i + 1 - windowSize) minimizerOrder.pop_front();
		while (minimizerOrder.size() > 0 && std::get<1>(minimizerOrder.back()) > hash) minimizerOrder.pop_back();
		if (minimizerOrder.size() > 0 && oldMinimizer != std::get<1>(minimizerOrder.front()))
		{
			auto pos = minimizerOrder.begin();
			while (pos != minimizerOrder.end() && std::get<1>(*pos) == std::get<1>(minimizerOrder.front()))
			{
				callback(std::get<0>(*pos), std::get<2>(*pos), std::get<3>(*pos));
				++pos;
			}
		}
		if (minimizerOrder.size() == 0 || hash == std::get<1>(minimizerOrder.front())) callback(i+1, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
		minimizerOrder.emplace_back(i+1, hash, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	}
}

template <typename F>
void findSyncmerPositions(const std::string& sequence, size_t kmerSize, size_t smerSize, F callback)
{
	if (sequence.size() < kmerSize) return;
	assert(smerSize < kmerSize);
	size_t windowSize = kmerSize - smerSize + 1;
	assert(windowSize >= 1);
	FastHasher fwkmerHasher { smerSize };
	for (size_t i = 0; i < smerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	std::deque<std::tuple<size_t, uint64_t>> smerOrder;
	smerOrder.emplace_back(0, fwkmerHasher.hash());
	for (size_t i = 1; i < windowSize; i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		size_t hash = fwkmerHasher.hash();
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
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i - windowSize) smerOrder.pop_front();
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
		if ((std::get<0>(smerOrder.front()) == i-windowSize+1) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i))
		{
			callback(i-windowSize+1);
		}
	}
}

void cleanTransitiveEdges(HashList& result, size_t kmerSize)
{
	TransitiveCleaner cleaner { kmerSize, result };
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> addEdgeCoverage;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> removeEdgeCoverage;

	size_t transitiveEdgesBroken = 0;
	for (size_t node = 0; node < result.edgeCoverage.size(); node++)
	{
		std::pair<size_t, bool> fw { node, true };
		for (auto target : result.edgeCoverage[fw])
		{
			std::vector<std::pair<size_t, bool>> vec;
			vec.push_back(fw);
			vec.push_back(target.first);
			vec = cleaner.insertMiddles(vec);
			if (vec.size() == 2) continue;
			transitiveEdgesBroken += 1;
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(vec[0], vec.back());
			removeEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			for (size_t i = 1; i < vec.size(); i++)
			{
				std::tie(canonFrom, canonTo) = canon(vec[i-1], vec[i]);
				addEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				result.coverage[vec[i].first] += target.second;
			}
		}
		std::pair<size_t, bool> bw { node, false };
		for (auto target : result.edgeCoverage[bw])
		{
			std::vector<std::pair<size_t, bool>> vec;
			vec.push_back(bw);
			vec.push_back(target.first);
			vec = cleaner.insertMiddles(vec);
			if (vec.size() == 2) continue;
			transitiveEdgesBroken += 1;
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(vec[0], vec.back());
			removeEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			for (size_t i = 1; i < vec.size(); i++)
			{
				std::tie(canonFrom, canonTo) = canon(vec[i-1], vec[i]);
				addEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				result.coverage[vec[i].first] += target.second;
			}
		}
	}
	for (auto t : cleaner.newSequenceOverlaps)
	{
		result.sequenceOverlap[std::get<0>(t)][std::get<1>(t)] = std::get<2>(t);
	}
	for (auto t : addEdgeCoverage)
	{
		result.edgeCoverage[std::get<0>(t)][std::get<1>(t)] += std::get<2>(t);
	}
	for (auto t : removeEdgeCoverage)
	{
		assert(result.edgeCoverage.at(std::get<0>(t)).count(std::get<1>(t)) == 1);
		assert(result.edgeCoverage.at(std::get<0>(t)).at(std::get<1>(t)) >= std::get<2>(t));
		result.edgeCoverage[std::get<0>(t)][std::get<1>(t)] -= std::get<2>(t);
	}
	std::cerr << transitiveEdgesBroken << " transitive edges cleaned" << std::endl;
}

HashList loadReadsAsHashes(const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const bool hpc, const bool collapseRunLengths)
{
	HashList result { kmerSize, collapseRunLengths };
	size_t totalNodes = 0;
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&result, &totalNodes, kmerSize, windowSize, hpc](const FastQ& read){
			if (read.sequence.size() == 0) return;
			std::vector<std::pair<std::string, std::vector<uint16_t>>> parts;
			if (hpc)
			{
				parts = runLengthEncode(read.sequence);
			}
			else
			{
				parts = noRunLengthEncode(read.sequence);
			}
			for (size_t i = 0; i < parts.size(); i++)
			{
				std::string& seq = parts[i].first;
				std::vector<uint16_t>& lens = parts[i].second;
				if (seq.size() <= kmerSize + windowSize) return;
				std::string revSeq = revCompRLE(seq);
				size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
				std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
				HashType lastHash = 0;
				findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, [kmerSize, windowSize, &lastHash, &last, &lens, &seq, &revSeq, &result, &lastMinimizerPosition, &totalNodes](size_t pos)
				{
					assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos > lastMinimizerPosition);
					assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= windowSize);
					assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
					std::string_view minimizerSequence { seq.data() + pos, kmerSize };
					size_t revPos = seq.size() - (pos + kmerSize);
					std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
					std::pair<size_t, bool> current;
					size_t overlap = lastMinimizerPosition + kmerSize - pos;
					std::tie(current, lastHash) = result.addNode(minimizerSequence, revMinimizerSequence, lens, pos, pos + kmerSize, lastHash, overlap, 0, 0);
					if (last.first != std::numeric_limits<size_t>::max() && pos - lastMinimizerPosition < kmerSize)
					{
						assert(lastMinimizerPosition + kmerSize >= pos);
						result.addSequenceOverlap(last, current, overlap);
						auto pair = canon(last, current);
						result.edgeCoverage[pair.first][pair.second] += 1;
					}
					lastMinimizerPosition = pos;
					result.coverage[current.first] += 1;
					last = current;
					totalNodes += 1;
				});
			}
		});
	}
	result.buildReverseCompHashSequences();
	std::cerr << totalNodes << " nodes" << std::endl;
	std::cerr << result.size() << " distinct fw/bw sequence nodes" << std::endl;
	return result;
}

void startUnitig(UnitigGraph& result, const UnitigGraph& old, std::pair<size_t, bool> start, const VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>>& edges, std::vector<std::pair<size_t, bool>>& belongsToUnitig)
{
	size_t currentUnitig = result.unitigs.size();
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	std::pair<size_t, bool> pos = start;
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
}

void startUnitig(UnitigGraph& result, std::pair<size_t, bool> start, const SparseEdgeContainer& edges, std::vector<bool>& belongsToUnitig, const HashList& hashlist, size_t minCoverage)
{
	size_t currentUnitig = result.unitigs.size();
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
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		auto bw = std::make_pair(node, false);
		if (edges.at(fw).size() != 1)
		{
			for (auto start : edges.at(fw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, bw, edges, belongsToUnitig);
		}
		if (edges.at(bw).size() != 1)
		{
			for (auto start : edges.at(bw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig);
		}
	}
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig);
	}
	for (size_t i = 0; i < oldgraph.edges.size(); i++)
	{
		auto fw = std::make_pair(i, true);
		for (auto curr : oldgraph.edges.at(fw))
		{
			auto from = belongsToUnitig.at(fw.first);
			bool prevFw = fw.second;
			auto to = belongsToUnitig.at(curr.first);
			bool currFw = curr.second;
			if (from.first == to.first) continue;
			from.second = !(from.second ^ prevFw);
			to.second = !(to.second ^ currFw);
			result.edges[from].emplace(to);
			result.edgeCoverage(from, to) = oldgraph.edgeCoverage(fw, curr);
		}
		auto bw = std::make_pair(i, false);
		for (auto curr : oldgraph.edges.at(bw))
		{
			auto from = belongsToUnitig.at(bw.first);
			bool prevFw = bw.second;
			auto to = belongsToUnitig.at(curr.first);
			bool currFw = curr.second;
			if (from.first == to.first) continue;
			from.second = !(from.second ^ prevFw);
			to.second = !(to.second ^ currFw);
			result.edges[from].emplace(to);
			result.edgeCoverage(from, to) = oldgraph.edgeCoverage(bw, curr);
		}
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

size_t getOverlapFromRLE(const HashList& hashlist, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	size_t overlap = hashlist.getOverlap(from, to);
	size_t RLEoverlap = 0;
	auto lens = hashlist.getHashCharacterLength(to.first);
	assert(lens.size() > overlap);
	for (size_t offset = 0; offset < overlap; offset++)
	{
		size_t i;
		if (to.second)
		{
			i = offset;
		}
		else
		{
			i = lens.size()-offset-1;
		}
		RLEoverlap += lens[i];
	}
	return RLEoverlap;
}

void writeGraph(const BluntGraph& graph, const std::string& filename)
{
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		file << "S\t" << (i+1) << "\t" << graph.nodes[i] << "\tll:f:" << graph.nodeAvgCoverage[i] << "\tFC:f:" << (graph.nodes[i].size() * graph.nodeAvgCoverage[i]) << std::endl;
	}
	for (auto edge : graph.edges)
	{
		file << "L\t" << (std::get<0>(edge)+1) << "\t" << (std::get<1>(edge) ? "+" : "-") << "\t" << (std::get<2>(edge)+1) << "\t" << (std::get<3>(edge) ? "+" : "-") << "\t0M\tec:i:" << std::get<4>(edge) << std::endl;
	}
}

void writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist)
{
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		file << "S\t" << (i+1) << "\t";
		size_t length = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			auto to = unitigs.unitigs[i][j];
			std::string sequenceRLE;
			std::vector<uint16_t> sequenceCharacterLength = hashlist.getHashCharacterLength(to.first);
			if (to.second)
			{
				sequenceRLE = hashlist.getHashSequenceRLE(to.first).toString();
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first).toString();
				std::reverse(sequenceCharacterLength.begin(), sequenceCharacterLength.end());
			}
			if (j > 0)
			{
				auto from = unitigs.unitigs[i][j-1];
				size_t overlap = hashlist.getOverlap(from, to);
				assert(overlap < sequenceRLE.size());
				sequenceRLE = sequenceRLE.substr(overlap);
				sequenceCharacterLength.erase(sequenceCharacterLength.begin(), sequenceCharacterLength.begin() + overlap);
			}
			std::string sequence = getSequence(sequenceRLE, sequenceCharacterLength);
			file << sequence;
			length += sequence.size();
		}
		file << "\tll:f:" << unitigs.averageCoverage(i);
		file << "\tFC:f:" << (unitigs.averageCoverage(i) * length);
		file << std::endl;
	}
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : unitigs.edges[fw])
		{
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
			size_t overlap = getOverlapFromRLE(hashlist, last, first);
			file << "L\t" << (fw.first+1) << "\t" << (fw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		for (auto to : unitigs.edges[bw])
		{
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
			size_t overlap = getOverlapFromRLE(hashlist, last, first);
			file << "L\t" << (bw.first+1) << "\t" << (bw.second ? "+" : "-") << "\t" << (to.first+1) << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(bw, to) << std::endl;
		}
	}
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

struct AssemblyStats
{
public:
	size_t nodes;
	size_t edges;
	size_t size;
	size_t N50;
	size_t approxKmers;
};

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

AssemblyStats getSizeAndN50(const HashList& hashlist, const UnitigGraph& graph, const size_t kmerSize, const size_t windowSize)
{
	size_t total = 0;
	std::vector<size_t> sizes;
	for (const auto& unitig : graph.unitigs)
	{
		size_t length = 0;
		for (size_t j = 0; j < unitig.size(); j++)
		{
			auto to = unitig[j];
			std::string sequenceRLE;
			std::vector<uint16_t> sequenceCharacterLength = hashlist.getHashCharacterLength(to.first);
			if (to.first)
			{
				sequenceRLE = hashlist.getHashSequenceRLE(to.first).toString();
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first).toString();
				std::reverse(sequenceCharacterLength.begin() ,sequenceCharacterLength.end());
			}
			if (j > 0)
			{
				auto from = unitig[j-1];
				size_t overlap = hashlist.getOverlap(from, to);
				assert(overlap < sequenceRLE.size());
				sequenceRLE = sequenceRLE.substr(overlap);
				sequenceCharacterLength.erase(sequenceCharacterLength.begin(), sequenceCharacterLength.begin() + overlap);
			}
			length += getSequence(sequenceRLE, sequenceCharacterLength).size();
		}
		total += length;
		sizes.push_back(length);
	}
	std::sort(sizes.begin(), sizes.end());
	size_t partialSum = 0;
	size_t N50 = 0;
	for (size_t i = sizes.size()-1; i < sizes.size(); i--)
	{
		partialSum += sizes[i];
		if (partialSum >= total * 0.5)
		{
			N50 = sizes[i];
			break;
		}
	}
	AssemblyStats result;
	result.nodes = graph.numNodes();
	result.edges = graph.numEdges();
	result.size = total;
	result.N50 = N50;
	result.approxKmers = total - graph.numNodes() * (kmerSize - windowSize/2 - 1);
	return result;
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool collapseRunLengths, const bool blunt)
{
	auto beforeReading = getTime();
	auto reads = loadReadsAsHashes(inputReads, kmerSize, windowSize, hpc, collapseRunLengths);
	auto beforeCleaning = getTime();
	// cleanTransitiveEdges(reads, kmerSize);
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage) unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	auto beforeWrite = getTime();
	AssemblyStats stats;
	if (blunt)
	{
		BluntGraph blunt { reads, unitigs };
		writeGraph(blunt, outputGraph);
		stats = getSizeAndN50(blunt);
	}
	else
	{
		writeGraph(unitigs, outputGraph, reads);
		stats = getSizeAndN50(reads, unitigs, kmerSize, windowSize);
	}
	auto afterWrite = getTime();
	std::cerr << "reading and hashing sequences took " << formatTime(beforeReading, beforeCleaning) << std::endl;
	std::cerr << "cleaning transitive edges took " << formatTime(beforeCleaning, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeWrite) << std::endl;
	std::cerr << "writing the graph and calculating stats took " << formatTime(beforeWrite, afterWrite) << std::endl;
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
