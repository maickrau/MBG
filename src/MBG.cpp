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
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
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

template <typename F>
uint64_t findSyncmerPositions(const std::string& sequence, size_t kmerSize, size_t smerSize, std::vector<std::tuple<size_t, uint64_t>>& smerOrder, F callback)
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
	smerOrder.emplace_back(0, fwkmerHasher.hash());
	for (size_t i = 1; i < windowSize; i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		size_t hash = fwkmerHasher.hash();
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
		minHash = std::min(minHash, hash);
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
		minHash = std::min(minHash, hash);
	}
	return minHash;
}

void loadReadsAsHashesMultithread(HashList& result, const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const bool hpc, const bool collapseRunLengths, const size_t numThreads)
{
	// check that all files actually exist
	for (const std::string& name : files)
	{
		std::ifstream file { name };
		if (!file.good())
		{
			std::cerr << "Input file " << name << " can't be read!" << std::endl;
			std::exit(1);
		}
	}
	std::atomic<size_t> totalNodes = 0;
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> sequenceQueue;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &result, &sequenceQueue, &totalNodes, hpc, kmerSize, windowSize, collapseRunLengths]()
		{
			// keep the same smerOrder to reduce mallocs which destroy multithreading performance
			std::vector<std::tuple<size_t, uint64_t>> smerOrder;
			smerOrder.reserve(windowSize);
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
				if (read->sequence.size() == 0) continue;
				std::vector<std::pair<std::string, std::vector<uint16_t>>> parts;
				if (hpc)
				{
					parts = runLengthEncode(read->sequence);
				}
				else
				{
					parts = noRunLengthEncode(read->sequence);
				}
				for (size_t i = 0; i < parts.size(); i++)
				{
					std::string& seq = parts[i].first;
					std::vector<uint16_t>& lens = parts[i].second;
					if (seq.size() <= kmerSize + windowSize) continue;
					std::string revSeq = revCompRLE(seq);
					size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
					std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
					HashType lastHash = 0;
					smerOrder.resize(0);
					std::vector<size_t> positions;
					uint64_t minHash = findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [&positions](size_t pos)
					{
						positions.push_back(pos);
					});
					for (auto pos : positions)
					{
						assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos > lastMinimizerPosition);
						assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= windowSize);
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
				}
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
			sequenceQueue.enqueue(ptr);
		});
	}
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	std::cerr << totalNodes << " nodes" << std::endl;
	std::cerr << result.size() << " distinct fw/bw sequence nodes" << std::endl;
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

void writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist, const size_t kmerSize)
{
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		file << "S\t" << (i+1) << "\t";
		size_t length = 0;
		size_t printBeforeLast = 0;
		for (size_t j = 1; j < unitigs.unitigs[i].size(); j++)
		{
			size_t overlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			printBeforeLast += kmerSize - overlap;
		}
		size_t prunt = 0;
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
				sequenceRLE = revCompRLE(hashlist.getHashSequenceRLE(to.first)).toString();
				std::reverse(sequenceCharacterLength.begin(), sequenceCharacterLength.end());
			}
			if (j > 0 && j < unitigs.unitigs[i].size()-1)
			{
				auto from = unitigs.unitigs[i][j-1];
				size_t overlap = hashlist.getOverlap(from, to);
				assert(overlap < sequenceRLE.size());
				sequenceRLE = sequenceRLE.substr(overlap);
				sequenceCharacterLength.erase(sequenceCharacterLength.begin(), sequenceCharacterLength.begin() + overlap);
			}
			assert(prunt <= printBeforeLast);
			if (j < unitigs.unitigs[i].size()-1 && sequenceRLE.size() > printBeforeLast - prunt)
			{
				sequenceRLE.erase(sequenceRLE.begin() + printBeforeLast - prunt, sequenceRLE.end());
				sequenceCharacterLength.erase(sequenceCharacterLength.begin() + printBeforeLast - prunt, sequenceCharacterLength.end());
			}
			prunt += sequenceRLE.size();
			std::string sequence = getSequence(sequenceRLE, sequenceCharacterLength);
			file << sequence;
			length += sequence.size();
		}
		assert(prunt == printBeforeLast + kmerSize);
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
				sequenceRLE = revCompRLE(hashlist.getHashSequenceRLE(to.first)).toString();
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

void forceEdgeConsistency(const UnitigGraph& unitigs, HashList& hashlist, const size_t kmerSize)
{
	std::unordered_map<NodeType, std::vector<std::pair<NodeType, size_t>>> parent;
	std::unordered_map<NodeType, std::vector<size_t>> rank;
	for (const auto& unitig : unitigs.unitigs)
	{
		assert(unitig.size() > 0);
		assert(parent.count(unitig[0].first) == 0);
		assert(parent.count(unitig.back().first) == 0);
		assert(rank.count(unitig[0].first) == 0);
		assert(rank.count(unitig.back().first) == 0);
		for (size_t i = 0; i < kmerSize; i++)
		{
			parent[unitig[0].first].emplace_back(unitig[0].first, i);
			rank[unitig[0].first].emplace_back(0);
			if (unitig.size() > 1)
			{
				parent[unitig.back().first].emplace_back(unitig.back().first, i);
				rank[unitig.back().first].emplace_back(0);
			}
		}
		if (unitig.size() == 1) continue;
		size_t firstToLastOverlap = kmerSize;
		for (size_t i = 1; i < unitig.size(); i++)
		{
			size_t overlap = hashlist.getOverlap(unitig[i-1], unitig[i]);
			assert(overlap < kmerSize);
			size_t removeOverlap = kmerSize - overlap;
			if (removeOverlap >= firstToLastOverlap)
			{
				firstToLastOverlap = 0;
				break;
			}
			firstToLastOverlap -= removeOverlap;
			if (i == unitig.size()-1) break;
		}
		if (firstToLastOverlap > 0)
		{
			assert(firstToLastOverlap < kmerSize);
			for (size_t i = 0; i < firstToLastOverlap; i++)
			{
				size_t firstOffset = kmerSize - firstToLastOverlap + i;
				assert(firstOffset < kmerSize);
				if (!unitig[0].second) firstOffset = kmerSize - 1 - firstOffset;
				assert(firstOffset < kmerSize);
				size_t secondOffset = i;
				if (!unitig.back().second) secondOffset = kmerSize - 1 - secondOffset;
				assert(secondOffset < kmerSize);
				merge(parent, rank, unitig[0].first, firstOffset, unitig.back().first, secondOffset);
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
				size_t firstOffset = kmerSize - overlap + i;
				assert(firstOffset < kmerSize);
				if (!from.second) firstOffset = kmerSize - 1 - firstOffset;
				assert(firstOffset < kmerSize);
				size_t secondOffset = i;
				if (!to.second) secondOffset = kmerSize - 1 - secondOffset;
				assert(secondOffset < kmerSize);
				merge(parent, rank, from.first, firstOffset, to.first, secondOffset);
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
				size_t firstOffset = kmerSize - overlap + i;
				assert(firstOffset < kmerSize);
				if (!from.second) firstOffset = kmerSize - 1 - firstOffset;
				assert(firstOffset < kmerSize);
				size_t secondOffset = i;
				if (!to.second) secondOffset = kmerSize - 1 - secondOffset;
				assert(secondOffset < kmerSize);
				merge(parent, rank, from.first, firstOffset, to.first, secondOffset);
			}
		}
	}
	std::unordered_set<NodeType> keys;
	for (auto pair : parent)
	{
		keys.emplace(pair.first);
	}
	for (auto node : keys)
	{
		for (size_t i = 0; i < kmerSize; i++)
		{
			auto found = find(parent, std::make_pair(node, i));
			size_t length = hashlist.getRunLength(found.first, found.second);
			assert(length > 0);
			hashlist.setRunLength(node, i, length);
		}
	}
}

void collectApproxNeighbors(std::unordered_map<uint64_t, unsigned char>& approxNeighbors, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
{
	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
	if (!toKmer.second) toSequence = revCompRLE(toSequence);
	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
	assert(edgeSequence.size() > kmerSize);
	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
	assert(edgeSequence.size() < 2 * kmerSize);
	if (edgeSequence.size() == kmerSize + 1) return;
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
		if (approxNeighbors.count(newBwHash) == 1 && approxNeighbors.at(newBwHash) != revCompRLE(edgeSequence[i-kmerSize])) approxNeighbors[newBwHash] = std::numeric_limits<unsigned char>::max();
		if (approxNeighbors.count(newBwHash) == 0) approxNeighbors[newBwHash] = revCompRLE(edgeSequence[i-kmerSize]);
	}
}

void collectExactNeighbors(const std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::unordered_map<HashType, unsigned char>& exactNeighbors, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
{
	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
	if (!toKmer.second) toSequence = revCompRLE(toSequence);
	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
	assert(edgeSequence.size() > kmerSize);
	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
	assert(edgeSequence.size() < 2 * kmerSize);
	if (edgeSequence.size() == kmerSize + 1) return;
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
		assert(approxNeighbors.count(oldFwHash) == 1);
		assert(approxNeighbors.count(newBwHash) == 1);
		// neither is a junction node if both not in approx
		if (approxNeighbors.at(oldFwHash) != std::numeric_limits<unsigned char>::max() && approxNeighbors.at(newBwHash) != std::numeric_limits<unsigned char>::max()) continue;
		HashType oldFwExactHash = hash(edgeSequence.substr(i-kmerSize, kmerSize));
		HashType newBwExactHash = hash(revCompRLE(edgeSequence.substr(i-kmerSize+1, kmerSize)));
		if (exactNeighbors.count(oldFwExactHash) == 1 && exactNeighbors.at(oldFwExactHash) != edgeSequence[i]) exactNeighbors[oldFwExactHash] = std::numeric_limits<unsigned char>::max();
		if (exactNeighbors.count(oldFwExactHash) == 0) exactNeighbors[oldFwExactHash] = edgeSequence[i];
		if (exactNeighbors.count(newBwExactHash) == 1 && exactNeighbors.at(newBwExactHash) != revCompRLE(edgeSequence[i-kmerSize])) exactNeighbors[newBwExactHash] = std::numeric_limits<unsigned char>::max();
		if (exactNeighbors.count(newBwExactHash) == 0) exactNeighbors[newBwExactHash] = revCompRLE(edgeSequence[i-kmerSize]);
	}
}

void collectJunctionHashes(const std::unordered_map<uint64_t, unsigned char>& approxNeighbors, std::unordered_set<uint64_t>& approxHashes, HashList& hashlist, const size_t kmerSize, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
{
	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
	if (!toKmer.second) toSequence = revCompRLE(toSequence);
	size_t overlap = hashlist.getOverlap(fromKmer, toKmer);
	std::string edgeSequence = fromSequence + toSequence.substr(overlap);
	assert(edgeSequence.size() > kmerSize);
	assert(edgeSequence.size() == kmerSize + kmerSize - overlap);
	assert(edgeSequence.size() < 2 * kmerSize);
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

bool addJunctionsToHashes(const std::unordered_set<uint64_t>& approxHashes, std::unordered_map<NodeType, std::pair<size_t, bool>>& belongsToUnitig, std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>& addEdges, HashList& hashlist, UnitigGraph& unitigs, const size_t kmerSize, std::pair<size_t, bool> unitigPosition, std::pair<size_t, bool> toUnitigPosition, std::pair<NodeType, bool> fromKmer, std::pair<NodeType, bool> toKmer)
{
	std::string fromSequence = hashlist.getHashSequenceRLE(fromKmer.first).toString();
	if (!fromKmer.second) fromSequence = revCompRLE(fromSequence);
	std::string toSequence = hashlist.getHashSequenceRLE(toKmer.first).toString();
	if (!toKmer.second) toSequence = revCompRLE(toSequence);
	size_t bigOverlap = hashlist.getOverlap(fromKmer, toKmer); // actually the smallest overlap in the entire function
	assert(bigOverlap < kmerSize);
	std::string edgeSequence = fromSequence + toSequence.substr(bigOverlap);
	assert(edgeSequence.size() > kmerSize);
	assert(edgeSequence.size() < 2 * kmerSize);
	if (edgeSequence.size() == kmerSize + 1) return false;
	std::string revSeq = revCompRLE(edgeSequence);
	std::vector<uint16_t> lens;
	lens = hashlist.getHashCharacterLength(fromKmer.first);
	if (!fromKmer.second) std::reverse(lens.begin(), lens.end());
	std::vector<uint16_t> tolens;
	tolens = hashlist.getHashCharacterLength(toKmer.first);
	if (!toKmer.second) std::reverse(tolens.begin(), tolens.end());
	lens.insert(lens.end(), tolens.begin() + bigOverlap, tolens.end());
	assert(lens.size() == edgeSequence.size());
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(edgeSequence[i]);
	}
	std::pair<NodeType, bool> lastKmer = fromKmer;
	std::pair<NodeType, bool> lastUnitig = unitigPosition;
	size_t lastPosition = 0;
	bool addedAny = false;
	size_t edgeCoverage = unitigs.edgeCoverage(unitigPosition, toUnitigPosition);
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

		std::string_view minimizerSequence { edgeSequence.data() + kmerStartPos, kmerSize };
		size_t revPos = edgeSequence.size() - (kmerStartPos + kmerSize);
		std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
		std::pair<size_t, bool> current;
		std::tie(current, lastHash) = hashlist.addNode(minimizerSequence, revMinimizerSequence, lens, kmerStartPos, kmerStartPos+kmerSize, lastHash, overlap, minHash);
		hashlist.addSequenceOverlap(lastKmer, current, overlap);
		std::pair<size_t, bool> newUnitigPosition { std::numeric_limits<size_t>::max(), true };
		if (belongsToUnitig.count(current.first) == 0)
		{
			unitigs.unitigs.emplace_back();
			unitigs.unitigs.back().emplace_back(current.first, true);
			unitigs.unitigCoverage.emplace_back();
			unitigs.unitigCoverage.back().emplace_back(0);
			unitigs.unitigCoverage.back()[0] += edgeCoverage;
			unitigs.edges.emplace_back();
			unitigs.edgeCov.emplace_back();
			belongsToUnitig[current.first] = std::make_pair(unitigs.unitigs.size()-1, true);
			newUnitigPosition = std::make_pair(unitigs.unitigs.size()-1, true);
			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
		}
		else
		{
			newUnitigPosition = belongsToUnitig.at(current.first);
			if (!current.second) newUnitigPosition.second = !newUnitigPosition.second;
			assert(unitigs.unitigCoverage[newUnitigPosition.first].size() == 1);
			unitigs.unitigCoverage[newUnitigPosition.first][0] += edgeCoverage;
		}
		assert(newUnitigPosition.first != std::numeric_limits<size_t>::max());
		unitigs.edgeCoverage(lastUnitig, newUnitigPosition) += edgeCoverage;
		auto key = canon(lastUnitig, newUnitigPosition);
		addEdges.insert(key);
		lastUnitig = newUnitigPosition;
		lastKmer = current;
		lastPosition = kmerStartPos;
		addedAny = true;
	}
	if (lastPosition != 0)
	{
		assert(lastKmer != fromKmer);
		assert(lastUnitig != unitigPosition);
		assert(lastPosition != 0);
		assert(addedAny);
		size_t overlap = kmerSize - (edgeSequence.size() - kmerSize) + lastPosition;
		assert(overlap < kmerSize);
		hashlist.addSequenceOverlap(lastKmer, toKmer, overlap);
		unitigs.edgeCoverage(lastUnitig, toUnitigPosition) += edgeCoverage;
		auto key = canon(lastUnitig, toUnitigPosition);
		addEdges.insert(key);
	}
	return addedAny;
}

void forceEdgeDeterminism(HashList& reads, UnitigGraph& unitigs, const size_t kmerSize, const double minUnitigCoverage)
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
			collectApproxNeighbors(approxNeighbors, reads, kmerSize, fromKmer, toKmer);
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			collectApproxNeighbors(approxNeighbors, reads, kmerSize, fromKmer, toKmer);
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
			collectJunctionHashes(approxNeighbors, approxJunctionHashes, reads, kmerSize, fromKmer, toKmer);
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			collectJunctionHashes(approxNeighbors, approxJunctionHashes, reads, kmerSize, fromKmer, toKmer);
		}
	}
	std::unordered_map<NodeType, std::pair<size_t, bool>> belongsToUnitig;
	std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> removeEdges;
	std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> addEdges;
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		auto fromKmer = unitigs.unitigs[i].back();
		for (auto edge : unitigs.edges[fw])
		{
			if (canon(fw, edge) != std::make_pair(fw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, addEdges, reads, unitigs, kmerSize, fw, edge, fromKmer, toKmer))
			{
				removeEdges.emplace(fw, edge);
			}
		}
		std::pair<size_t, bool> bw { i, false };
		fromKmer = reverse(unitigs.unitigs[i][0]);
		for (auto edge : unitigs.edges[bw])
		{
			if (canon(bw, edge) != std::make_pair(bw, edge)) continue;
			auto toKmer = unitigs.unitigs[edge.first][0];
			if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
			if (addJunctionsToHashes(approxJunctionHashes, belongsToUnitig, addEdges, reads, unitigs, kmerSize, bw, edge, fromKmer, toKmer))
			{
				removeEdges.emplace(bw, edge);
			}
		}
	}
	assert((addEdges.size() > 0) == (removeEdges.size() > 0));
	bool addedAny = removeEdges.size() > 0;
	for (auto pair : addEdges)
	{
		unitigs.edges[pair.first].insert(pair.second);
		unitigs.edges[reverse(pair.second)].insert(reverse(pair.first));
	}
	for (auto pair : removeEdges)
	{
		if (unitigs.edges[pair.first].count(pair.second) == 1)
		{
			unitigs.edges[pair.first].erase(pair.second);
		}
		if (unitigs.edges[reverse(pair.second)].count(reverse(pair.first)) == 1)
		{
			unitigs.edges[reverse(pair.second)].erase(reverse(pair.first));
		}
	}
	if (addedAny)
	{
		unitigs = getUnitigs(unitigs);
		if (minUnitigCoverage > 1)
		{
			size_t oldSize = unitigs.unitigs.size();
			unitigs = unitigs.filterUnitigsByCoverage(minUnitigCoverage);
			if (unitigs.unitigs.size() != oldSize) unitigs = getUnitigs(unitigs);
		}
	}
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool collapseRunLengths, const bool blunt, const size_t numThreads)
{
	auto beforeReading = getTime();
	HashList reads { kmerSize, collapseRunLengths, numThreads };
	loadReadsAsHashesMultithread(reads, inputReads, kmerSize, windowSize, hpc, collapseRunLengths, numThreads);
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage) unitigs = getUnitigs(unitigs.filterUnitigsByCoverage(minUnitigCoverage));
	auto beforeDeterminism = getTime();
	auto beforeConsistency = getTime();
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
		if (windowSize > 1) forceEdgeDeterminism(reads, unitigs, kmerSize, minUnitigCoverage);
		beforeConsistency = getTime();
		if (!collapseRunLengths) forceEdgeConsistency(unitigs, reads, kmerSize);
		beforeWrite = getTime();
		writeGraph(unitigs, outputGraph, reads, kmerSize);
		stats = getSizeAndN50(reads, unitigs, kmerSize, windowSize);
	}
	auto afterWrite = getTime();
	std::cerr << "reading and hashing sequences took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeDeterminism) << std::endl;
	if (!blunt && windowSize > 1) std::cerr << "forcing edge determinism took " << formatTime(beforeDeterminism, beforeConsistency) << std::endl;
	if (!blunt && !collapseRunLengths) std::cerr << "forcing edge consistency took " << formatTime(beforeConsistency, beforeWrite) << std::endl;
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
