#include <iostream>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <phmap.h>
#include "MBGCommon.h"
#include "UnitigResolver.h"
#include "RankBitvector.h"
#include "ReadHelper.h"

class ResolvableUnitigGraph
{
public:
	ResolvableUnitigGraph(const HashList& hashlist, const size_t kmerSize) :
	lastTippableChecked(0),
	hashlist(hashlist),
	kmerSize(kmerSize)
	{
	}
	std::vector<std::vector<std::pair<size_t, bool>>> unitigs;
	std::vector<size_t> unitigLeftClipBp;
	std::vector<size_t> unitigRightClipBp;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> overlaps;
	std::vector<bool> unitigRemoved;
	std::vector<std::unordered_set<size_t>> readsCrossingNode;
	std::vector<size_t> everTippable;
	size_t lastTippableChecked;
	size_t getBpOverlap(const std::pair<size_t, bool> from, const std::pair<size_t, bool> to) const
	{
		size_t kmerOverlap = overlaps.at(canon(from, to));
		size_t result;
		if (kmerOverlap == 0)
		{
			std::pair<size_t, bool> fromKmer = from.second ? unitigs[from.first].back() : reverse(unitigs[from.first][0]);
			std::pair<size_t, bool> toKmer = to.second ? unitigs[to.first][0] : reverse(unitigs[to.first].back());
			result = hashlist.getOverlap(fromKmer, toKmer);
		}
		else
		{
			result = kmerOverlap * kmerSize;
			for (size_t i = 1; i < kmerOverlap; i++)
			{
				std::pair<size_t, bool> oldKmer = unitigs[to.first][i-1];
				std::pair<size_t, bool> newKmer = unitigs[to.first][i];
				if (!to.second)
				{
					oldKmer = reverse(unitigs[to.first][unitigs[to.first].size()-i]);
					newKmer = reverse(unitigs[to.first][unitigs[to.first].size()-i-1]);
				}
				result -= hashlist.getOverlap(oldKmer, newKmer);
			}
		}
		size_t fromClip = from.second ? unitigRightClipBp[from.first] : unitigLeftClipBp[from.first];
		size_t toClip = to.second ? unitigLeftClipBp[to.first] : unitigRightClipBp[to.first];
		assert(result >= fromClip + toClip);
		result -= fromClip + toClip;
		return result;
	}
	size_t unitigLength(size_t i) const
	{
		size_t result = unitigs[i].size() * kmerSize;
		for (size_t j = 1; j < unitigs[i].size(); j++)
		{
			assert(result > hashlist.getOverlap(unitigs[i][j-1], unitigs[i][j]));
			result -= hashlist.getOverlap(unitigs[i][j-1], unitigs[i][j]);
		}
		assert(result >= kmerSize);
		assert(result > unitigLeftClipBp[i] + unitigRightClipBp[i]);
		result -= unitigLeftClipBp[i] + unitigRightClipBp[i];
		assert(result >= kmerSize);
		return result;
	}
private:
	const HashList& hashlist;
	const size_t kmerSize;
};

void compact(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& paths, std::vector<size_t>& queueNodes)
{
	{
		std::vector<std::unordered_set<size_t>> tmp;
		std::swap(resolvableGraph.readsCrossingNode, tmp);
	}
	RankBitvector kept { resolvableGraph.unitigs.size() };
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		kept.set(i, !resolvableGraph.unitigRemoved[i]);
	}
	kept.buildRanks();
	size_t newSize = kept.getRank(kept.size()-1) + (kept.get(kept.size()-1) ? 1 : 0);
	for (size_t i = resolvableGraph.everTippable.size()-1; i < resolvableGraph.everTippable.size(); i--)
	{
		if (!kept.get(resolvableGraph.everTippable[i]))
		{
			std::swap(resolvableGraph.everTippable[i], resolvableGraph.everTippable.back());
			resolvableGraph.everTippable.pop_back();
			continue;
		}
		resolvableGraph.everTippable[i] = kept.getRank(resolvableGraph.everTippable[i]);
	}
	resolvableGraph.lastTippableChecked = kept.getRank(resolvableGraph.lastTippableChecked);
	resolvableGraph.unitigRemoved.resize(newSize);
	for (size_t i = 0; i < resolvableGraph.unitigRemoved.size(); i++)
	{
		resolvableGraph.unitigRemoved[i] = false;
	}
	{
		std::vector<std::vector<std::pair<size_t, bool>>> newUnitigs;
		std::vector<size_t> newUnitigLeftClipBp;
		std::vector<size_t> newUnitigRightClipBp;
		newUnitigs.resize(newSize);
		newUnitigRightClipBp.resize(newSize);
		newUnitigLeftClipBp.resize(newSize);
		for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			assert(newIndex < newUnitigs.size());
			assert(newIndex < newUnitigLeftClipBp.size());
			assert(newIndex < newUnitigRightClipBp.size());
			std::swap(newUnitigs[newIndex], resolvableGraph.unitigs[i]);
			newUnitigLeftClipBp[newIndex] = resolvableGraph.unitigLeftClipBp[i];
			newUnitigRightClipBp[newIndex] = resolvableGraph.unitigRightClipBp[i];
		}
		std::swap(resolvableGraph.unitigs, newUnitigs);
		std::swap(resolvableGraph.unitigRightClipBp, newUnitigRightClipBp);
		std::swap(resolvableGraph.unitigLeftClipBp, newUnitigLeftClipBp);
	}
	{
		VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> newEdges;
		newEdges.resize(newSize);
		for (size_t i = 0; i < resolvableGraph.edges.size(); i++)
		{
			if (!kept.get(i)) continue;
			size_t newIndex = kept.getRank(i);
			assert(newIndex < newEdges.size());
			for (auto target : resolvableGraph.edges[std::make_pair(i, true)])
			{
				assert(kept.get(target.first));
				assert(resolvableGraph.edges[reverse(target)].count(reverse(std::make_pair(i, true))) == 1);
				newEdges[std::make_pair(newIndex, true)].emplace(kept.getRank(target.first), target.second);
			}
			for (auto target : resolvableGraph.edges[std::make_pair(i, false)])
			{
				assert(kept.get(target.first));
				assert(resolvableGraph.edges[reverse(target)].count(reverse(std::make_pair(i, false))) == 1);
				newEdges[std::make_pair(newIndex, false)].emplace(kept.getRank(target.first), target.second);
			}
		}
		std::swap(resolvableGraph.edges, newEdges);
	}
	{
		phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newOverlaps;
		for (auto pair : resolvableGraph.overlaps)
		{
			auto from = pair.first.first;
			auto to = pair.first.second;
			if (!kept.get(from.first) || !kept.get(to.first)) continue;
			auto overlap = pair.second;
			from.first = kept.getRank(from.first);
			to.first = kept.getRank(to.first);
			assert(from.first < resolvableGraph.unitigs.size());
			assert(to.first < resolvableGraph.unitigs.size());
			newOverlaps[canon(from, to)] = overlap;
		}
		std::swap(resolvableGraph.overlaps, newOverlaps);
	}
	for (size_t i = paths.size()-1; i < paths.size(); i--)
	{
		if (paths[i].path.size() == 0)
		{
			std::swap(paths[i], paths.back());
			paths.pop_back();
			continue;
		}
		for (size_t j = 0; j < paths[i].path.size(); j++)
		{
			assert(kept.get(paths[i].path[j].first));
			paths[i].path[j].first = kept.getRank(paths[i].path[j].first);
			assert(paths[i].path[j].first < resolvableGraph.unitigs.size());
		}
	}
	for (size_t i = queueNodes.size()-1; i < queueNodes.size(); i--)
	{
		if (!kept.get(queueNodes[i]))
		{
			std::swap(queueNodes[i], queueNodes.back());
			queueNodes.pop_back();
			continue;
		}
		queueNodes[i] = kept.getRank(queueNodes[i]);
	}
	resolvableGraph.readsCrossingNode.resize(newSize);
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::unordered_set<size_t> crossesNodes;
		for (auto pair : paths[i].path)
		{
			crossesNodes.insert(pair.first);
		}
		for (auto node : crossesNodes)
		{
			resolvableGraph.readsCrossingNode[node].insert(i);
		}
	}
	for (size_t i = 0; i < resolvableGraph.unitigRemoved.size(); i++)
	{
		assert(!resolvableGraph.unitigRemoved[i]);
	}
}

class UnitigLengthComparer
{
public:
	UnitigLengthComparer(const ResolvableUnitigGraph& resolvableGraph) :
	resolvableGraph(resolvableGraph)
	{
	}
	bool operator()(const size_t left, const size_t right) const
	{
		return resolvableGraph.unitigLength(left) > resolvableGraph.unitigLength(right);
	}
private:
	const ResolvableUnitigGraph& resolvableGraph;
};

ResolvableUnitigGraph getUnitigs(const UnitigGraph& initial, size_t minCoverage, const HashList& hashlist, const size_t kmerSize)
{
	ResolvableUnitigGraph result { hashlist, kmerSize };
	result.unitigs.resize(initial.unitigs.size());
	result.unitigLeftClipBp.resize(initial.unitigs.size(), 0);
	result.unitigRightClipBp.resize(initial.unitigs.size(), 0);
	result.unitigRemoved.resize(initial.unitigs.size(), false);
	result.edges.resize(initial.unitigs.size());
	result.readsCrossingNode.resize(initial.unitigs.size());
	for (size_t i = 0; i < initial.unitigs.size(); i++)
	{
		result.unitigs[i].insert(result.unitigs[i].end(), initial.unitigs[i].begin(), initial.unitigs[i].end());
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : initial.edgeCov[fw])
		{
			if (pair.second < minCoverage) continue;
			auto canonpair = canon(fw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[fw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(fw));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : initial.edgeCov[bw])
		{
			if (pair.second < minCoverage) continue;
			auto canonpair = canon(bw, pair.first);
			result.overlaps[canonpair] = 0;
			result.edges[bw].emplace(pair.first);
			result.edges[reverse(pair.first)].emplace(reverse(bw));
		}
	}
	return result;
}

std::vector<ReadPath> getUnitigPaths(const ResolvableUnitigGraph& graph, const HashList& hashlist, const std::vector<std::string>& readFiles, const size_t numThreads, const ReadpartIterator& partIterator, const size_t kmerSize)
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
			size_t lastReadPos = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
			std::tuple<size_t, size_t, bool> lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			for (const size_t readPos : positions)
			{
				VectorView<CharType> minimizerSequence { seq, readPos, readPos + kmerSize };
				std::pair<size_t, bool> kmer = hashlist.getNodeOrNull(minimizerSequence);
				if (kmer.first == std::numeric_limits<size_t>::max()) continue;
				size_t readPosExpandedStart = poses[readPos];
				size_t readPosExpandedEnd = poses[readPos+kmerSize];
				if (kmer.first >= kmerLocator.size() || std::get<0>(kmerLocator[kmer.first]) == std::numeric_limits<size_t>::max())
				{
					if (current.path.size() > 0)
					{
						std::lock_guard<std::mutex> lock { resultMutex };
						result.emplace_back();
						std::swap(result.back(), current);
					}
					current.path.clear();
					current.readName = read.seq_id;
					current.readLength = rawSeq.size();
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
					current.readPosesExpandedStart.push_back(readPosExpandedStart);
					current.readPosesExpandedEnd.push_back(readPosExpandedEnd);
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
							std::lock_guard<std::mutex> lock { resultMutex };
							result.emplace_back();
							std::swap(result.back(), current);
						}
						current.path.clear();
						current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
						current.readPoses.clear();
						current.readPosesExpandedStart.clear();
						current.readPosesExpandedEnd.clear();
						current.readPoses.push_back(readPos);
						current.readPosesExpandedStart.push_back(readPosExpandedStart);
						current.readPosesExpandedEnd.push_back(readPosExpandedEnd);
						current.leftClip = std::get<1>(pos);
						current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
						current.readName = read.seq_id;
						current.readLength = rawSeq.size();
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
					current.readPosesExpandedStart.push_back(readPosExpandedStart);
					current.readPosesExpandedEnd.push_back(readPosExpandedEnd);
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
					current.readPosesExpandedStart.push_back(readPosExpandedStart);
					current.readPosesExpandedEnd.push_back(readPosExpandedEnd);
					current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1;
					lastPos = pos;
					lastReadPos = readPos;
					lastKmer = kmer;
					continue;
				}
				assert(current.path.size() > 0);
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					result.emplace_back();
					std::swap(result.back(), current);
				}
				current.path.clear();
				current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
				current.readPoses.clear();
				current.readPosesExpandedStart.clear();
				current.readPosesExpandedEnd.clear();
				current.readPoses.push_back(readPos);
				current.readPosesExpandedStart.push_back(readPosExpandedStart);
				current.readPosesExpandedEnd.push_back(readPosExpandedEnd);
				current.leftClip = std::get<1>(pos);
				current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
				current.readName = read.seq_id;
				current.readLength = rawSeq.size();
				assert(current.leftClip + current.rightClip + 1 == graph.unitigs[std::get<0>(pos)].size());
				lastPos = pos;
				lastReadPos = readPos;
				lastKmer = kmer;
				continue;
			}
			if (current.path.size() > 0)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				result.emplace_back();
				std::swap(result.back(), current);
			}
		});
	});

	return result;
}

std::pair<UnitigGraph, std::vector<ReadPath>> resolvableToUnitigs(const ResolvableUnitigGraph& resolvableGraph, const std::vector<ReadPath>& readPaths)
{
	UnitigGraph result;
	RankBitvector newIndex { resolvableGraph.unitigs.size() };
	assert(resolvableGraph.unitigs.size() == resolvableGraph.unitigRemoved.size());
	assert(resolvableGraph.unitigs.size() == resolvableGraph.edges.size());
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		newIndex.set(i, !resolvableGraph.unitigRemoved[i]);
	}
	newIndex.buildRanks();
	size_t newSize = newIndex.getRank(newIndex.size()-1) + (newIndex.get(newIndex.size()-1) ? 1 : 0);
	result.unitigs.resize(newSize);
	result.leftClip.resize(newSize);
	result.rightClip.resize(newSize);
	result.unitigCoverage.resize(newSize);
	result.edges.resize(newSize);
	result.edgeCov.resize(newSize);
	result.edgeOvlp.resize(newSize);
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (!newIndex.get(i)) continue;
		size_t unitig = newIndex.getRank(i);
		result.unitigs[unitig].insert(result.unitigs[unitig].end(), resolvableGraph.unitigs[i].begin(), resolvableGraph.unitigs[i].end());
		for (size_t j = 0; j < resolvableGraph.unitigs[i].size(); j++)
		{
			result.unitigCoverage[unitig].push_back(0);
		}
		result.leftClip[unitig] = resolvableGraph.unitigLeftClipBp[i];
		result.rightClip[unitig] = resolvableGraph.unitigRightClipBp[i];
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> newfw { newIndex.getRank(i), true };
		for (auto edge : resolvableGraph.edges[fw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			result.edges[newfw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newfw));
			result.edgeCoverage(newfw, newEdge) = 0;
			result.edgeOverlap(newfw, newEdge) = resolvableGraph.overlaps.at(canon(fw, edge));
		}
		std::pair<size_t, bool> bw { i, false };
		std::pair<size_t, bool> newbw { newIndex.getRank(i), false };
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(newIndex.get(edge.first));
			std::pair<size_t, bool> newEdge { newIndex.getRank(edge.first), edge.second };
			result.edges[newbw].emplace(newEdge);
			result.edges[reverse(newEdge)].emplace(reverse(newbw));
			result.edgeCoverage(newbw, newEdge) = 0;
			result.edgeOverlap(newbw, newEdge) = resolvableGraph.overlaps.at(canon(bw, edge));
		}
	}
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		if (path.path.size() == 1)
		{
			assert(newIndex.get(path.path[0].first));
			size_t unitig = newIndex.getRank(path.path[0].first);
			for (size_t i = path.leftClip; i < result.unitigCoverage[unitig].size() - path.rightClip; i++)
			{
				size_t index = i;
				if (!path.path[0].second) index = result.unitigCoverage[unitig].size() - 1 - i;
				result.unitigCoverage[unitig][index] += 1;
			}
		}
		else
		{
			for (size_t i = 1; i < path.path.size()-1; i++)
			{
				assert(newIndex.get(path.path[i].first));
				size_t unitig = newIndex.getRank(path.path[i].first);
				for (size_t j = 0; j < result.unitigCoverage[unitig].size(); j++)
				{
					result.unitigCoverage[unitig][j] += 1;
				}
			}
			assert(newIndex.get(path.path[0].first));
			size_t unitig = newIndex.getRank(path.path[0].first);
			for (size_t i = path.leftClip; i < result.unitigCoverage[unitig].size(); i++)
			{
				size_t index = i;
				if (!path.path[0].second) index = result.unitigCoverage[unitig].size()-1-i;
				result.unitigCoverage[unitig][index] += 1;
			}
			assert(newIndex.get(path.path.back().first));
			unitig = newIndex.getRank(path.path.back().first);
			for (size_t i = 0; i < result.unitigCoverage[unitig].size() - path.rightClip; i++)
			{
				size_t index = i;
				if (!path.path.back().second) index = result.unitigCoverage[unitig].size()-1-i;
				result.unitigCoverage[unitig][index] += 1;
			}
		}
	}
	for (size_t i = 0; i < result.unitigCoverage.size(); i++)
	{
		for (size_t j = 0; j < result.unitigCoverage[i].size(); j++)
		{
			assert(result.unitigCoverage[i][j] > 0);
		}
	}
	std::vector<ReadPath> resultReads;
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		resultReads.push_back(path);
		for (size_t i = 0; i < resultReads.back().path.size(); i++)
		{
			assert(newIndex.get(resultReads.back().path[i].first));
			resultReads.back().path[i].first = newIndex.getRank(resultReads.back().path[i].first);
		}
	}
	return std::make_pair(result, resultReads);
}

std::vector<std::pair<size_t, bool>> extend(const ResolvableUnitigGraph& resolvableGraph, const std::pair<size_t, bool> start)
{
	std::pair<size_t, bool> pos = start;
	std::vector<std::pair<size_t, bool>> result;
	result.emplace_back(pos);
	assert(!resolvableGraph.unitigRemoved[start.first]);
	while (true)
	{
		assert(!resolvableGraph.unitigRemoved[pos.first]);
		if (resolvableGraph.edges[pos].size() != 1) break;
		auto newPos = *resolvableGraph.edges[pos].begin();
		auto revNewPos = reverse(newPos);
		if (resolvableGraph.edges[revNewPos].size() != 1) break;
		if (newPos == start)
		{
			// include the start to note that it's circular
			result.emplace_back(newPos);
			break;
		}
		if (newPos.first == pos.first) break;
		pos = newPos;
		result.emplace_back(pos);
	}
	return result;
}

size_t getNumberOfHashes(const ResolvableUnitigGraph& resolvableGraph, size_t leftClip, size_t rightClip, const std::vector<std::pair<size_t, bool>>& path)
{
	if (path.size() == 0) return 0;
	size_t result = 0;
	assert(path.size() < 2 || leftClip < resolvableGraph.unitigs[path[0].first].size());
	assert(path.size() < 2 || rightClip < resolvableGraph.unitigs[path.back().first].size());
	assert(path.size() != 1 || leftClip + rightClip < resolvableGraph.unitigs[path[0].first].size());
	for (size_t i = 0; i < path.size(); i++)
	{
		result += resolvableGraph.unitigs[path[i].first].size();
		if (i > 0)
		{
			assert(resolvableGraph.unitigs[path[i].first].size() >= resolvableGraph.overlaps.at(canon(path[i-1], path[i])));
			result -= resolvableGraph.overlaps.at(canon(path[i-1], path[i]));
		}
	}
	assert(result > leftClip + rightClip);
	result -= leftClip;
	result -= rightClip;
	return result;
}

void erasePath(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, size_t i)
{
	std::unordered_set<size_t> nodes;
	for (auto pos : readPaths[i].path)
	{
		nodes.emplace(pos.first);
	}
	for (auto node : nodes)
	{
		assert(resolvableGraph.readsCrossingNode[node].count(i) == 1);
		resolvableGraph.readsCrossingNode[node].erase(i);
	}
	readPaths[i].path.clear();
}

void addPath(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, ReadPath&& newPath)
{
	if (newPath.path.size() == 0) return;
	assert(newPath.readPosesExpandedStart.size() == newPath.readPoses.size());
	assert(newPath.readPosesExpandedEnd.size() == newPath.readPoses.size());
	assert(getNumberOfHashes(resolvableGraph, newPath.leftClip, newPath.rightClip, newPath.path) == newPath.readPoses.size());
	if (newPath.path.size() == 1)
	{
		assert(resolvableGraph.unitigs[newPath.path[0].first].size() > newPath.leftClip + newPath.rightClip);
	}
	else
	{
		assert(resolvableGraph.unitigs[newPath.path[0].first].size() > newPath.leftClip);
		assert(resolvableGraph.unitigs[newPath.path.back().first].size() > newPath.rightClip);
	}
	for (size_t i = 0; i < newPath.path.size(); i++)
	{
		assert(!resolvableGraph.unitigRemoved[newPath.path[i].first]);
		if (i > 0)
		{
			assert(resolvableGraph.edges[newPath.path[i-1]].count(newPath.path[i]) == 1);
			assert(resolvableGraph.edges[reverse(newPath.path[i])].count(reverse(newPath.path[i-1])) == 1);
		}
	}
	size_t newIndex = readPaths.size();
	std::unordered_set<size_t> nodes;
	for (auto pos : newPath.path)
	{
		nodes.insert(pos.first);
	}
	readPaths.emplace_back(std::move(newPath));
	for (auto node : nodes)
	{
		assert(resolvableGraph.readsCrossingNode[node].count(newIndex) == 0);
		resolvableGraph.readsCrossingNode[node].emplace(newIndex);
	}
}

void replacePathNodes(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const std::vector<std::pair<size_t, bool>>& newUnitig, size_t newUnitigIndex)
{
	assert(newUnitig.size() >= 2);
	std::unordered_set<size_t> relevantReads;
	for (auto pos : newUnitig)
	{
		relevantReads.insert(resolvableGraph.readsCrossingNode[pos.first].begin(), resolvableGraph.readsCrossingNode[pos.first].end());
	}
	std::unordered_map<size_t, size_t> leftClip;
	std::unordered_map<size_t, size_t> rightClip;
	assert(resolvableGraph.unitigs[newUnitig[0].first].size() <= resolvableGraph.unitigs[newUnitigIndex].size());
	assert(resolvableGraph.unitigLength(newUnitig[0].first) < resolvableGraph.unitigLength(newUnitigIndex));
	size_t leftClipSum = 0;
	size_t rightClipSum = resolvableGraph.unitigs[newUnitigIndex].size() - resolvableGraph.unitigs[newUnitig[0].first].size();
	leftClip[newUnitig[0].first] = leftClipSum;
	rightClip[newUnitig[0].first] = rightClipSum;
	assert(leftClipSum + rightClipSum + resolvableGraph.unitigs[newUnitig[0].first].size() == resolvableGraph.unitigs[newUnitigIndex].size());
	for (size_t i = 1; i < newUnitig.size(); i++)
	{
		assert(leftClip.count(newUnitig[i].first) == 0);
		assert(rightClip.count(newUnitig[i].first) == 0);
		assert(resolvableGraph.unitigs[newUnitig[i-1].first].size() >= resolvableGraph.overlaps[canon(newUnitig[i-1], newUnitig[i])]);
		assert(resolvableGraph.unitigs[newUnitig[i].first].size() >= resolvableGraph.overlaps[canon(newUnitig[i-1], newUnitig[i])]);
		assert(resolvableGraph.unitigs[newUnitig[i].first].size() - resolvableGraph.overlaps[canon(newUnitig[i-1], newUnitig[i])] <= rightClipSum);
		leftClipSum = leftClipSum + resolvableGraph.unitigs[newUnitig[i-1].first].size() - resolvableGraph.overlaps[canon(newUnitig[i-1], newUnitig[i])];
		rightClipSum = rightClipSum - resolvableGraph.unitigs[newUnitig[i].first].size() + resolvableGraph.overlaps[canon(newUnitig[i-1], newUnitig[i])];
		leftClip[newUnitig[i].first] = leftClipSum;
		rightClip[newUnitig[i].first] = rightClipSum;
		assert(leftClipSum + rightClipSum + resolvableGraph.unitigs[newUnitig[i].first].size() == resolvableGraph.unitigs[newUnitigIndex].size());
	}
	assert(rightClipSum == 0);
	std::unordered_map<size_t, bool> nodeForwardInUnitig;
	std::unordered_set<size_t> nodeInUnitig;
	for (auto pos : newUnitig)
	{
		assert(nodeInUnitig.count(pos.first) == 0);
		nodeInUnitig.emplace(pos.first);
		nodeForwardInUnitig[pos.first] = pos.second;
	}
	for (size_t i : relevantReads)
	{
		ReadPath newPath;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (nodeInUnitig.count(readPaths[i].path[j].first) == 0)
			{
				newPath.path.emplace_back(readPaths[i].path[j]);
				continue;
			}
			if (j == 0)
			{
				bool fw = nodeForwardInUnitig.at(readPaths[i].path[j].first);
				if (!readPaths[i].path[j].second) fw = !fw;
				newPath.path.emplace_back(newUnitigIndex, fw);
			}
			else if (readPaths[i].path[j] == newUnitig[0])
			{
				newPath.path.emplace_back(newUnitigIndex, true);
			}
			else if (readPaths[i].path[j] == reverse(newUnitig.back()))
			{
				newPath.path.emplace_back(newUnitigIndex, false);
			}
		}
		newPath.leftClip = readPaths[i].leftClip;
		newPath.rightClip = readPaths[i].rightClip;
		newPath.readPoses = readPaths[i].readPoses;
		newPath.readPosesExpandedStart = readPaths[i].readPosesExpandedStart;
		newPath.readPosesExpandedEnd = readPaths[i].readPosesExpandedEnd;
		newPath.readName = readPaths[i].readName;
		newPath.readLength = readPaths[i].readLength;
		if (nodeInUnitig.count(readPaths[i].path[0].first) == 1)
		{
			bool fw = nodeForwardInUnitig[readPaths[i].path[0].first];
			if (!readPaths[i].path[0].second) fw = !fw;
			if (fw)
			{
				newPath.leftClip += leftClip[readPaths[i].path[0].first];
			}
			else
			{
				newPath.leftClip += rightClip[readPaths[i].path[0].first];
			}
		}
		if (nodeInUnitig.count(readPaths[i].path.back().first) == 1)
		{
			bool fw = nodeForwardInUnitig[readPaths[i].path.back().first];
			if (!readPaths[i].path.back().second) fw = !fw;
			if (fw)
			{
				newPath.rightClip += rightClip[readPaths[i].path.back().first];
			}
			else
			{
				newPath.rightClip += leftClip[readPaths[i].path.back().first];
			}
		}
		addPath(resolvableGraph, readPaths, std::move(newPath));
		erasePath(resolvableGraph, readPaths, i);
	}
}

size_t unitigifyOne(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const size_t unitig)
{
	std::vector<std::pair<size_t, bool>> fwExtension = extend(resolvableGraph, std::make_pair(unitig, true));
	std::vector<std::pair<size_t, bool>> bwExtension = extend(resolvableGraph, std::make_pair(unitig, false));
	assert(bwExtension.size() >= 1);
	assert(fwExtension.size() >= 1);
	std::vector<std::pair<size_t, bool>> newUnitig = revCompPath(bwExtension);
	if (newUnitig.size() >= 2 && newUnitig[0] == newUnitig.back())
	{
		// circular
		assert(newUnitig.size() == fwExtension.size());
		for (size_t i = 0; i < fwExtension.size(); i++)
		{
			assert(newUnitig[i] == fwExtension[i]);
		}
		newUnitig.pop_back();
	}
	else
	{
		assert(newUnitig.back() == std::make_pair(unitig, true));
		newUnitig.insert(newUnitig.end(), fwExtension.begin()+1, fwExtension.end());
	}
	assert(newUnitig.size() >= 1);
	if (newUnitig.size() == 1) return 1;
	size_t newIndex = resolvableGraph.unitigs.size();
	size_t leftClip, rightClip;
	if (newUnitig[0].second)
	{
		leftClip = resolvableGraph.unitigLeftClipBp[newUnitig[0].first];
	}
	else
	{
		leftClip = resolvableGraph.unitigRightClipBp[newUnitig[0].first];
	}
	if (newUnitig.back().second)
	{
		rightClip = resolvableGraph.unitigRightClipBp[newUnitig.back().first];
	}
	else
	{
		rightClip = resolvableGraph.unitigLeftClipBp[newUnitig.back().first];
	}
	resolvableGraph.unitigLeftClipBp.push_back(leftClip);
	resolvableGraph.unitigRightClipBp.push_back(rightClip);
	resolvableGraph.unitigs.emplace_back();
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	for (size_t i = 0; i < newUnitig.size(); i++)
	{
		assert(i == 0 || resolvableGraph.edges[reverse(newUnitig[i])].size() == 1);
		assert(i == newUnitig.size()-1 || resolvableGraph.edges[newUnitig[i]].size() == 1);
		std::vector<std::pair<size_t, bool>> add = resolvableGraph.unitigs[newUnitig[i].first];
		if (!newUnitig[i].second) add = revCompPath(add);
		size_t overlap = 0;
		if (i > 0) overlap = resolvableGraph.overlaps.at(canon(newUnitig[i-1], newUnitig[i]));
		resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin()+overlap, add.end());
	}
	std::unordered_set<size_t> nodesInUnitig;
	for (auto pos : newUnitig)
	{
		assert(nodesInUnitig.count(pos.first) == 0);
		nodesInUnitig.insert(pos.first);
	}
	std::pair<size_t, bool> bw { newIndex, false };
	std::pair<size_t, bool> fw { newIndex, true };
	for (auto edge : resolvableGraph.edges[reverse(newUnitig[0])])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(bw, bw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
		else if (edge == newUnitig[0])
		{
			resolvableGraph.overlaps[canon(bw, fw)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(fw);
			resolvableGraph.edges[reverse(fw)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(bw, edge)] = resolvableGraph.overlaps.at(canon(reverse(newUnitig[0]), edge));
			resolvableGraph.edges[bw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(bw));
			assert(resolvableGraph.edges[reverse(edge)].count(newUnitig[0]) == 1);
			resolvableGraph.edges[reverse(edge)].erase(newUnitig[0]);
		}
	}
	for (auto edge : resolvableGraph.edges[newUnitig.back()])
	{
		assert(!resolvableGraph.unitigRemoved[edge.first]);
		if (edge == reverse(newUnitig.back()))
		{
			resolvableGraph.overlaps[canon(fw, bw)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(bw);
			resolvableGraph.edges[reverse(bw)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(newUnitig.back()));
		}
		else
		{
			assert(nodesInUnitig.count(edge.first) == 0);
			resolvableGraph.overlaps[canon(fw, edge)] = resolvableGraph.overlaps.at(canon(newUnitig.back(), edge));
			resolvableGraph.edges[fw].emplace(edge);
			resolvableGraph.edges[reverse(edge)].emplace(reverse(fw));
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(newUnitig.back())) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(newUnitig.back()));
		}
	}
	for (size_t i = 0; i < newUnitig.size(); i++)
	{
		resolvableGraph.unitigRemoved[newUnitig[i].first] = true;
		resolvableGraph.edges[newUnitig[i]].clear();
		resolvableGraph.edges[reverse(newUnitig[i])].clear();
	}
	replacePathNodes(resolvableGraph, readPaths, newUnitig, newIndex);
	return newUnitig.size();
}

void unitigifyAll(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths)
{
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		unitigifyOne(resolvableGraph, readPaths, i);
	}
}

void unresolveRecursively(const ResolvableUnitigGraph& resolvableGraph, const std::unordered_set<size_t>& resolvables, std::unordered_set<size_t>& unresolvables, const size_t node)
{
	unresolvables.emplace(node);
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		if (resolvables.count(edge.first) == 0) continue;
		if (unresolvables.count(edge.first) == 1) continue;
		unresolveRecursively(resolvableGraph, resolvables, unresolvables, edge.first);
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		if (resolvables.count(edge.first) == 0) continue;
		if (unresolvables.count(edge.first) == 1) continue;
		unresolveRecursively(resolvableGraph, resolvables, unresolvables, edge.first);
	}
}

void createEdgeNode(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, const size_t kmerSize, std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes, const std::unordered_set<size_t>& resolvables, const std::unordered_set<size_t>& unresolvables, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	size_t newIndex = resolvableGraph.unitigs.size();
	newEdgeNodes[std::make_pair(from, to)] = newIndex;
	resolvableGraph.unitigs.emplace_back();
	std::vector<std::pair<size_t, bool>> add = resolvableGraph.unitigs[from.first];
	size_t leftClipBp = resolvableGraph.unitigLeftClipBp[from.first];
	size_t rightClipBp = resolvableGraph.unitigRightClipBp[from.first];
	if (!from.second)
	{
		add = revCompPath(add);
		std::swap(leftClipBp, rightClipBp);
	}
	resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin(), add.end());
	size_t overlapIncrement = 0;
	if (resolvables.count(to.first) == 1 && unresolvables.count(to.first) == 0)
	{
		size_t start = resolvableGraph.overlaps.at(canon(from, to));
		assert(start <= resolvableGraph.unitigs[to.first].size());
		add = resolvableGraph.unitigs[to.first];
		if (!to.second) add = revCompPath(add);
		resolvableGraph.unitigs.back().insert(resolvableGraph.unitigs.back().end(), add.begin()+start, add.end());
		if (to.second)
		{
			rightClipBp = resolvableGraph.unitigRightClipBp[to.first];
		}
		else
		{
			rightClipBp = resolvableGraph.unitigLeftClipBp[to.first];
		}
	}
	else
	{
		size_t availableIncrease = resolvableGraph.unitigLength(to.first) - resolvableGraph.getBpOverlap(from, to) - 1;
		assert(availableIncrease >= 1);
		if (rightClipBp == 0)
		{
			size_t start = resolvableGraph.overlaps.at(canon(from, to));
			assert(start < resolvableGraph.unitigs[to.first].size());
			add = resolvableGraph.unitigs[to.first];
			if (!to.second) add = revCompPath(add);
			overlapIncrement = 1;
			for (size_t i = 0; i < start; i++)
			{
				assert(resolvableGraph.unitigs.back()[resolvableGraph.unitigs.back().size() - start+i] == add[i]);
			}
			assert(add.size() > start);
			rightClipBp = kmerSize - hashlist.getOverlap(resolvableGraph.unitigs.back().back(), add[start]) - 1;
			assert(rightClipBp < kmerSize);
			if (rightClipBp > availableIncrease - 1)
			{
				rightClipBp -= availableIncrease - 1;
			}
			else if (rightClipBp > 1)
			{
				rightClipBp = 1;
			}
			resolvableGraph.unitigs.back().emplace_back(add[start]);
		}
		else
		{
			assert(rightClipBp < kmerSize);
			if (rightClipBp > availableIncrease)
			{
				rightClipBp -= availableIncrease;
			}
			else
			{
				rightClipBp = 0;
			}
		}
	}
	resolvableGraph.unitigRightClipBp.push_back(rightClipBp);
	resolvableGraph.unitigLeftClipBp.push_back(leftClipBp);
	resolvableGraph.edges.emplace_back();
	resolvableGraph.unitigRemoved.emplace_back(false);
	resolvableGraph.readsCrossingNode.emplace_back();
	if (resolvables.count(to.first) == 0 || unresolvables.count(to.first) == 1)
	{
		resolvableGraph.edges[std::make_pair(newIndex, true)].emplace(to);
		resolvableGraph.edges[reverse(to)].emplace(std::make_pair(newIndex, false));
		assert(resolvableGraph.unitigs[newIndex].size() >= resolvableGraph.overlaps.at(canon(from, to)) + overlapIncrement);
		assert(resolvableGraph.unitigs[to.first].size() >= resolvableGraph.overlaps.at(canon(from, to)) + overlapIncrement);
		resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = resolvableGraph.overlaps.at(canon(from, to)) + overlapIncrement;
		// resolvableGraph.overlaps[canon(std::make_pair(newIndex, true), to)] = resolvableGraph.unitigs[from.first].size();
		assert(resolvableGraph.getBpOverlap(std::make_pair(newIndex, true), to) < resolvableGraph.unitigLength(to.first));
	}
	assert(resolvableGraph.unitigLength(newIndex) > resolvableGraph.unitigLength(from.first));
}

std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> getValidTriplets(const ResolvableUnitigGraph& resolvableGraph, const std::unordered_set<size_t>& resolvables, const std::vector<ReadPath>& readPaths, size_t node, size_t minCoverage)
{
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> empty;
	std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> tripletCoverage;
	for (size_t i : resolvableGraph.readsCrossingNode[node])
	{
		for (size_t j = 1; j < readPaths[i].path.size()-1; j++)
		{
			if (readPaths[i].path[j].first != node) continue;
			std::pair<size_t, bool> left;
			std::pair<size_t, bool> right;
			if (readPaths[i].path[j].second)
			{
				left = readPaths[i].path[j-1];
				right = readPaths[i].path[j+1];
			}
			else
			{
				left = reverse(readPaths[i].path[j+1]);
				right = reverse(readPaths[i].path[j-1]);
			}
			tripletCoverage[std::make_pair(left, right)] += 1;
			assert(!resolvableGraph.unitigRemoved[node]);
			assert(!resolvableGraph.unitigRemoved[left.first]);
			assert(!resolvableGraph.unitigRemoved[right.first]);
			assert(resolvableGraph.edges[left].count(std::make_pair(node, true)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, false)].count(reverse(left)) == 1);
			assert(resolvableGraph.edges[reverse(right)].count(std::make_pair(node, false)) == 1);
			assert(resolvableGraph.edges[std::make_pair(node, true)].count(right) == 1);
		}
	}
	std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> coveredTriplets;
	for (auto pair : tripletCoverage)
	{
		if (pair.second < minCoverage) continue;
		coveredTriplets.push_back(pair.first);
	}
	std::unordered_set<std::pair<size_t, bool>> coveredInNeighbors;
	std::unordered_set<std::pair<size_t, bool>> coveredOutNeighbors;
	for (auto pair : coveredTriplets)
	{
		coveredInNeighbors.emplace(reverse(pair.first));
		coveredOutNeighbors.emplace(pair.second);
	}
	assert(coveredInNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, false)].size());
	assert(coveredOutNeighbors.size() <= resolvableGraph.edges[std::make_pair(node, true)].size());
	if (coveredInNeighbors.size() < resolvableGraph.edges[std::make_pair(node, false)].size()) return empty;
	if (coveredOutNeighbors.size() < resolvableGraph.edges[std::make_pair(node, true)].size()) return empty;
	return coveredTriplets;
}

void replacePaths(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const std::unordered_set<size_t>& resolvables, const std::unordered_set<size_t>& unresolvables, const std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t>& newEdgeNodes)
{
	std::unordered_set<size_t> relevantReads;
	std::unordered_set<size_t> nodesInNewEdgeNodes;
	for (auto pair : newEdgeNodes)
	{
		assert(nodesInNewEdgeNodes.count(pair.second) == 0);
		assert(resolvables.count(pair.first.first.first) == 1);
		assert(unresolvables.count(pair.first.first.first) == 0);
		nodesInNewEdgeNodes.insert(pair.second);
	}
	for (auto node : resolvables)
	{
		if (unresolvables.count(node) == 1) continue;
		relevantReads.insert(resolvableGraph.readsCrossingNode[node].begin(), resolvableGraph.readsCrossingNode[node].end());
	}
	for (const size_t i : relevantReads)
	{
		ReadPath newPath;
		newPath.leftClip = readPaths[i].leftClip;
		newPath.rightClip = readPaths[i].rightClip;
		newPath.readPoses = readPaths[i].readPoses;
		newPath.readPosesExpandedStart = readPaths[i].readPosesExpandedStart;
		newPath.readPosesExpandedEnd = readPaths[i].readPosesExpandedEnd;
		newPath.readName = readPaths[i].readName;
		newPath.readLength = readPaths[i].readLength;
		std::vector<size_t> nodePosStarts;
		std::vector<size_t> nodePosEnds;
		size_t kmerPathLength = getNumberOfHashes(resolvableGraph, 0, 0, readPaths[i].path);
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (resolvables.count(readPaths[i].path[j].first) == 0 || unresolvables.count(readPaths[i].path[j].first) == 1)
			{
				newPath.path.push_back(readPaths[i].path[j]);
				size_t start = 0;
				if (j > 0)
				{
					start = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j });
					assert(start >= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j])));
					start -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
				}
				nodePosStarts.push_back(start);
				size_t end = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 });
				nodePosEnds.push_back(end);
				continue;
			}
			if (j > 0 && newEdgeNodes.count(std::make_pair(reverse(readPaths[i].path[j]), reverse(readPaths[i].path[j-1]))) == 1)
			{
				size_t start = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j });
				if (resolvables.count(readPaths[i].path[j-1].first) == 0 || unresolvables.count(readPaths[i].path[j-1].first) == 1)
				{
					assert(start >= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j])));
					start -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
				}
				else
				{
					assert(start >= resolvableGraph.unitigs[readPaths[i].path[j-1].first].size());
					start -= resolvableGraph.unitigs[readPaths[i].path[j-1].first].size();
				}
				nodePosStarts.push_back(start);
				size_t end = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 });
				nodePosEnds.push_back(end);
				newPath.path.emplace_back(newEdgeNodes.at(std::make_pair(reverse(readPaths[i].path[j]), reverse(readPaths[i].path[j-1]))), false);
			}
			if (j < readPaths[i].path.size()-1 && newEdgeNodes.count(std::make_pair(readPaths[i].path[j], readPaths[i].path[j+1])) == 1)
			{
				size_t start = 0;
				if (j > 0)
				{
					start = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j });
					assert(start >= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j])));
					start -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
				}
				nodePosStarts.push_back(start);
				size_t end = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 1 });
				if (resolvables.count(readPaths[i].path[j+1].first) == 1 && unresolvables.count(readPaths[i].path[j+1].first) == 0)
				{
					end = getNumberOfHashes(resolvableGraph, 0, 0, std::vector<std::pair<size_t, bool>> { readPaths[i].path.begin(), readPaths[i].path.begin() + j + 2 });
				}
				nodePosEnds.push_back(end);
				newPath.path.emplace_back(newEdgeNodes.at(std::make_pair(readPaths[i].path[j], readPaths[i].path[j+1])), true);
			}
		}
		if (newPath.path.size() == 0)
		{
			erasePath(resolvableGraph, readPaths, i);
			continue;
		}
		size_t startRemove = nodePosStarts[0];
		if (nodePosStarts[0] > newPath.leftClip)
		{
			newPath.readPoses.eraseRange(0, nodePosStarts[0] - newPath.leftClip);
			newPath.leftClip = 0;
		}
		else
		{
			assert(newPath.leftClip >= nodePosStarts[0]);
			newPath.leftClip -= nodePosStarts[0];
		}
		if (nodePosEnds.back() < kmerPathLength - newPath.rightClip)
		{
			size_t extraClip = (kmerPathLength - newPath.rightClip) - nodePosEnds.back();
			newPath.readPoses.eraseRange(newPath.readPoses.size() - extraClip, newPath.readPoses.size());
			newPath.rightClip = 0;
		}
		else
		{
			assert(newPath.rightClip >= (kmerPathLength - nodePosEnds.back()));
			newPath.rightClip -= (kmerPathLength - nodePosEnds.back());
		}
		assert(nodePosEnds.back() - nodePosStarts[0] == newPath.readPoses.size() + newPath.leftClip + newPath.rightClip);
		assert(nodePosStarts.size() == nodePosEnds.size());
		assert(nodePosEnds.size() == newPath.path.size());
		for (size_t j = 0; j < nodePosStarts.size(); j++)
		{
			assert(nodePosStarts[j] >= startRemove);
			nodePosStarts[j] -= startRemove;
			if (nodePosStarts[j] < newPath.leftClip)
			{
				nodePosStarts[j] = 0;
			}
			else
			{
				nodePosStarts[j] -= newPath.leftClip;
			}
			assert(nodePosEnds[j] >= startRemove);
			nodePosEnds[j] -= startRemove;
			assert(nodePosEnds[j] >= newPath.leftClip);
			nodePosEnds[j] -= newPath.leftClip;
			if (nodePosEnds[j] > newPath.readPoses.size()) nodePosEnds[j] = newPath.readPoses.size();
		}
		assert(nodePosStarts[0] == 0);
		assert(nodePosEnds.back() == newPath.readPoses.size());
		size_t lastStart = 0;
		for (size_t j = 1; j < newPath.path.size(); j++)
		{
			assert(resolvableGraph.edges[newPath.path[j-1]].count(newPath.path[j]) == resolvableGraph.edges[reverse(newPath.path[j])].count(reverse(newPath.path[j-1])));
			if (resolvableGraph.edges[newPath.path[j-1]].count(newPath.path[j]) == 0)
			{
				ReadPath path;
				path.leftClip = newPath.leftClip;
				if (lastStart > 0) path.leftClip = 0;
				path.rightClip = 0;
				path.path.insert(path.path.end(), newPath.path.begin() + lastStart, newPath.path.begin() + j);
				size_t posesStart = nodePosStarts[lastStart];
				size_t posesEnd = nodePosEnds[j-1];
				assert(posesStart < posesEnd);
				assert(posesEnd <= newPath.readPoses.size());
				for (size_t k = posesStart; k < posesEnd; k++)
				{
					path.readPoses.push_back(newPath.readPoses[k]);
					path.readPosesExpandedStart.push_back(newPath.readPosesExpandedStart[k]);
					path.readPosesExpandedEnd.push_back(newPath.readPosesExpandedEnd[k]);
				}
				path.readName = newPath.readName;
				path.readLength = newPath.readLength;
				addPath(resolvableGraph, readPaths, std::move(path));
				lastStart = j;
			}
		}
		ReadPath path;
		path.leftClip = newPath.leftClip;
		if (lastStart > 0) path.leftClip = 0;
		path.rightClip = newPath.rightClip;
		path.path.insert(path.path.end(), newPath.path.begin() + lastStart, newPath.path.end());
		size_t posesStart = nodePosStarts[lastStart];
		size_t posesEnd = nodePosEnds.back();
		assert(posesStart < posesEnd);
		assert(posesEnd <= newPath.readPoses.size());
		for (size_t k = posesStart; k < posesEnd; k++)
		{
			path.readPoses.push_back(newPath.readPoses[k]);
			path.readPosesExpandedStart.push_back(newPath.readPosesExpandedStart[k]);
			path.readPosesExpandedEnd.push_back(newPath.readPosesExpandedEnd[k]);
		}
		path.readName = newPath.readName;
		path.readLength = newPath.readLength;
		addPath(resolvableGraph, readPaths, std::move(path));
		erasePath(resolvableGraph, readPaths, i);
	}
}

struct ResolutionResult
{
public:
	size_t nodesResolved;
	size_t nodesAdded;
	std::unordered_set<size_t> maybeUnitigifiable;
};

ResolutionResult resolve(ResolvableUnitigGraph& resolvableGraph, const HashList& hashlist, const size_t kmerSize, std::vector<ReadPath>& readPaths, const std::unordered_set<size_t>& resolvables, const size_t minCoverage)
{
	ResolutionResult result;
	result.nodesResolved = 0;
	result.nodesAdded = 0;
	std::unordered_set<size_t> unresolvables;
	std::unordered_set<size_t> actuallyResolvables = resolvables;
	for (auto node : resolvables)
	{
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
		if (triplets.size() == 0)
		{
			unresolvables.insert(node);
			actuallyResolvables.erase(node);
			// unresolveRecursively(resolvableGraph, resolvables, unresolvables, node);
		}
	}
	std::vector<size_t> check;
	check.insert(check.end(), actuallyResolvables.begin(), actuallyResolvables.end());
	while (true)
	{
		std::vector<size_t> removeThese;
		std::unordered_set<size_t> newCheck;
		for (auto node : check)
		{
			if (unresolvables.count(node) == 1) continue;
			auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
			for (auto triplet : triplets)
			{
				if (resolvables.count(triplet.first.first) == 0 || unresolvables.count(triplet.first.first) == 1)
				{
					if (resolvables.count(triplet.second.first) == 0 || unresolvables.count(triplet.second.first) == 1)
					{
						if (resolvableGraph.unitigLength(triplet.first.first) == resolvableGraph.getBpOverlap(std::make_pair(node, false), reverse(triplet.first))+1)
						{
							if (resolvableGraph.unitigLength(triplet.second.first) == resolvableGraph.getBpOverlap(std::make_pair(node, true), triplet.second)+1)
							{
								assert(unresolvables.count(node) == 0);
								unresolvables.insert(node);
								removeThese.emplace_back(node);
								for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
								{
									if (actuallyResolvables.count(edge.first) == 1) newCheck.insert(edge.first);
								}
								for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
								{
									if (actuallyResolvables.count(edge.first) == 1) newCheck.insert(edge.first);
								}
								break;
							}
						}
					}
				}
			}
		}
		check.clear();
		check.insert(check.end(), newCheck.begin(), newCheck.end());
		for (auto node : removeThese) actuallyResolvables.erase(node);
		if (newCheck.size() == 0) break;
	}
	assert(resolvables.size() - unresolvables.size() == actuallyResolvables.size());
	if (unresolvables.size() == resolvables.size()) return result;
	std::unordered_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> newEdgeNodes;
	for (auto node : actuallyResolvables)
	{
		assert(getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage).size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto edge : resolvableGraph.edges[pos])
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1) continue;
			assert(actuallyResolvables.count(edge.first) == 1 || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, hashlist, kmerSize, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
		}
		pos = std::make_pair(node, false);
		for (auto edge : resolvableGraph.edges[pos])
		{
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 0);
			if ((resolvables.count(edge.first) == 0 || unresolvables.count(edge.first) == 1) && resolvableGraph.unitigLength(edge.first) == resolvableGraph.getBpOverlap(pos, edge) + 1) continue;
			assert(actuallyResolvables.count(edge.first) == 1 || resolvableGraph.unitigLength(edge.first) > resolvableGraph.getBpOverlap(pos, edge) + 1);
			// assert(unresolvables.count(edge.first) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(edge), reverse(pos))) == 1) continue;
			createEdgeNode(resolvableGraph, hashlist, kmerSize, newEdgeNodes, resolvables, unresolvables, pos, edge);
			assert(newEdgeNodes.count(std::make_pair(pos, edge)) == 1);
		}
	}
	for (auto node : actuallyResolvables)
	{
		auto triplets = getValidTriplets(resolvableGraph, resolvables, readPaths, node, minCoverage);
		assert(triplets.size() > 0);
		std::pair<size_t, bool> pos { node, true };
		for (auto triplet : triplets)
		{
			const std::pair<size_t, bool> before = triplet.first;
			const std::pair<size_t, bool> after = triplet.second;
			std::pair<size_t, bool> leftNode = before;
			std::pair<size_t, bool> rightNode = after;
			size_t overlap = resolvableGraph.unitigs[node].size();
			bool foundLeft = false;
			bool foundRight = false;
			assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 0 || newEdgeNodes.count(std::make_pair(before, pos)) == 0);
			assert(newEdgeNodes.count(std::make_pair(pos, after)) == 0 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
			if (newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(before, pos)) == 0);
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(pos), reverse(before))), false);
				foundLeft = true;
			}
			else if (newEdgeNodes.count(std::make_pair(before, pos)) == 1)
			{
				leftNode = std::make_pair(newEdgeNodes.at(std::make_pair(before, pos)), true);
				foundLeft = true;
			}
			else
			{
				assert(resolvables.count(before.first) == 0 || unresolvables.count(before.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(pos, after)) == 1 || newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1);
				assert(resolvableGraph.unitigLength(before.first) == resolvableGraph.getBpOverlap(before, pos)+1);
				overlap = resolvableGraph.overlaps.at(canon(before, pos));
			}
			if (newEdgeNodes.count(std::make_pair(pos, after)) == 1)
			{
				assert(newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 0);
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(pos, after)), true);
				foundRight = true;
			}
			else if (newEdgeNodes.count(std::make_pair(reverse(after), reverse(pos))) == 1)
			{
				rightNode = std::make_pair(newEdgeNodes.at(std::make_pair(reverse(after), reverse(pos))), false);
				foundRight = true;
			}
			else
			{
				assert(resolvables.count(after.first) == 0 || unresolvables.count(after.first) == 1);
				assert(newEdgeNodes.count(std::make_pair(reverse(pos), reverse(before))) == 1 || newEdgeNodes.count(std::make_pair(before, pos)) == 1);
				assert(resolvableGraph.unitigLength(after.first) == resolvableGraph.getBpOverlap(pos, after)+1);
				overlap = resolvableGraph.overlaps.at(canon(pos, after));
			}
			result.maybeUnitigifiable.insert(leftNode.first);
			result.maybeUnitigifiable.insert(rightNode.first);
			assert(foundLeft || foundRight);
			assert(!resolvableGraph.unitigRemoved[leftNode.first]);
			assert(!resolvableGraph.unitigRemoved[rightNode.first]);
			assert(resolvables.count(leftNode.first) == 0 || unresolvables.count(leftNode.first) == 1);
			assert(resolvables.count(rightNode.first) == 0 || unresolvables.count(rightNode.first) == 1);
			assert(resolvableGraph.edges[leftNode].count(rightNode) == 0);
			assert(resolvableGraph.edges[reverse(rightNode)].count(reverse(leftNode)) == 0);
			resolvableGraph.edges[leftNode].emplace(rightNode);
			resolvableGraph.edges[reverse(rightNode)].emplace(reverse(leftNode));
			resolvableGraph.overlaps[canon(leftNode, rightNode)] = overlap;
			auto leftSeq = resolvableGraph.unitigs[leftNode.first];
			if (!leftNode.second) leftSeq = revCompPath(leftSeq);
			auto rightSeq = resolvableGraph.unitigs[rightNode.first];
			if (!rightNode.second) rightSeq = revCompPath(rightSeq);
			assert(leftSeq.size() >= overlap);
			assert(rightSeq.size() >= overlap);
			for (size_t i = 0; i < overlap; i++)
			{
				assert(leftSeq[leftSeq.size()-overlap+i] == rightSeq[i]);
			}
			assert(resolvableGraph.edges[leftNode].count(rightNode) == 1);
			assert(resolvableGraph.edges[reverse(rightNode)].count(reverse(leftNode)) == 1);
			assert(resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(leftNode.first));
			assert(resolvableGraph.getBpOverlap(leftNode, rightNode) < resolvableGraph.unitigLength(rightNode.first));
		}
	}
	for (auto node : actuallyResolvables)
	{
		resolvableGraph.unitigRemoved[node] = true;
		std::pair<size_t, bool> fw { node, true };
		std::pair<size_t, bool> bw { node, false };
		for (auto edge : resolvableGraph.edges[fw])
		{
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(fw)) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(fw));
		}
		for (auto edge : resolvableGraph.edges[bw])
		{
			assert(resolvableGraph.edges[reverse(edge)].count(reverse(bw)) == 1);
			resolvableGraph.edges[reverse(edge)].erase(reverse(bw));
		}
		resolvableGraph.edges[fw].clear();
		resolvableGraph.edges[bw].clear();
	}
	replacePaths(resolvableGraph, readPaths, resolvables, unresolvables, newEdgeNodes);
	assert(resolvables.size() > unresolvables.size());
	result.nodesResolved = resolvables.size() - unresolvables.size();
	result.nodesAdded = newEdgeNodes.size();
	return result;
}

void checkValidity(const ResolvableUnitigGraph& graph, const std::vector<ReadPath>& readPaths, const size_t kmerSize)
{
	// return;
	assert(graph.unitigs.size() == graph.edges.size());
	assert(graph.unitigs.size() == graph.unitigRightClipBp.size());
	assert(graph.unitigs.size() == graph.unitigLeftClipBp.size());
	assert(graph.unitigs.size() == graph.unitigRemoved.size());
	assert(graph.unitigs.size() == graph.readsCrossingNode.size());
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		if (graph.unitigRemoved[i])
		{
			assert(graph.edges[std::make_pair(i, true)].size() == 0);
			assert(graph.edges[std::make_pair(i, false)].size() == 0);
			continue;
		}
		assert(graph.unitigLeftClipBp[i] < kmerSize);
		assert(graph.unitigRightClipBp[i] < kmerSize);
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges[fw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(fw)) == 1);
			assert(graph.edges[fw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2 || edge == fw);
		}
		for (auto edge : graph.edges[bw])
		{
			assert(!graph.unitigRemoved[edge.first]);
			assert(graph.edges[reverse(edge)].count(reverse(bw)) == 1);
			assert(graph.edges[bw].size() >= 2 || graph.edges[reverse(edge)].size() >= 2 || edge == bw);
		}
	}
	for (const auto& path : readPaths)
	{
		for (size_t i = 0; i < path.path.size(); i++)
		{
			assert(!graph.unitigRemoved[path.path[i].first]);
			if (i > 0) assert(graph.edges[path.path[i-1]].count(path.path[i]) == 1);
		}
	}
	std::vector<std::vector<size_t>> coverages;
	coverages.resize(graph.unitigs.size());
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		coverages[i].resize(graph.unitigs[i].size(), 0);
	}
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		if (path.path.size() == 1)
		{
			assert(path.leftClip + path.rightClip < graph.unitigs[path.path[0].first].size());
			for (size_t i = path.leftClip; i < graph.unitigs[path.path[0].first].size() - path.rightClip; i++)
			{
				size_t index = i;
				if (!path.path[0].second) index = graph.unitigs[path.path[0].first].size() - i - 1;
				coverages[path.path[0].first][index] += 1;
			}
			continue;
		}
		for (size_t i = 1; i < path.readPoses.size(); i++)
		{
			assert(path.readPoses[i] <= path.readPoses[i-1] + kmerSize);
		}
		assert(path.path.size() >= 2);
		assert(path.leftClip < graph.unitigs[path.path[0].first].size());
		for (size_t i = path.leftClip; i < graph.unitigs[path.path[0].first].size(); i++)
		{
			size_t index = i;
			if (!path.path[0].second) index = graph.unitigs[path.path[0].first].size() - i - 1;
			coverages[path.path[0].first][index] += 1;
		}
		for (size_t i = 1; i < path.path.size()-1; i++)
		{
			for (size_t j = 0; j < graph.unitigs[path.path[i].first].size(); j++)
			{
				coverages[path.path[i].first][j] += 1;
			}
		}
		assert(path.rightClip < graph.unitigs[path.path.back().first].size());
		for (size_t i = 0; i < graph.unitigs[path.path.back().first].size() - path.rightClip; i++)
		{
			size_t index = i;
			if (!path.path.back().second) index = graph.unitigs[path.path.back().first].size() - i - 1;
			coverages[path.path.back().first][index] += 1;
		}
	}
	for (size_t i = 0; i < coverages.size(); i++)
	{
		if (graph.unitigRemoved[i]) continue;
		for (size_t j = 0; j < coverages[i].size(); j++)
		{
			assert(coverages[i][j] > 0);
		}
	}
}

double getCoverage(const ResolvableUnitigGraph& resolvableGraph, const std::vector<ReadPath>& readPaths, const size_t unitig)
{
	double result = 0;
	for (auto i : resolvableGraph.readsCrossingNode[unitig])
	{
		if (readPaths[i].path.size() == 1)
		{
			assert(readPaths[i].path[0].first == unitig);
			assert(readPaths[i].leftClip + readPaths[i].rightClip < resolvableGraph.unitigs[readPaths[i].path[0].first].size());
			result += (double)(resolvableGraph.unitigs[readPaths[i].path[0].first].size() - readPaths[i].leftClip - readPaths[i].rightClip) / (double)(resolvableGraph.unitigs[readPaths[i].path[0].first].size());
		}
		else
		{
			for (size_t j = 0; j < readPaths[i].path.size(); j++)
			{
				if (readPaths[i].path[j].first != unitig) continue;
				assert(j != 0 || readPaths[i].leftClip < resolvableGraph.unitigs[unitig].size());
				assert(j != readPaths[i].path.size()-1 || readPaths[i].rightClip < resolvableGraph.unitigs[unitig].size());
				if (j == 0) result += (double)(resolvableGraph.unitigs[unitig].size() - readPaths[i].leftClip) / (double)(resolvableGraph.unitigs[unitig].size());
				if (j > 0 && j < readPaths[i].path.size()-1) result += 1;
				if (j == readPaths[i].path.size()-1) result += (double)(resolvableGraph.unitigs[unitig].size() - readPaths[i].rightClip) / (double)(resolvableGraph.unitigs[unitig].size());
			}
		}
	}
	return result;
}

size_t getEdgeCoverage(const ResolvableUnitigGraph& resolvableGraph, const std::vector<ReadPath>& readPaths, const std::pair<size_t, bool> from, const std::pair<size_t, bool> to)
{
	std::unordered_set<size_t> relevantReads;
	for (auto read : resolvableGraph.readsCrossingNode[from.first])
	{
		if (resolvableGraph.readsCrossingNode[to.first].count(read) == 1)
		{
			relevantReads.insert(read);
		}
	}
	size_t result = 0;
	for (size_t i : relevantReads)
	{
		for (size_t j = 1; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j-1] == from && readPaths[i].path[j] == to)
			{
				result += 1;
			}
			else if (readPaths[i].path[j-1] == reverse(to) && readPaths[i].path[j] == reverse(from))
			{
				result += 1;
			}
		}
	}
	return result;
}

void removeNode(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, size_t node)
{
	assert(!resolvableGraph.unitigRemoved[node]);
	for (auto edge : resolvableGraph.edges[std::make_pair(node, true)])
	{
		if (edge.first == node) continue;
		assert(resolvableGraph.edges[reverse(edge)].count(std::make_pair(node, false)) == 1);
		resolvableGraph.edges[reverse(edge)].erase(std::make_pair(node, false));
	}
	for (auto edge : resolvableGraph.edges[std::make_pair(node, false)])
	{
		if (edge.first == node) continue;
		assert(resolvableGraph.edges[reverse(edge)].count(std::make_pair(node, true)) == 1);
		resolvableGraph.edges[reverse(edge)].erase(std::make_pair(node, true));
	}
	resolvableGraph.edges[std::make_pair(node, true)].clear();
	resolvableGraph.edges[std::make_pair(node, false)].clear();
	std::unordered_set<size_t> relevantReads = resolvableGraph.readsCrossingNode[node];
	for (size_t i : relevantReads)
	{
		std::vector<size_t> nodePosStarts;
		std::vector<size_t> nodePosEnds;
		size_t pathKmerLength = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (j > 0) pathKmerLength -= resolvableGraph.overlaps.at(canon(readPaths[i].path[j-1], readPaths[i].path[j]));
			nodePosStarts.push_back(pathKmerLength);
			pathKmerLength += resolvableGraph.unitigs[readPaths[i].path[j].first].size();
			nodePosEnds.push_back(pathKmerLength);
		}
		assert(pathKmerLength == readPaths[i].readPoses.size() + readPaths[i].leftClip + readPaths[i].rightClip);
		assert(nodePosStarts[0] == 0);
		for (size_t j = 0; j < nodePosStarts.size(); j++)
		{
			if (nodePosStarts[j] < readPaths[i].leftClip)
			{
				nodePosStarts[j] = 0;
			}
			else
			{
				nodePosStarts[j] -= readPaths[i].leftClip;
			}
			assert(nodePosEnds[j] >= readPaths[i].leftClip);
			nodePosEnds[j] -= readPaths[i].leftClip;
			if (nodePosEnds[j] > readPaths[i].readPoses.size()) nodePosEnds[j] = readPaths[i].readPoses.size();
		}
		assert(nodePosStarts[0] == 0);
		assert(nodePosEnds.back() == readPaths[i].readPoses.size());
		size_t lastStart = 0;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (readPaths[i].path[j].first == node)
			{
				if (j == lastStart)
				{
					lastStart = j+1;
					continue;
				}
				ReadPath path;
				path.leftClip = readPaths[i].leftClip;
				if (lastStart > 0) path.leftClip = 0;
				path.rightClip = 0;
				path.path.insert(path.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.begin() + j);
				size_t posesStart = nodePosStarts[lastStart];
				size_t posesEnd = nodePosEnds[j-1];
				assert(posesStart < posesEnd);
				assert(posesEnd <= readPaths[i].readPoses.size());
				for (size_t k = posesStart; k < posesEnd; k++)
				{
					path.readPoses.push_back(readPaths[i].readPoses[k]);
					path.readPosesExpandedStart.push_back(readPaths[i].readPosesExpandedStart[k]);
					path.readPosesExpandedEnd.push_back(readPaths[i].readPosesExpandedEnd[k]);
				}
				path.readName = readPaths[i].readName;
				path.readLength = readPaths[i].readLength;
				addPath(resolvableGraph, readPaths, std::move(path));
				lastStart = j+1;
			}
		}
		if (lastStart < readPaths[i].path.size())
		{
			ReadPath path;
			path.leftClip = readPaths[i].leftClip;
			if (lastStart > 0) path.leftClip = 0;
			path.rightClip = readPaths[i].rightClip;
			path.path.insert(path.path.end(), readPaths[i].path.begin() + lastStart, readPaths[i].path.end());
			size_t posesStart = nodePosStarts[lastStart];
			size_t posesEnd = nodePosEnds.back();
			assert(posesStart < posesEnd);
			assert(posesEnd <= readPaths[i].readPoses.size());
			for (size_t k = posesStart; k < posesEnd; k++)
			{
				path.readPoses.push_back(readPaths[i].readPoses[k]);
				path.readPosesExpandedStart.push_back(readPaths[i].readPosesExpandedStart[k]);
				path.readPosesExpandedEnd.push_back(readPaths[i].readPosesExpandedEnd[k]);
			}
			path.readName = readPaths[i].readName;
			path.readLength = readPaths[i].readLength;
			addPath(resolvableGraph, readPaths, std::move(path));
		}
		erasePath(resolvableGraph, readPaths, i);
	}
	assert(resolvableGraph.readsCrossingNode[node].size() == 0);
	resolvableGraph.unitigRemoved[node] = true;
}

struct UntippingResult
{
	size_t nodesRemoved;
	std::unordered_set<size_t> maybeUnitigifiable;
};

void tryRemoveTip(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const size_t kmerSize, const size_t i, UntippingResult& result)
{
	assert(!resolvableGraph.unitigRemoved[i]);
	size_t maxEdgeCoverage = 0;
	size_t fwHasSafeEdge = false;
	size_t fwHasSafeNode = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)])
	{
		maxEdgeCoverage = std::max(maxEdgeCoverage, getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, true), edge));
		if (getCoverage(resolvableGraph, readPaths, edge.first) >= minSafeCoverage) fwHasSafeNode = true;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) fwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, true)].size() > 0 && (!fwHasSafeNode || !fwHasSafeEdge)) return;
	bool bwHasSafeNode = false;
	bool bwHasSafeEdge = false;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)])
	{
		maxEdgeCoverage = std::max(maxEdgeCoverage, getEdgeCoverage(resolvableGraph, readPaths, std::make_pair(i, false), edge));
		if (getCoverage(resolvableGraph, readPaths, edge.first) >= minSafeCoverage) bwHasSafeNode = true;
		for (auto edge2 : resolvableGraph.edges[reverse(edge)])
		{
			if (getEdgeCoverage(resolvableGraph, readPaths, reverse(edge), edge2) >= minSafeCoverage) bwHasSafeEdge = true;
		}
	}
	if (resolvableGraph.edges[std::make_pair(i, false)].size() > 0 && (!bwHasSafeNode || !bwHasSafeEdge)) return;
	if (maxEdgeCoverage > maxRemovableCoverage) return;
	for (auto edge : resolvableGraph.edges[std::make_pair(i, true)]) result.maybeUnitigifiable.insert(edge.first);
	for (auto edge : resolvableGraph.edges[std::make_pair(i, false)]) result.maybeUnitigifiable.insert(edge.first);
	removeNode(resolvableGraph, readPaths, i);
	result.nodesRemoved += 1;
}

UntippingResult removeLowCoverageTips(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const HashList& hashlist, const double maxRemovableCoverage, const double minSafeCoverage, const size_t maxRemovableLength, const size_t kmerSize, const std::unordered_set<size_t>& maybeUntippable)
{
	for (size_t i = resolvableGraph.lastTippableChecked; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = getCoverage(resolvableGraph, readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	for (auto i : maybeUntippable)
	{
		if (i >= resolvableGraph.lastTippableChecked) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, false)].size() >= 2) continue;
		if (resolvableGraph.edges[std::make_pair(i, true)].size() == 0 && resolvableGraph.edges[std::make_pair(i, false)].size() == 0) continue;
		size_t unitigLength = resolvableGraph.unitigLength(i);
		if (unitigLength > 10000) continue;
		double coverage = getCoverage(resolvableGraph, readPaths, i);
		if (coverage > 3) continue;
		resolvableGraph.everTippable.push_back(i);
	}
	resolvableGraph.lastTippableChecked = resolvableGraph.unitigs.size();
	UntippingResult result;
	result.nodesRemoved = 0;
	for (size_t index = resolvableGraph.everTippable.size()-1; index < resolvableGraph.everTippable.size(); index--)
	{
		size_t i = resolvableGraph.everTippable[index];
		if (resolvableGraph.unitigRemoved[i])
		{
			std::swap(resolvableGraph.everTippable[index], resolvableGraph.everTippable.back());
			resolvableGraph.everTippable.pop_back();
			continue;
		}
		tryRemoveTip(resolvableGraph, readPaths, hashlist, maxRemovableCoverage, minSafeCoverage, maxRemovableLength, kmerSize, i, result);
	}
	return result;
}

void resolveRound(ResolvableUnitigGraph& resolvableGraph, std::vector<ReadPath>& readPaths, const HashList& hashlist, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength)
{
	checkValidity(resolvableGraph, readPaths, kmerSize);
	std::priority_queue<size_t, std::vector<size_t>, UnitigLengthComparer> queue { UnitigLengthComparer { resolvableGraph } };
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		if (resolvableGraph.unitigRemoved[i]) continue;
		queue.emplace(i);
	}
	size_t lastTopSize = 0;
	size_t nodesRemoved = 0;
	while (queue.size() > 0)
	{
		size_t topSize = resolvableGraph.unitigLength(queue.top());
		if (topSize >= maxResolveLength) break;
		// assert(topSize >= lastTopSize);
		lastTopSize = topSize;
		std::unordered_set<size_t> resolvables;
		while (queue.size() > 0 && resolvableGraph.unitigLength(queue.top()) == topSize)
		{
			if (!resolvableGraph.unitigRemoved[queue.top()])
			{
				if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 1 && resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 1)
				{
					if (resolvableGraph.edges[std::make_pair(queue.top(), true)].size() >= 2 || resolvableGraph.edges[std::make_pair(queue.top(), false)].size() >= 2)
					{
						resolvables.emplace(queue.top());
					}
					else
					{
						assert(resolvableGraph.edges[std::make_pair(queue.top(), true)].size() == 1);
						assert(resolvableGraph.edges[std::make_pair(queue.top(), false)].size() == 1);
						if (resolvableGraph.unitigLength(resolvableGraph.edges[std::make_pair(queue.top(), true)].begin()->first) == topSize || resolvableGraph.unitigLength(resolvableGraph.edges[std::make_pair(queue.top(), false)].begin()->first) == topSize)
						{
							if (resolvableGraph.edges[std::make_pair(queue.top(), true)].begin()->first != queue.top())
							{
								resolvables.emplace(queue.top());
							}
						}
					}
				}
			}
			queue.pop();
		}
		if (resolvables.size() == 0) continue;
		assert(resolvables.size() > 0);
		size_t oldSize = resolvableGraph.unitigs.size();
		checkValidity(resolvableGraph, readPaths, kmerSize);
		std::cerr << "try resolve k=" << topSize;
		auto resolutionResult = resolve(resolvableGraph, hashlist, kmerSize, readPaths, resolvables, minCoverage);
		size_t newSize = resolvableGraph.unitigs.size();
		std::cerr << ", replaced " << resolutionResult.nodesResolved << " nodes with " << resolutionResult.nodesAdded << " nodes";
		nodesRemoved += resolutionResult.nodesResolved;
		if (resolutionResult.nodesResolved == 0)
		{
			std::cerr << std::endl;
			continue;
		}
		size_t unitigified = 0;
		size_t unitigifiedTo = 0;
		for (auto i : resolutionResult.maybeUnitigifiable)
		// for (size_t i = oldSize; i < newSize; i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			size_t unitigifiedHere = unitigifyOne(resolvableGraph, readPaths, i);
			nodesRemoved += unitigifiedHere;
			if (unitigifiedHere > 1)
			{
				unitigified += unitigifiedHere;
				unitigifiedTo += 1;
			}
			else
			{
				assert(!resolvableGraph.unitigRemoved[i]);
				queue.emplace(i);
			}
		}
		if (unitigified > 0)
		{
			std::cerr << ", unitigified " << unitigified << " nodes to " << unitigifiedTo << " nodes";
		}
		std::cerr << std::endl;
		auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, kmerSize, resolutionResult.maybeUnitigifiable);
		resolutionResult.maybeUnitigifiable.insert(removed.maybeUnitigifiable.begin(), removed.maybeUnitigifiable.end());
		auto removed2 = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 2, 5, 10000, kmerSize, resolutionResult.maybeUnitigifiable);
		nodesRemoved += removed.nodesRemoved + removed2.nodesRemoved;
		if (removed.nodesRemoved + removed2.nodesRemoved > 0)
		{
			std::cerr << "removed " << removed.nodesRemoved + removed2.nodesRemoved << " tips" << std::endl;
			removed.maybeUnitigifiable.insert(removed2.maybeUnitigifiable.begin(), removed2.maybeUnitigifiable.end());
			for (auto i : removed.maybeUnitigifiable)
			{
				if (resolvableGraph.unitigRemoved[i]) continue;
				size_t unitigifiedHere = unitigifyOne(resolvableGraph, readPaths, i);
				nodesRemoved += unitigifiedHere;
				if (unitigifiedHere == 0)
				{
					assert(!resolvableGraph.unitigRemoved[i]);
					queue.emplace(i);
				}
			}
		}
		checkValidity(resolvableGraph, readPaths, kmerSize);
		for (size_t i = oldSize; i < resolvableGraph.unitigs.size(); i++)
		{
			if (resolvableGraph.unitigRemoved[i]) continue;
			queue.emplace(i);
		}
		if (nodesRemoved > resolvableGraph.unitigs.size() / 2)
		{
			std::vector<size_t> queueNodes;
			while (queue.size() > 0)
			{
				auto top = queue.top();
				queue.pop();
				if (resolvableGraph.unitigRemoved[top]) continue;
				queueNodes.emplace_back(top);
			}
			size_t oldSize = resolvableGraph.unitigs.size();
			compact(resolvableGraph, readPaths, queueNodes);
			nodesRemoved = 0;
			size_t newSize = resolvableGraph.unitigs.size();
			std::cerr << "compacted from " << oldSize << " to " << newSize << std::endl;
			checkValidity(resolvableGraph, readPaths, kmerSize);
			for (auto node : queueNodes)
			{
				assert(node < resolvableGraph.unitigs.size());
				assert(!resolvableGraph.unitigRemoved[node]);
				queue.emplace(node);
			}
		}
	}
}

std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, const std::vector<std::string>& readFiles, const size_t numThreads, const ReadpartIterator& partIterator, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength)
{
	auto resolvableGraph = getUnitigs(initial, minCoverage, hashlist, kmerSize);
	auto readPaths = getUnitigPaths(resolvableGraph, hashlist, readFiles, numThreads, partIterator, kmerSize);
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (auto pos : readPaths[i].path)
		{
			resolvableGraph.readsCrossingNode[pos.first].emplace(i);
		}
		assert(readPaths[i].path.size() > 0);
		assert(getNumberOfHashes(resolvableGraph, readPaths[i].leftClip, readPaths[i].rightClip, readPaths[i].path) == readPaths[i].readPoses.size());
		if (readPaths[i].path.size() == 1)
		{
			assert(resolvableGraph.unitigs[readPaths[i].path[0].first].size() > readPaths[i].leftClip + readPaths[i].rightClip);
		}
		else
		{
			assert(resolvableGraph.unitigs[readPaths[i].path[0].first].size() > readPaths[i].leftClip);
			assert(resolvableGraph.unitigs[readPaths[i].path.back().first].size() > readPaths[i].rightClip);
		}
	}
	for (size_t i = 0; i < resolvableGraph.unitigs.size(); i++)
	{
		assert(resolvableGraph.readsCrossingNode[i].size() >= 1);
	}
	unitigifyAll(resolvableGraph, readPaths);
	auto removed = removeLowCoverageTips(resolvableGraph, readPaths, hashlist, 3, 10, 10000, kmerSize, std::unordered_set<size_t> {});
	if (removed.nodesRemoved > 0)
	{
		std::cerr << "removed " << removed.nodesRemoved << " tips" << std::endl;
		unitigifyAll(resolvableGraph, readPaths);
	}
	resolveRound(resolvableGraph, readPaths, hashlist, minCoverage, kmerSize, maxResolveLength);
	resolveRound(resolvableGraph, readPaths, hashlist, 1, kmerSize, maxResolveLength);
	checkValidity(resolvableGraph, readPaths, kmerSize);
	return resolvableToUnitigs(resolvableGraph, readPaths);
}
