#ifndef KmerMatcher_h
#define KmerMatcher_h

#include <vector>
#include <tuple>
#include "HashList.h"
#include "UnitigGraph.h"
#include "UnitigResolver.h"

std::vector<std::tuple<size_t, size_t, bool>> getKmerLocator(const UnitigGraph& graph);

template <typename F>
void iterateReadPaths(const UnitigGraph& graph, const HashList& hashlist, const size_t kmerSize, const std::vector<std::tuple<size_t, size_t, bool>>& kmerLocator, const ReadInfo& read, const std::vector<size_t>& positions, const std::vector<HashType>& hashes, F callback)
{
	ReadPath current;
	current.readName = read.readName;
	current.readLength = read.readLength;
	current.readLengthHPC = read.readLengthHpc;
	size_t lastReadPos = std::numeric_limits<size_t>::max();
	std::pair<size_t, bool> lastKmer { std::numeric_limits<size_t>::max(), true };
	std::tuple<size_t, size_t, bool> lastPos = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
	assert(positions.size() == hashes.size());
	for (size_t i = 0; i < positions.size(); i++)
	{
		const size_t readPos = positions[i];
		const HashType fwHash = hashes[i];
		assert(readPos + kmerSize <= read.readLengthHpc);
		std::pair<size_t, bool> kmer = hashlist.getNodeOrNull(fwHash);
		if (kmer.first == std::numeric_limits<size_t>::max()) continue;
		if (kmer.first >= kmerLocator.size() || std::get<0>(kmerLocator[kmer.first]) == std::numeric_limits<size_t>::max())
		{
			if (current.path.size() > 0)
			{
				assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
				callback(current);
			}
			current.path.clear();
			current.readName = read.readName;
			current.readLength = read.readLength;
			current.readLengthHPC = read.readLengthHpc;
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
					callback(current);
				}
				current.path.clear();
				current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
				current.readPoses.clear();
				current.readPoses.push_back(readPos);
				current.leftClip = std::get<1>(pos);
				current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
				current.readName = read.readName;
				current.readLength = read.readLength;
				current.readLengthHPC = read.readLengthHpc;
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
		assert(graph.edges.hasEdge(fromEdge, toEdge) == graph.edges.hasEdge(reverse(toEdge), reverse(fromEdge)));
		if (graph.edges.hasEdge(fromEdge, toEdge) == 1 && std::get<1>(pos) == 0 && std::get<1>(lastPos) == graph.unitigs[std::get<0>(lastPos)].size()-1)
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
			callback(current);
		}
		current.path.clear();
		current.path.emplace_back(std::get<0>(pos), std::get<2>(pos));
		current.readPoses.clear();
		current.readPoses.push_back(readPos);
		current.leftClip = std::get<1>(pos);
		current.rightClip = graph.unitigs[std::get<0>(pos)].size() - 1 - std::get<1>(pos);
		current.readName = read.readName;
		current.readLength = read.readLength;
		current.readLengthHPC = read.readLengthHpc;
		assert(current.leftClip + current.rightClip + 1 == graph.unitigs[std::get<0>(pos)].size());
		lastPos = pos;
		lastReadPos = readPos;
		lastKmer = kmer;
		continue;
	}
	if (current.path.size() > 0)
	{
		assert(current.readPoses.back() + kmerSize <= current.readLengthHPC);
		callback(current);
	}
}

#endif
