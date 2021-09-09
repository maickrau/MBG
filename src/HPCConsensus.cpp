#include <limits>
#include "HPCConsensus.h"
#include "MBGCommon.h"
#include "VectorView.h"
#include "ConsensusMaker.h"

class KmerLocator
{
public:
	KmerLocator(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs)
	{
		nodes.emplace_back();
		nodes.back().parent = 0;
		nodes.back().depth = 0;
		std::get<0>(nodes.back().uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<1>(nodes.back().uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<2>(nodes.back().uniqueLocation) = true;
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			for (size_t j = 0; j < unitigs[i].size(); j++)
			{
				addPath(unitigs, i, j, true);
				addPath(unitigs, i, j, false);
			}
		}
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			for (size_t j = 0; j < unitigs[i].size(); j++)
			{
				auto found = locate(unitigs[i], j);
				if (std::get<0>(found) == std::numeric_limits<size_t>::max()) continue;
				assert(std::get<0>(found) == i);
				assert(std::get<1>(found) == j);
				assert(std::get<2>(found) == true);
			}
			auto rev = revCompPath(unitigs[i]);
			for (size_t j = 0; j < rev.size(); j++)
			{
				auto found = locate(rev, j);
				if (std::get<0>(found) == std::numeric_limits<size_t>::max()) continue;
				assert(std::get<0>(found) == i);
				assert(std::get<1>(found) == unitigs[i].size()-1-j);
				assert(std::get<2>(found) == false);
			}
		}
	}
	std::tuple<size_t, size_t, bool> locate(const std::vector<std::pair<size_t, bool>>& path, size_t start) const
	{
		std::pair<size_t, bool> pos = path[start];
		size_t node = 0;
		size_t depth = 0;
		assert(nodes[0].children.count(pos) == 1);
		while (true)
		{
			if (nodes[node].children.count(pos) == 0)
			{
				return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			}
			node = nodes[node].children.at(pos);
			if (std::get<0>(nodes[node].uniqueLocation) != std::numeric_limits<size_t>::max())
			{
				size_t i = std::get<0>(nodes[node].uniqueLocation);
				size_t j = std::get<1>(nodes[node].uniqueLocation);
				bool fw = std::get<2>(nodes[node].uniqueLocation);
				return std::make_tuple(i, j, fw);
			}
			depth += 1;
			if (start+depth >= path.size()) return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			pos = path[start+depth];
		}
		assert(false);
		return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
	}
private:
	void addPath(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs, size_t i, size_t j, bool fw)
	{
		size_t node = 0;
		size_t depth = 0;
		while (true)
		{
			size_t checkIndex = j;
			if (fw)
			{
				checkIndex += depth;
			}
			else
			{
				checkIndex -= depth;
			}
			if (checkIndex >= unitigs[i].size()) return;
			std::pair<size_t, bool> checkKmer = unitigs[i][checkIndex];
			if (!fw) checkKmer = reverse(checkKmer);
			if (nodes[node].children.count(checkKmer) == 0) break;
			node = nodes[node].children.at(checkKmer);
			depth += 1;
			if (std::get<0>(nodes[node].uniqueLocation) != std::numeric_limits<size_t>::max()) knockdown(unitigs, node);
		}
		size_t checkIndex = j;
		if (fw)
		{
			checkIndex += depth;
		}
		else
		{
			checkIndex -= depth;
		}
		assert(checkIndex < unitigs[i].size());
		size_t newNode = nodes.size();
		nodes.emplace_back();
		nodes[newNode].parent = node;
		nodes[newNode].depth = depth;
		nodes[newNode].uniqueLocation = std::make_tuple(i, j, fw);
		std::pair<size_t, bool> checkKmer = unitigs[i][checkIndex];
		if (!fw) checkKmer = reverse(checkKmer);
		nodes[node].children[checkKmer] = newNode;
	}
	void knockdown(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs, size_t node)
	{
		size_t i = std::get<0>(nodes[node].uniqueLocation);
		size_t j = std::get<1>(nodes[node].uniqueLocation);
		bool fw = std::get<2>(nodes[node].uniqueLocation);
		assert(i != std::numeric_limits<size_t>::max());
		size_t depth = nodes[node].depth;
		size_t next;
		if (fw)
		{
			next = j + depth + 1;
		}
		else
		{
			next = j - depth - 1;
		}
		std::get<0>(nodes[node].uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<1>(nodes[node].uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<2>(nodes[node].uniqueLocation) = true;
		assert(j < unitigs[i].size());
		if (next >= unitigs[i].size()) return;
		std::pair<size_t, bool> key;
		key = unitigs[i][next];
		if (!fw) key = reverse(key);
		assert(nodes[node].children.count(key) == 0);
		size_t nextIndex = nodes.size();
		nodes[node].children[key] = nextIndex;
		nodes.emplace_back();
		nodes[nextIndex].parent = node;
		nodes[nextIndex].depth = depth+1;
		std::get<0>(nodes[nextIndex].uniqueLocation) = i;
		std::get<1>(nodes[nextIndex].uniqueLocation) = j;
		std::get<2>(nodes[nextIndex].uniqueLocation) = fw;
	}
	class TreeNode
	{
	public:
		size_t parent;
		size_t depth;
		std::unordered_map<std::pair<size_t, bool>, size_t> children;
		std::tuple<size_t, size_t, bool> uniqueLocation;
	};
	std::vector<TreeNode> nodes;
};

void addCounts(ConsensusMaker& consensusMaker, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const size_t seqStart, const size_t seqEnd, const size_t unitig, const size_t unitigStart, const size_t unitigEnd, const bool fw)
{
	consensusMaker.addStrings(unitig, unitigStart, unitigEnd, [seqStart, seqEnd, &consensusMaker, &seq, &poses, &rawSeq, fw](size_t i)
	{
		size_t seqOff = seqStart + i;
		if (!fw) seqOff = seqEnd - 1 - i;
		assert(seqOff < seq.size());
		uint16_t compressed = seq[seqOff];
		if (!fw) compressed = complement(compressed);
		std::string seq;
		if (fw)
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			seq = rawSeq.substr(expandedStart, expandedEnd - expandedStart);
		}
		else
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			seq = revCompRaw(rawSeq.substr(expandedStart, expandedEnd - expandedStart));
		}
		return std::make_pair(compressed, seq);
	});
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> getHPCUnitigSequences(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const std::vector<ReadPath>& readPaths, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	ConsensusMaker consensusMaker;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>>> matchBlocks;
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		if (path.path.size() == 1)
		{
			size_t readStart = path.readPoses[0];
			size_t readEnd = path.readPoses.back() + kmerSize;
			size_t startKmerIndex = path.leftClip;
			size_t endKmerIndex = unitigs.unitigs[path.path[0].first].size() - 1 - path.rightClip;
			if (!path.path[0].second)
			{
				startKmerIndex = unitigs.unitigs[path.path[0].first].size() - 1 - startKmerIndex;
				endKmerIndex = unitigs.unitigs[path.path[0].first].size() - 1 - endKmerIndex;
				std::swap(startKmerIndex, endKmerIndex);
			}
			assert(endKmerIndex >= startKmerIndex);
			assert(endKmerIndex < unitigs.unitigs[path.path[0].first].size());
			size_t unitigStart = 0;
			size_t unitigEnd = 0;
			for (size_t i = 1; i <= startKmerIndex; i++)
			{
				unitigStart += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[0].first][i-1], unitigs.unitigs[path.path[0].first][i]);
			}
			unitigEnd = unitigStart;
			for (size_t i = startKmerIndex+1; i <= endKmerIndex; i++)
			{
				unitigEnd += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[0].first][i-1], unitigs.unitigs[path.path[0].first][i]);
			}
			unitigEnd += kmerSize;
			assert(readEnd > readStart);
			assert(unitigEnd > unitigStart);
			assert(readEnd - readStart == unitigEnd - unitigStart);
			matchBlocks[path.readName].emplace_back(path.path[0].first, readStart, unitigStart, readEnd - readStart, path.path[0].second);
			continue;
		}
		assert(path.path.size() >= 2);
		size_t startPosIndex = 0;
		size_t endPosIndex = unitigs.unitigs[path.path[0].first].size() - 1 - path.leftClip;
		size_t readStart = path.readPoses[startPosIndex];
		size_t readEnd = path.readPoses[endPosIndex] + kmerSize;
		size_t startKmerIndex = path.leftClip;
		size_t endKmerIndex = unitigs.unitigs[path.path[0].first].size()-1;
		if (!path.path[0].second)
		{
			startKmerIndex = unitigs.unitigs[path.path[0].first].size() - 1 - startKmerIndex;
			endKmerIndex = unitigs.unitigs[path.path[0].first].size() - 1 - endKmerIndex;
			std::swap(startKmerIndex, endKmerIndex);
		}
		assert(endKmerIndex >= startKmerIndex);
		assert(endKmerIndex < unitigs.unitigs[path.path[0].first].size());
		size_t unitigStart = 0;
		size_t unitigEnd = 0;
		for (size_t i = 1; i <= startKmerIndex; i++)
		{
			unitigStart += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[0].first][i-1], unitigs.unitigs[path.path[0].first][i]);
		}
		unitigEnd = unitigStart;
		for (size_t i = startKmerIndex+1; i <= endKmerIndex; i++)
		{
			unitigEnd += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[0].first][i-1], unitigs.unitigs[path.path[0].first][i]);
		}
		unitigEnd += kmerSize;
		assert(readEnd > readStart);
		assert(unitigEnd > unitigStart);
		assert(readEnd - readStart == unitigEnd - unitigStart);
		matchBlocks[path.readName].emplace_back(path.path[0].first, readStart, unitigStart, readEnd - readStart, path.path[0].second);
		for (size_t i = 1; i < path.path.size()-1; i++)
		{
			startPosIndex = endPosIndex + 1 - unitigs.edgeOverlap(path.path[i-1], path.path[i]);
			endPosIndex = startPosIndex + unitigs.unitigs[path.path[i].first].size() - 1;
			readStart = path.readPoses[startPosIndex];
			readEnd = path.readPoses[endPosIndex] + kmerSize;
			startKmerIndex = 0;
			endKmerIndex = unitigs.unitigs[path.path[i].first].size()-1;
			unitigStart = 0;
			unitigEnd = 0;
			for (size_t j = 1; j <= startKmerIndex; j++)
			{
				unitigStart += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[i].first][j-1], unitigs.unitigs[path.path[i].first][j]);
			}
			unitigEnd = unitigStart;
			for (size_t j = startKmerIndex+1; j <= endKmerIndex; j++)
			{
				unitigEnd += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path[i].first][j-1], unitigs.unitigs[path.path[i].first][j]);
			}
			unitigEnd += kmerSize;
			assert(readEnd > readStart);
			assert(unitigEnd > unitigStart);
			assert(readEnd - readStart == unitigEnd - unitigStart);
			matchBlocks[path.readName].emplace_back(path.path[i].first, readStart, unitigStart, readEnd - readStart, path.path[i].second);
		}
		startPosIndex = endPosIndex + 1 - unitigs.edgeOverlap(path.path[path.path.size()-2], path.path.back());
		endPosIndex = startPosIndex + unitigs.unitigs[path.path.back().first].size() - 1 - path.rightClip;
		assert(endPosIndex == path.readPoses.size()-1);
		readStart = path.readPoses[startPosIndex];
		readEnd = path.readPoses[endPosIndex] + kmerSize;
		startKmerIndex = 0;
		endKmerIndex = unitigs.unitigs[path.path.back().first].size()-1 - path.rightClip;
		if (!path.path.back().second)
		{
			startKmerIndex = unitigs.unitigs[path.path.back().first].size() - 1 - startKmerIndex;
			endKmerIndex = unitigs.unitigs[path.path.back().first].size() - 1 - endKmerIndex;
			std::swap(startKmerIndex, endKmerIndex);
		}
		assert(endKmerIndex >= startKmerIndex);
		assert(endKmerIndex < unitigs.unitigs[path.path.back().first].size());
		unitigStart = 0;
		unitigEnd = 0;
		for (size_t j = 1; j <= startKmerIndex; j++)
		{
			unitigStart += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path.back().first][j-1], unitigs.unitigs[path.path.back().first][j]);
		}
		unitigEnd = unitigStart;
		for (size_t j = startKmerIndex+1; j <= endKmerIndex; j++)
		{
			unitigEnd += kmerSize - hashlist.getOverlap(unitigs.unitigs[path.path.back().first][j-1], unitigs.unitigs[path.path.back().first][j]);
		}
		unitigEnd += kmerSize;
		assert(readEnd > readStart);
		assert(unitigEnd > unitigStart);
		assert(readEnd - readStart == unitigEnd - unitigStart);
		matchBlocks[path.readName].emplace_back(path.path.back().first, readStart, unitigStart, readEnd - readStart, path.path.back().second);
	}
	std::vector<size_t> unitigLengths;
	std::vector<std::vector<size_t>> bpOffsets;
	bpOffsets.resize(unitigs.unitigs.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t offset = 0;
		bpOffsets[i].push_back(0);
		for (size_t j = 1; j < unitigs.unitigs[i].size(); j++)
		{
			size_t rleOverlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			assert(rleOverlap < kmerSize);
			offset += kmerSize - rleOverlap;
			bpOffsets[i].push_back(offset);
		}
		unitigLengths.push_back(offset + kmerSize);
	}
	consensusMaker.init(unitigLengths);
	iterateReadsMultithreaded(filenames, numThreads, [&consensusMaker, &unitigLengths, &bpOffsets, &unitigs, &partIterator, &hashlist, &matchBlocks, kmerSize](size_t thread, FastQ& read)
	{
		if (matchBlocks.count(read.seq_id) == 0) return;
		partIterator.iterateParts(read, [read, &consensusMaker, &hashlist, &unitigLengths, &unitigs, &bpOffsets, &matchBlocks, kmerSize](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq)
		{
			for (auto block : matchBlocks.at(read.seq_id))
			{
				size_t unitig = std::get<0>(block);
				size_t readStartPos = std::get<1>(block);
				size_t unitigStartPos = std::get<2>(block);
				size_t matchLength = std::get<3>(block);
				bool forward = std::get<4>(block);
				addCounts(consensusMaker, seq, poses, rawSeq, readStartPos, readStartPos + matchLength, unitig, unitigStartPos, unitigStartPos + matchLength, forward);
			}
		});
	});
	return consensusMaker.getSequences();
}
