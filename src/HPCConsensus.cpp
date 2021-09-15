#include <limits>
#include "HPCConsensus.h"
#include "MBGCommon.h"
#include "VectorView.h"
#include "ConsensusMaker.h"

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
		assert(offset + kmerSize > unitigs.leftClip[i] + unitigs.rightClip[i]);
		unitigLengths.push_back(offset + kmerSize - unitigs.leftClip[i] - unitigs.rightClip[i]);
	}
	ConsensusMaker consensusMaker;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>>> matchBlocks;
	for (const auto& path : readPaths)
	{
		if (path.path.size() == 0) continue;
		size_t pathKmerCount = 0;
		std::vector<size_t> pathKmerStarts;
		std::vector<size_t> pathKmerEnds;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			if (i > 0) pathKmerCount -= unitigs.edgeOverlap(path.path[i-1], path.path[i]);
			pathKmerStarts.push_back(pathKmerCount);
			pathKmerCount += unitigs.unitigs[path.path[i].first].size();
			pathKmerEnds.push_back(pathKmerCount);
		}
		assert(pathKmerStarts[0] == 0);
		assert(pathKmerEnds.back() == pathKmerCount);
		assert(pathKmerCount == path.readPoses.size() + path.leftClip + path.rightClip);
		for (size_t i = 0; i < path.path.size(); i++)
		{
			size_t wantedStart = pathKmerStarts[i];
			size_t wantedEnd = pathKmerEnds[i];
			size_t possibleStart = path.leftClip;
			size_t possibleEnd = pathKmerCount - path.rightClip;
			size_t startKmerIndex = 0;
			size_t endKmerIndex = 0;
			size_t startPosIndex = 0;
			size_t endPosIndex = 0;
			if (wantedStart < possibleStart)
			{
				startPosIndex = 0;
				startKmerIndex = possibleStart - wantedStart;
			}
			else
			{
				startPosIndex = wantedStart - possibleStart;
				startKmerIndex = 0;
			}
			assert(wantedEnd > wantedStart + startKmerIndex);
			assert(possibleEnd > possibleStart + startPosIndex);
			size_t maxMatch = std::min(wantedEnd - wantedStart - startKmerIndex, possibleEnd - possibleStart - startPosIndex);
			endKmerIndex = startKmerIndex + maxMatch - 1;
			endPosIndex = startPosIndex + maxMatch - 1;
			if (!path.path[i].second)
			{
				startKmerIndex = bpOffsets[path.path[i].first].size() - 1 - startKmerIndex;
				endKmerIndex = bpOffsets[path.path[i].first].size() - 1 - endKmerIndex;
				std::swap(startKmerIndex, endKmerIndex);
			}
			assert(i != 0 || startPosIndex == 0);
			assert(i != path.path.size()-1 || endPosIndex == path.readPoses.size()-1);
			assert(startKmerIndex <= endKmerIndex);
			assert(startPosIndex <= endPosIndex);
			assert(endPosIndex < path.readPoses.size());
			assert(endKmerIndex < bpOffsets[path.path[i].first].size());
			size_t readStart = path.readPoses[startPosIndex];
			size_t readEnd = path.readPoses[endPosIndex] + kmerSize;
			size_t unitigStart = bpOffsets[path.path[i].first][startKmerIndex];
			size_t unitigEnd = bpOffsets[path.path[i].first][endKmerIndex] + kmerSize;
			size_t unitig = path.path[i].first;
			bool fw = path.path[i].second;
			if (unitigStart < unitigs.leftClip[unitig])
			{
				if (fw)
				{
					readStart += unitigs.leftClip[unitig] - unitigStart;
				}
				else
				{
					readEnd -= unitigs.leftClip[unitig] - unitigStart;
				}
				unitigStart = 0;
			}
			else
			{
				unitigStart -= unitigs.leftClip[unitig];
			}
			assert(unitigEnd >= unitigs.leftClip[unitig]);
			unitigEnd -= unitigs.leftClip[unitig];
			if (unitigEnd > unitigLengths[unitig])
			{
				if (fw)
				{
					readEnd -= unitigEnd - unitigLengths[unitig];
				}
				else
				{
					readStart += unitigEnd - unitigLengths[unitig];
				}
				unitigEnd = unitigLengths[unitig];
			}
			assert(readEnd > readStart);
			assert(unitigEnd > unitigStart);
			assert(readEnd - readStart == unitigEnd - unitigStart);
			matchBlocks[path.readName].emplace_back(unitig, readStart, unitigStart, readEnd - readStart, fw);
		}
	}
	consensusMaker.init(unitigLengths);
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::pair<size_t, bool> pos { i, true };
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[edge.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.addEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
		pos = std::make_pair(i, false);
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[pos.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.addEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
	}
	consensusMaker.findParentLinks();
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
