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
		unitigLengths.push_back(offset + kmerSize);
	}
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
			size_t unitigStart = bpOffsets[path.path[0].first][startKmerIndex];
			size_t unitigEnd = bpOffsets[path.path[0].first][endKmerIndex] + kmerSize;
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
		size_t unitigStart = bpOffsets[path.path[0].first][startKmerIndex];
		size_t unitigEnd = bpOffsets[path.path[0].first][endKmerIndex] + kmerSize;
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
			unitigStart = bpOffsets[path.path[i].first][startKmerIndex];
			unitigEnd = bpOffsets[path.path[i].first][endKmerIndex] + kmerSize;
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
		unitigStart = bpOffsets[path.path.back().first][startKmerIndex];
		unitigEnd = bpOffsets[path.path.back().first][endKmerIndex] + kmerSize;
		assert(readEnd > readStart);
		assert(unitigEnd > unitigStart);
		assert(readEnd - readStart == unitigEnd - unitigStart);
		matchBlocks[path.readName].emplace_back(path.path.back().first, readStart, unitigStart, readEnd - readStart, path.path.back().second);
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
